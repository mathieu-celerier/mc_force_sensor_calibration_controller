#include "calibrate.h"

#include <ceres/ceres.h>

#include <mc_rtc/io_utils.h>
#include <RBDyn/FK.h>
#include <SpaceVecAlg/SpaceVecAlg>

#include "CostFunctor.h"

struct Minimize
{
  template<typename T>
  bool operator()(const T * const x, T * residual) const
  {
    residual[0] = x[0];
    residual[1] = x[1];
    residual[2] = x[2];
    return true;
  }
};

inline std::vector<std::string> getSuccessorBodies(const mc_rbdyn::Robot & robot_,
                                                   const std::string & rootBody,
                                                   bool includeRoot = false)
{
  auto & robot = const_cast<mc_rbdyn::Robot &>(robot_); // For successorJoints,
                                                        // should be const
  auto bIdx = robot.bodyIndexByName(rootBody);
  // Graph of successor joint built using robot's root as the root of the graph
  // This returns a map of
  // - BODY1 -> [SUCCESSOR JOINT1 of BODY1, SUCCESSOR JOINT2 of BODY1.. ]
  // - BODYN -> [SUCCESSOR JOINT1 of BODYN, SUCCESSOR JOINT2 of BODYN.. ]
  auto successorJointsGraph = robot.mbg().successorJoints(robot.mb().body(0).name());

  std::function<std::vector<std::string>(const std::vector<std::string> & succJoints)> computeSuccBodyNames;

  std::vector<std::string> successorBodyNames;
  if(includeRoot)
  {
    successorBodyNames.push_back(rootBody);
  }
  computeSuccBodyNames = [&successorBodyNames, &successorJointsGraph, &robot,
                          &computeSuccBodyNames](const std::vector<std::string> & succJoints)
  {
    for(const auto & joint : succJoints)
    {
      auto successorBodyIdx = robot.mb().successor(robot.mb().jointIndexByName(joint));
      const auto & successorBodyName = robot.mb().body(successorBodyIdx).name();
      const auto & successorJoints = successorJointsGraph.at(successorBodyName);
      // Name of the successor body of this joint + name of all its successors
      successorBodyNames.push_back(successorBodyName);
      computeSuccBodyNames(successorJoints);
    }
    return successorBodyNames;
  };

  return computeSuccBodyNames(successorJointsGraph[rootBody]);
}

CalibrationResult calibrate(const mc_rbdyn::Robot & robot,
                            const std::string & sensorN,
                            const Measurements & measurements,
                            const InitialGuess & initialGuess,
                            const Eigen::Matrix3d & I_com,
                            bool verbose)
{
  CalibrationResult result{initialGuess};
  const auto & sensor = robot.forceSensor(sensorN);

  // Build the problem.
  ceres::Problem problem;

  auto & mass = result.mass;
  auto * com = result.com.data();
  auto * rpy = result.rpy.data();
  auto * offset = result.offset.data();

  // For each measurement add a residual block
  for(const auto & measurement : measurements)
  {
    const auto & pos = measurement.X_0_p;
    const auto & f = measurement.measure;
    ceres::CostFunction * cost_function = new ceres::AutoDiffCostFunction<CostFunctor, 6, 1, 3, 3, 6>(
        new CostFunctor(pos, f, sensor.X_p_f(), measurement.V_p, measurement.A_p, I_com));
    problem.AddResidualBlock(cost_function, new ceres::CauchyLoss(0.5), &mass, rpy, com, offset);
  }
  ceres::CostFunction * min_rpy = new ceres::AutoDiffCostFunction<Minimize, 3, 3>(new Minimize());
  problem.AddResidualBlock(min_rpy, nullptr, rpy);
  problem.SetParameterLowerBound(&mass, 0, 0);
  problem.SetParameterLowerBound(rpy, 0, -2 * M_PI);
  problem.SetParameterLowerBound(rpy, 1, -2 * M_PI);
  problem.SetParameterLowerBound(rpy, 2, -2 * M_PI);
  problem.SetParameterUpperBound(rpy, 0, 2 * M_PI);
  problem.SetParameterUpperBound(rpy, 1, 2 * M_PI);
  problem.SetParameterUpperBound(rpy, 2, 2 * M_PI);

  // Run the solver!
  ceres::Solver::Options options;
  options.linear_solver_type = ceres::DENSE_QR;
  options.minimizer_progress_to_stdout = verbose;
  options.max_num_iterations = 1000;
  ceres::Solver::Summary summary;
  Solve(options, &problem, &summary);
  result.success = summary.IsSolutionUsable();

  // clang-format off
  mc_rtc::log::info(
R"({},
success     : {}
mass        : {}
rpy         : {}, {}, {}
com         : {}, {}, {}
force offset: {}, {}, {}, {}, {}, {})",
      summary.BriefReport(),
      result.success,
      result.mass,
      result.rpy[0], result.rpy[1], result.rpy[2],
      result.com[0], result.com[1], result.com[2],
      result.offset[0], result.offset[1], result.offset[2], result.offset[3], result.offset[4], result.offset[5]);
  // clang-format on
  return result;
}

InitialGuess computeInitialGuessFromModel(const mc_rbdyn::Robot & robot,
                                          const std::string & sensorN,
                                          bool includeParent,
                                          bool verbose)
{
  InitialGuess guess;
  const auto & sensor = robot.forceSensor(sensorN);
  const auto & X_0_parent = robot.bodyPosW(sensor.parentBody());

  // Compute initial guess from the robot model:
  // - mass: mass of all links under the force sensor. This assumes:
  //   1/ That the force sensor is attached to its parent link in the model (and
  //   thus that its mass is accounted for in the parent link)
  //   2/ That the mass of all child links under the force sensor is correct
  // FIXME for now include the parent body in the computation as HRP2 model is
  // buggy. Remove "true" argument once fixed
  auto successorBodies = getSuccessorBodies(robot, sensor.parentBody(), includeParent);
  double totalMass = 0;
  Eigen::Vector3d com = Eigen::Vector3d::Zero();
  if(successorBodies.size())
  {
    auto & mb = robot.mb();
    auto & mbc = robot.mbc();
    for(const auto & bodyName : successorBodies)
    {
      auto bodyIndex = robot.mb().bodyIndexByName(bodyName);
      const auto & body = robot.mb().body(bodyIndex);
      double mass = body.inertia().mass();
      totalMass += mass;
      auto X_parent_body = mbc.bodyPosW[bodyIndex] * X_0_parent.inv();
      sva::PTransformd scaledBobyPosW(X_parent_body.rotation(), mass * X_parent_body.translation());
      com += (sva::PTransformd(body.inertia().momentum()) * scaledBobyPosW).translation();
      if(verbose)
      {
        mc_rtc::log::info("[Initial Guess] Body: {}, Mass: {}", body.name(), mass);
      }
    }
    if(totalMass > 0)
    {
      com /= totalMass;
    }
    else
    {
      mc_rtc::log::warning(
          "Warning: no mass provided for the bodies attached to force sensor \"{}\", assuming mass = 0 and CoM=[0,0,0]",
          sensor.name());
    }
  }

  guess.mass = totalMass;
  guess.com[0] = com[0];
  guess.com[1] = com[1];
  guess.com[2] = com[2];

  return guess;
}

Eigen::Matrix3d computeToolInertia(const mc_rbdyn::Robot & robot, const std::string & sensorN, bool includeParent)
{
  const auto & sensor = robot.forceSensor(sensorN);
  const auto & X_0_parent = robot.bodyPosW(sensor.parentBody());

  auto successorBodies = getSuccessorBodies(robot, sensor.parentBody(), includeParent);
  sva::RBInertiad combined(0, Eigen::Vector3d::Zero(), Eigen::Matrix3d::Zero());
  auto & mbc = robot.mbc();
  for(const auto & bodyName : successorBodies)
  {
    auto bodyIndex = robot.mb().bodyIndexByName(bodyName);
    const auto & body = robot.mb().body(bodyIndex);
    // Same X_parent_body convention as computeInitialGuessFromModel, used here to bring each
    // body's own spatial (rigid body) inertia into the sensor parent body's frame so they can
    // be summed directly.
    auto X_parent_body = mbc.bodyPosW[bodyIndex] * X_0_parent.inv();
    combined = combined + X_parent_body.transMul(body.inertia());
  }

  if(combined.mass() <= 0)
  {
    return Eigen::Matrix3d::Zero();
  }
  // Re-express about the combined CoM (rather than the parent body's origin) since that is
  // what the calibration's rotational inertial term needs.
  sva::PTransformd X_p_vb(Eigen::Vector3d(combined.momentum() / combined.mass()));
  sva::RBInertiad aboutCoM = X_p_vb.dualMul(combined);
  return aboutCoM.inertia();
}
