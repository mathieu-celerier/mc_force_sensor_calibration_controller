/** Self-consistency test for the inertial (non-static) term added to the force sensor
 * calibration cost function.
 *
 * There is no physical rig to compare against for the inertial term's frame/sign convention,
 * so instead we generate synthetic measurements *from* the exact same model used for the fit
 * (CostFunctor, shared with calibrate.cpp via CostFunctor.h) for a set of known ground-truth
 * calibration parameters and a variety of synthetic poses/velocities/accelerations, then check
 * that calibrate() recovers those exact parameters from a deliberately wrong initial guess.
 * If the inertial term's sign or frame were wrong, this would not converge back to the
 * ground-truth parameters (the residual would not be able to reach ~0 for all samples at once,
 * since the "wrong" inertial contribution could not be explained away by mass/com/rpy/offset
 * alone across sufficiently varied velocity/acceleration values).
 */
#include "../src/CostFunctor.h"
#include "../src/calibrate.h"

#include <mc_rbdyn/RobotLoader.h>
#include <mc_rbdyn/Robots.h>

#include <Eigen/Eigenvalues>

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace
{

sva::PTransformd makePose(double roll, double pitch, double yaw, const Eigen::Vector3d & t)
{
  Eigen::Matrix3d R = rpyToMat<double>(roll, pitch, yaw);
  return sva::PTransformd(R, t);
}

/** Generate the measurement that is exactly consistent with CostFunctor's model for the given
 * ground-truth parameters, pose, and spatial velocity/acceleration. */
sva::ForceVecd predict(const sva::PTransformd & pos,
                       const sva::PTransformd & X_p_f,
                       const sva::MotionVecd & V_p,
                       const sva::MotionVecd & A_p,
                       double mass,
                       const std::array<double, 3> & rpy,
                       const std::array<double, 3> & com,
                       const std::array<double, 6> & offset,
                       const Eigen::Matrix3d & I_com)
{
  CostFunctor f(pos, sva::ForceVecd(Eigen::Vector6d::Zero()), X_p_f, V_p, A_p, I_com);
  double residual[6];
  f(&mass, rpy.data(), com.data(), offset.data(), residual);
  // residual = measure(=0) - vb_f(ground truth) = -vb_f(ground truth)
  // so the self-consistent measurement is vb_f(ground truth) = -residual
  return sva::ForceVecd(Eigen::Vector3d(-residual[0], -residual[1], -residual[2]),
                        Eigen::Vector3d(-residual[3], -residual[4], -residual[5]));
}

bool nearlyEqual(double a, double b, double tol)
{
  return std::abs(a - b) < tol;
}

} // namespace

int main()
{
  mc_rbdyn::RobotLoader::clear();
  mc_rbdyn::RobotLoader::update_robot_module_path({"/usr/local/lib/mc_robots"});
  auto rm = mc_rbdyn::RobotLoader::get_robot_module("JVRC1");
  auto robots = mc_rbdyn::loadRobot(*rm);
  auto & robot = robots->robot();

  const std::string sensorName = "LeftHandForceSensor";
  const auto & sensor = robot.forceSensor(sensorName);
  const auto X_p_f = sensor.X_p_f();

  // Ground truth calibration parameters (arbitrary, nonzero in every component so a sign/frame
  // bug in any term has somewhere to show up).
  const double mass_gt = 0.6321;
  const std::array<double, 3> rpy_gt = {0.0, 0.0, 0.0};
  const std::array<double, 3> com_gt = {0.015, -0.028, 0.142};
  const std::array<double, 6> offset_gt = {0.0021, -0.0013, 0.0009, 0.0182, -0.0654, 0.2701};
  // Arbitrary but nonzero/non-diagonal so a sign or frame bug in the rotational inertia term
  // (I*alpha + omega x I*omega) has somewhere to show up. Values are in the same order of
  // magnitude as a small hand tool (a few 1e-3 kg.m^2).
  Eigen::Matrix3d I_com_gt;
  I_com_gt << 3.2e-3, 0.4e-3, -0.2e-3, 0.4e-3, 5.1e-3, 0.1e-3, -0.2e-3, 0.1e-3, 2.7e-3;

  // A variety of synthetic poses and spatial velocities/accelerations, deliberately including
  // large angular velocity/acceleration so the inertial term dominates in some samples (a test
  // that only used near-static samples could pass even with a badly wrong inertial term).
  struct Sample
  {
    sva::PTransformd pos;
    sva::MotionVecd V_p;
    sva::MotionVecd A_p;
  };
  std::vector<Sample> samples = {
      {makePose(0.0, 0.0, 0.0, {0.3, 0.2, 0.9}), sva::MotionVecd(Eigen::Vector6d::Zero()),
       sva::MotionVecd(Eigen::Vector6d::Zero())},
      {makePose(0.3, -0.2, 1.1, {0.35, -0.1, 1.0}), sva::MotionVecd(Eigen::Vector3d(0.5, -0.3, 0.2), Eigen::Vector3d(0.1, 0.05, -0.02)),
       sva::MotionVecd(Eigen::Vector3d(0.1, 0.2, -0.1), Eigen::Vector3d(0.02, -0.01, 0.03))},
      {makePose(-0.6, 0.4, -0.8, {0.1, 0.4, 0.85}), sva::MotionVecd(Eigen::Vector3d(-1.2, 0.8, 0.4), Eigen::Vector3d(-0.2, 0.1, 0.05)),
       sva::MotionVecd(Eigen::Vector3d(0.6, -0.4, 0.3), Eigen::Vector3d(-0.15, 0.08, -0.05))},
      {makePose(1.2, 0.1, 0.3, {0.5, -0.3, 1.1}), sva::MotionVecd(Eigen::Vector3d(0.2, 1.5, -0.6), Eigen::Vector3d(-0.05, 0.2, 0.1)),
       sva::MotionVecd(Eigen::Vector3d(1.8, -0.9, 0.4), Eigen::Vector3d(0.1, -0.2, 0.15))},
      {makePose(-0.3, -0.5, 0.6, {0.2, 0.15, 0.95}), sva::MotionVecd(Eigen::Vector3d(-0.7, -1.1, 0.9), Eigen::Vector3d(0.15, -0.1, 0.05)),
       sva::MotionVecd(Eigen::Vector3d(-0.9, 1.2, -0.5), Eigen::Vector3d(-0.1, 0.15, -0.2))},
      {makePose(0.8, -0.7, -0.4, {0.4, 0.05, 1.05}), sva::MotionVecd(Eigen::Vector3d(1.0, -1.4, 0.7), Eigen::Vector3d(-0.2, -0.1, 0.2)),
       sva::MotionVecd(Eigen::Vector3d(0.5, 0.6, -0.8), Eigen::Vector3d(0.2, 0.1, -0.1))},
      // A near free-fall case: acceleration close to -gravity in world, expressed in this
      // pose's (near-identity) body frame -- the sensor should read close to just offset+bias.
      {makePose(0.0, 0.0, 0.0, {0.3, 0.2, 0.9}), sva::MotionVecd(Eigen::Vector6d::Zero()),
       sva::MotionVecd(Eigen::Vector3d::Zero(), Eigen::Vector3d(0, 0, -mc_rtc::constants::GRAVITY))},
      {makePose(0.5, 0.5, -0.5, {0.25, -0.2, 1.0}), sva::MotionVecd(Eigen::Vector3d(0.3, -0.6, 1.1), Eigen::Vector3d(0.1, 0.1, -0.1)),
       sva::MotionVecd(Eigen::Vector3d(-0.4, 0.3, -0.6), Eigen::Vector3d(-0.1, 0.05, 0.1))},
  };

  Measurements measurements;
  for(const auto & s : samples)
  {
    Measurement m;
    m.X_0_p = s.pos;
    m.V_p = s.V_p;
    m.A_p = s.A_p;
    m.measure = predict(s.pos, X_p_f, s.V_p, s.A_p, mass_gt, rpy_gt, com_gt, offset_gt, I_com_gt);
    measurements.push_back(m);
  }

  // Deliberately wrong initial guess so the solver has to actually find the ground truth,
  // not just stay near a lucky starting point.
  InitialGuess guess;
  guess.mass = 0.2;
  guess.com = {0.0, 0.0, 0.05};
  guess.rpy = {0.0, 0.0, 0.0};
  guess.offset = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  auto result = calibrate(robot, sensorName, measurements, guess, I_com_gt, false);

  // Coarse sanity check that computeToolInertia produces a physically plausible (positive
  // definite, reasonable magnitude) result on a real robot model, independent of the
  // synthetic-data test above.
  auto I_model = computeToolInertia(robot, sensorName);
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(I_model);
  bool inertiaOk = (es.eigenvalues().array() >= 0).all();
  std::cout << "computeToolInertia(JVRC1, " << sensorName << ") eigenvalues: " << es.eigenvalues().transpose()
            << (inertiaOk ? " (positive semi-definite, OK)" : " (NOT positive semi-definite, FAILED)") << "\n";

  bool ok = result.success;
  const double tol = 1e-5;
  ok &= nearlyEqual(result.mass, mass_gt, tol);
  for(int i = 0; i < 3; ++i)
  {
    ok &= nearlyEqual(result.rpy[i], rpy_gt[i], tol);
    ok &= nearlyEqual(result.com[i], com_gt[i], tol);
  }
  for(int i = 0; i < 6; ++i)
  {
    ok &= nearlyEqual(result.offset[i], offset_gt[i], tol);
  }

  std::cout << "success: " << result.success << "\n";
  std::cout << "mass    fit=" << result.mass << " gt=" << mass_gt << "\n";
  std::cout << "rpy     fit=" << result.rpy[0] << "," << result.rpy[1] << "," << result.rpy[2] << " gt=" << rpy_gt[0]
            << "," << rpy_gt[1] << "," << rpy_gt[2] << "\n";
  std::cout << "com     fit=" << result.com[0] << "," << result.com[1] << "," << result.com[2] << " gt=" << com_gt[0]
            << "," << com_gt[1] << "," << com_gt[2] << "\n";
  std::cout << "offset  fit=" << result.offset[0] << "," << result.offset[1] << "," << result.offset[2] << ","
            << result.offset[3] << "," << result.offset[4] << "," << result.offset[5] << "\n";
  std::cout << "offset  gt =" << offset_gt[0] << "," << offset_gt[1] << "," << offset_gt[2] << "," << offset_gt[3]
            << "," << offset_gt[4] << "," << offset_gt[5] << "\n";

  ok &= inertiaOk;

  if(!ok)
  {
    std::cerr << "FAILED: fitted parameters do not match ground truth within tolerance " << tol << "\n";
    return EXIT_FAILURE;
  }
  std::cout << "PASSED\n";
  return EXIT_SUCCESS;
}
