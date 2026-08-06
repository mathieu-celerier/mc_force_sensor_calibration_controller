#include "CalibrationMotion.h"
#include <mc_filter/utils/clamp.h>
#include "../ForceSensorCalibration.h"

void CalibrationMotion::start(mc_control::fsm::Controller & ctl)
{
  ctl.datastore().make_call("CalibrationMotion::Stop", [this]() { interrupted_ = true; });
  ctl.datastore().make<std::map<std::string, std::pair<double, double>>>("CalibrationMotion::JointDerivatives");
  auto & robot = ctl.robot();
  auto robotConf = ctl.config()("robots")(robot.name());
  if(!robotConf.has("motion"))
  {
    mc_rtc::log::error("[{}] Calibration controller expects a joints entry", name());
    output("FAILURE");
  }
  auto conf = robotConf("motion");
  conf("duration", duration_);
  conf("percentLimits", percentLimits_);
  mc_filter::utils::clampInPlace(percentLimits_, 0, 1);

  auto postureTask = ctl.getPostureTask(robot.name());
  savedStiffness_ = postureTask->stiffness();
  postureTask->stiffness(conf("stiffness", 10));
  constexpr double PI = mc_rtc::constants::PI;
  for(const auto & jConfig : conf("joints"))
  {
    std::string name = jConfig("name");
    if(!ctl.robot().hasJoint(name))
    {
      mc_rtc::log::error("[{}] No joint named \"{}\" in robot \"{}\"", this->name(), name, ctl.robot().name());
      output("FAILURE");
    }
    auto percentLimits = percentLimits_;
    jConfig("percentLimits", percentLimits);
    mc_filter::utils::clampInPlace(percentLimits, 0, 1);
    double period = jConfig("period");
    auto jidx = robot.jointIndexByName(name);
    auto start = robot.mbc().q[jidx][0];
    auto actualLower = robot.ql()[jidx][0];
    auto actualUpper = robot.qu()[jidx][0];
    auto actualRange = actualUpper - actualLower;

    // Reduced range, symmetric around the joint's range midpoint by default
    const auto range = percentLimits * actualRange;
    auto lower = actualLower + (actualRange - range) / 2;
    auto upper = actualUpper - (actualRange - range) / 2;

    // Optional explicit bounds (e.g. to keep an asymmetric collision-free window instead of
    // the symmetric percentLimits shrink above) override the computed lower/upper
    jConfig("lower", lower);
    jConfig("upper", upper);
    mc_filter::utils::clampInPlace(lower, actualLower, actualUpper);
    mc_filter::utils::clampInPlace(upper, actualLower, actualUpper);

    if(lower >= upper)
    {
      mc_rtc::log::error("[{}] Invalid motion bounds for joint {}: lower ({}) must be less than upper ({})",
                         this->name(), name, lower, upper);
      output("FAILURE");
    }

    if(start < lower || start > upper)
    {
      mc_rtc::log::error("[{}] Starting joint configuration of joint {} [{}] is outside of the reduced limit range "
                         "[{}, {}] (percentLimits: {}, actual joint limits: [{}, {}]",
                         this->name(), name, start, lower, upper, percentLimits, actualLower, actualUpper);
      output("FAILURE");
    }

    // compute the starting time such that the joint does not move initially
    // that is such that f(start_dt) = start
    // i.e start_dt = f^(-1)(start)
    double start_dt = period * (acos(sqrt(start - lower) / sqrt(upper - lower))) / PI;
    jointUpdates_.emplace_back(
        /* f(t): periodic function that moves the joint between its limits, along with its
         * analytically-known velocity/acceleration (exact, no differentiation of a measured
         * signal needed since we command this trajectory) exposed for CalibrationMotionLogging
         * to account for the tool's own inertial motion instead of assuming static
         * equilibrium */
        [this, &ctl, postureTask, lower, upper, start_dt, period, name]()
        {
          auto t = start_dt + dt_;
          auto w = (2 * PI) / period;
          auto q = lower + (upper - lower) * (1 + cos(w * t)) / 2;
          auto qd = -(upper - lower) / 2 * w * sin(w * t);
          auto qdd = -(upper - lower) / 2 * w * w * cos(w * t);
          postureTask->target({{name, {q}}});
          auto & derivatives = ctl.datastore().get<std::map<std::string, std::pair<double, double>>>(
              "CalibrationMotion::JointDerivatives");
          derivatives[name] = {qd, qdd};
        });
  }

  ctl.gui()->addElement(
      {},
      mc_rtc::gui::NumberSlider(
          "Progress", [this]() { return dt_; }, [](double) {}, 0, duration_),
      mc_rtc::gui::Button("Stop Motion",
                          [this]()
                          {
                            mc_rtc::log::warning(
                                "[{}] Motion was interrupted before it's planned duration ({:.2f}/{:.2f}s)", name(),
                                dt_, duration_);
                            interrupted_ = true;
                          }));
}

bool CalibrationMotion::run(mc_control::fsm::Controller & ctl_)
{
  if(output() == "FAILURE")
  {
    return true;
  }

  // Update all joint positions
  for(auto & updateJoint : jointUpdates_)
  {
    updateJoint();
  }

  if(dt_ > duration_)
  {
    dt_ = duration_;
    auto postureTask = ctl_.getPostureTask(ctl_.robot().name());
    postureTask->refVel(Eigen::VectorXd{ctl_.robot().mb().nrDof()}.setZero());
    if(postureTask->speed().norm() < 1e-5)
    {
      output("OK");
      return true;
    }
  }
  else if(interrupted_)
  {
    output("INTERRUPTED");
    return true;
  }

  dt_ += ctl_.timeStep;
  return false;
}

void CalibrationMotion::teardown(mc_control::fsm::Controller & ctl_)
{
  auto postureTask = ctl_.getPostureTask(ctl_.robot().name());
  postureTask->stiffness(savedStiffness_);
  ctl_.gui()->removeElement({}, "Progress");
  ctl_.gui()->removeElement({}, "Stop Motion");
  ctl_.datastore().remove("CalibrationMotion::Stop");
  ctl_.datastore().remove("CalibrationMotion::JointDerivatives");
}

EXPORT_SINGLE_STATE("CalibrationMotion", CalibrationMotion)
