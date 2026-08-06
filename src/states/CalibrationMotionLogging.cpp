#include "CalibrationMotionLogging.h"

#include <mc_rbdyn/rpy_utils.h>
#include <mc_rtc/io_utils.h>

#include <RBDyn/FA.h>
#include <RBDyn/FK.h>
#include <RBDyn/FV.h>

#include <boost/filesystem.hpp>

#include <map>
#include <utility>

#include "../ForceSensorCalibration.h"
#include "../Measurement.h"

namespace bfs = boost::filesystem;

void CalibrationMotionLogging::start(mc_control::fsm::Controller & ctl_)
{
  auto & ctl = static_cast<ForceSensorCalibration &>(ctl_);

  auto & robot = ctl.robot();
  auto robotConf = ctl.config()("robots")(robot.name());

  singularityThreshold_ = robotConf("SingularityThreshold", ctl.config()("SingularityThreshold"));

  if(!robotConf.has("forceSensors"))
  {
    mc_rtc::log::error_and_throw<std::runtime_error>("Calibration controller expects a forceSensors entry");
  }
  sensors_ = robotConf("forceSensors");
  auto outputPath_ = (bfs::temp_directory_path() / fmt::format("calib-force-sensors-data-{}", robot.name())).string();

  // Attempt to create the output files
  if(!boost::filesystem::exists(outputPath_))
  {
    if(!boost::filesystem::create_directory(outputPath_))
    {
      mc_rtc::log::error_and_throw<std::runtime_error>("[{}] Could not create output folder {}", name(), outputPath_);
    }
  }

  if(!ctl.datastore().has("measurements"))
  {
    ctl.datastore().make<SensorMeasurements>("measurements");
  }
  auto & measurements = ctl.datastore().get<SensorMeasurements>("measurements");
  measurements.clear();
  for(const auto & s : sensors_)
  {
    measurements[s] = {};
    measurementsCount_[s] = {0, 0};
    if(singularityThreshold_ > 0)
    {
      const auto & sensor = ctl_.robot().forceSensor(s);
      jacobians_[s] = {ctl_.robot(), sensor.parentBody()};
    }
  }
}

bool CalibrationMotionLogging::run(mc_control::fsm::Controller & ctl_)
{
  auto & robot = ctl_.robot();
  auto & real = ctl_.realRobot();

  // Spatial velocity/acceleration of each sensor's parent body, computed from the
  // analytically-known commanded joint velocity/acceleration (CalibrationMotion publishes the
  // exact derivatives of the trajectory it commands, so this is noise-free -- unlike
  // differentiating an actual measured signal) instead of assuming static equilibrium. Joints
  // not part of the calibration motion are held static (zero velocity/acceleration), matching
  // reality since CalibrationMotion only moves the excited joints.
  rbd::MultiBodyConfig mbc = real.mbc();
  for(auto & jv : mbc.alpha)
  {
    std::fill(jv.begin(), jv.end(), 0.0);
  }
  for(auto & jv : mbc.alphaD)
  {
    std::fill(jv.begin(), jv.end(), 0.0);
  }
  if(ctl_.datastore().has("CalibrationMotion::JointDerivatives"))
  {
    const auto & derivatives = ctl_.datastore().get<std::map<std::string, std::pair<double, double>>>(
        "CalibrationMotion::JointDerivatives");
    for(const auto & [jointName, qdqdd] : derivatives)
    {
      auto jidx = real.jointIndexByName(jointName);
      mbc.alpha[jidx][0] = qdqdd.first;
      mbc.alphaD[jidx][0] = qdqdd.second;
    }
  }
  rbd::forwardKinematics(real.mb(), mbc);
  rbd::forwardVelocity(real.mb(), mbc);
  rbd::forwardAcceleration(real.mb(), mbc);

  auto & measurements = ctl_.datastore().get<SensorMeasurements>("measurements");
  for(const auto & s : sensors_)
  {
    const auto & sensor = robot.forceSensor(s);
    auto bodyIdx = real.bodyIndexByName(sensor.parentBody());
    const auto & X_0_p = real.bodyPosW()[bodyIdx];
    const auto & measure = sensor.wrench();
    const auto & V_p = mbc.bodyVelB[bodyIdx];
    const auto & A_p = mbc.bodyAccB[bodyIdx];
    measurementsCount_[s].second++;
    if(singularityThreshold_ > 0)
    {
      auto & j = jacobians_[s];
      const auto & jacMat = j.jacobian.jacobian(robot.mb(), robot.mbc());
      j.svd.compute(jacMat);
      if(j.svd.singularValues().tail(1)(0) > singularityThreshold_)
      {
        measurements[s].push_back({X_0_p, measure, V_p, A_p});
        measurementsCount_[s].first++;
      }
    }
    else
    {
      measurements[s].push_back({X_0_p, measure, V_p, A_p});
      measurementsCount_[s].first++;
    }
  }
  output("OK");
  return true;
}

void CalibrationMotionLogging::teardown(mc_control::fsm::Controller &)
{
  for(const auto & m : measurementsCount_)
  {
    mc_rtc::log::info("[{}] {} records taken out of {} iterations", m.first, m.second.first, m.second.second);
  }
}

EXPORT_SINGLE_STATE("CalibrationMotionLogging", CalibrationMotionLogging)
