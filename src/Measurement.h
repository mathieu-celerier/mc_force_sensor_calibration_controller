#pragma once

#include <SpaceVecAlg/SpaceVecAlg>

#include <map>
#include <string>
#include <vector>

/** A measurement for the calibration */
struct Measurement
{
  /** Position of the force sensor's parent body */
  sva::PTransformd X_0_p;
  /** Raw force sensor reading */
  sva::ForceVecd measure;
  /** Spatial velocity of the force sensor's parent body, expressed in the body's own frame
   * (i.e. rbd::MultiBodyConfig::bodyVelB convention). Used to account for the tool's own
   * inertial motion during the calibration motion instead of assuming static equilibrium.
   * Zero if the motion providing this measurement was not tracked (e.g. static calibration). */
  sva::MotionVecd V_p = sva::MotionVecd(Eigen::Vector6d::Zero());
  /** Spatial (Featherstone) acceleration of the force sensor's parent body, expressed in the
   * body's own frame (i.e. rbd::MultiBodyConfig::bodyAccB convention, gravity excluded). */
  sva::MotionVecd A_p = sva::MotionVecd(Eigen::Vector6d::Zero());
};

using Measurements = std::vector<Measurement>;
using SensorMeasurements = std::map<std::string, Measurements>;
