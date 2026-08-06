#pragma once
#include <mc_rbdyn/Robot.h>
#include <fmt/format.h> // for fmt::format / format_to
#include <fmt/ostream.h> // enables fmt fallback to operator<<
#include <ostream>
#include <sstream>

#include "Measurement.h"

struct InitialGuess
{
  double mass = 0;
  std::array<double, 3> rpy = {0};
  std::array<double, 3> com = {0};
  std::array<double, 6> offset = {0};
};

inline std::ostream & operator<<(std::ostream & os, const InitialGuess & r)
{
  // clang-format off
  fmt::format_to(std::ostream_iterator<char>(os),
R"(mass        : {}
rpy         : {}, {}, {}
com         : {}, {}, {}
force offset: {}, {}, {}, {}, {}, {})",
      r.mass,
      r.rpy[0], r.rpy[1], r.rpy[2],
      r.com[0], r.com[1], r.com[2],
      r.offset[0], r.offset[1], r.offset[2], r.offset[3], r.offset[4], r.offset[5]);
  // clang-format on
  return os;
}

template<typename Char>
struct fmt::formatter<InitialGuess, Char> : fmt::formatter<std::basic_string_view<Char>, Char>
{
  template<typename FormatContext>
  auto format(const InitialGuess & v, FormatContext & ctx) const
  {
    std::ostringstream os;
    os << v;
    auto s = os.str();
    return fmt::formatter<std::basic_string_view<Char>, Char>::format(s, ctx);
  }
};

namespace mc_rtc
{
template<>
struct ConfigurationLoader<InitialGuess>
{
  static InitialGuess load(const mc_rtc::Configuration & config)
  {
    InitialGuess r;
    config("mass", r.mass);
    config("com", r.com);
    config("rpy", r.rpy);
    config("offset", r.offset);
    return r;
  }

  static mc_rtc::Configuration save(const InitialGuess & r)
  {
    mc_rtc::Configuration c;
    c.add("mass", r.mass);
    c.add("com", r.com);
    c.add("rpy", r.rpy);
    c.add("offset", r.offset);
    return c;
  }
};
} // namespace mc_rtc

struct CalibrationResult
{
  CalibrationResult(const InitialGuess & guess)
  {
    mass = guess.mass;
    rpy = guess.rpy;
    com = guess.com;
    offset = guess.offset;
  }

  bool success = false;
  double mass = 0;
  std::array<double, 3> rpy = {0};
  std::array<double, 3> com = {0};
  std::array<double, 6> offset = {0};
};

InitialGuess computeInitialGuessFromModel(const mc_rbdyn::Robot & robot,
                                          const std::string & sensor,
                                          bool includeParent = false,
                                          bool verbose = false);

/** Rotational inertia of the tool (everything downstream of the force sensor, same bodies as
 * computeInitialGuessFromModel), about its own center of mass, expressed in the sensor's
 * parent body frame axes. Used as a fixed (non-fit) input to account for the tool's rotational
 * inertial wrench (I*alpha + omega x I*omega) during the calibration motion instead of only
 * its translational inertia (mass*a_com). Only an approximation once the fitted CoM diverges
 * from the model's nominal CoM, but far better than omitting rotational inertia entirely. */
Eigen::Matrix3d computeToolInertia(const mc_rbdyn::Robot & robot,
                                   const std::string & sensor,
                                   bool includeParent = false);

CalibrationResult calibrate(const mc_rbdyn::Robot & robot,
                            const std::string & sensor,
                            const Measurements & measurements,
                            const InitialGuess & initialGuess,
                            const Eigen::Matrix3d & I_com = Eigen::Matrix3d::Zero(),
                            bool verbose = false);
