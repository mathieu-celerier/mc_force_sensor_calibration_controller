#pragma once

#include <mc_rtc/constants.h>
#include <SpaceVecAlg/SpaceVecAlg>

/** We redefined sva::Rot. functions to make them work with non-scalar types */

template<typename T>
inline Eigen::Matrix3<T> RotX(T theta)
{
  T s = sin(theta), c = cos(theta);
  return (Eigen::Matrix3<T>() << T(1), T(0), T(0), T(0), c, s, T(0), -s, c).finished();
}

template<typename T>
inline Eigen::Matrix3<T> RotY(T theta)
{
  T s = sin(theta), c = cos(theta);
  return (Eigen::Matrix3<T>() << c, T(0), -s, T(0), T(1), T(0), s, T(0), c).finished();
}

template<typename T>
inline Eigen::Matrix3<T> RotZ(T theta)
{
  T s = sin(theta), c = cos(theta);
  return (Eigen::Matrix3<T>() << c, s, T(0), -s, c, T(0), T(0), T(0), T(1)).finished();
}

template<typename T>
Eigen::Matrix3<T> rpyToMat(const T & r, const T & p, const T & y)
{
  return RotX<T>(r) * RotY<T>(p) * RotZ<T>(y);
}

/** Cost functor for the force sensor calibration problem.
 *
 * Exposed in its own header (rather than kept private to calibrate.cpp) so that it can be
 * reused, with T=double, to generate self-consistent synthetic measurements in tests: feeding
 * back the exact same model used for the fit is the only way to test the inertial-term
 * frame/sign conventions without a physical rig to compare against.
 */
struct CostFunctor
{
  CostFunctor(const sva::PTransformd & pos,
              const sva::ForceVecd & measure,
              const sva::PTransformd & X_p_f,
              const sva::MotionVecd & V_p,
              const sva::MotionVecd & A_p,
              const Eigen::Matrix3d & I_com = Eigen::Matrix3d::Zero())
  : pos_(pos), measure_(measure), X_p_f_(X_p_f), V_p_(V_p), A_p_(A_p), I_com_(I_com)
  {
  }

  template<typename T>
  bool operator()(const T * const mass, const T * const rpy, const T * const com, const T * const offset, T * residual)
      const
  {
    const T gravity(mc_rtc::constants::GRAVITY);
    sva::ForceVec<T> vf(Eigen::Vector3<T>::Zero(), Eigen::Vector3<T>(T(0), T(0), -mass[0] * gravity));
    sva::PTransform<T> X_s_ds = sva::PTransform<T>(rpyToMat(rpy[0], rpy[1], rpy[2]));
    sva::PTransform<T> X_p_vb = sva::PTransform<T>(Eigen::Vector3<T>(com[0], com[1], com[2]));
    sva::ForceVec<T> off = sva::ForceVec<T>(Eigen::Vector3<T>(offset[0], offset[1], offset[2]),
                                            Eigen::Vector3<T>(offset[3], offset[4], offset[5]));
    sva::PTransform<T> X_0_p = pos_.cast<T>();
    sva::PTransform<T> X_0_vb = sva::PTransform<T>(Eigen::Matrix3<T>::Identity(), (X_p_vb * X_0_p).translation());
    sva::PTransform<T> X_p_ds = X_s_ds * X_p_f_.cast<T>();
    sva::PTransform<T> X_ds_vb = X_0_vb * (X_p_ds * X_0_p).inv();
    sva::ForceVec<T> vb_f = off + X_ds_vb.transMul(vf);

    // Account for the tool's own inertial motion (translational only, the tool's rotational
    // inertia tensor is not a fit parameter) instead of assuming static equilibrium. V_p_/A_p_
    // are the sensor parent body's spatial velocity/(Featherstone) acceleration expressed in
    // its own frame. Since vb (the CoM) is rigidly attached to p at a fixed offset (com), its
    // spatial velocity/acceleration transform exactly via X_p_vb (no relative joint motion).
    sva::MotionVec<T> V_p_T = V_p_.cast<T>();
    sva::MotionVec<T> A_p_T = A_p_.cast<T>();
    sva::MotionVec<T> V_vb = X_p_vb * V_p_T;
    sva::MotionVec<T> A_vb = X_p_vb * A_p_T;
    // Featherstone spatial acceleration -> classical (physical) point acceleration:
    // a_classical = a_spatial.linear() + omega x v_linear (evaluated at the same point).
    Eigen::Vector3<T> a_vb_classical = A_vb.linear() + V_vb.angular().cross(V_vb.linear());
    // Re-express in the sensor's own (ds) frame axes to combine with off/measure_. Only the
    // rotation matters here since acceleration is a free vector.
    Eigen::Vector3<T> a_vb_classical_ds = X_p_ds.rotation() * a_vb_classical;
    // Rotational counterpart: I*alpha + omega x (I*omega), using the tool's own (fixed, not
    // fit) inertia tensor about its CoM. The angular part of a MotionVec is unaffected by a
    // pure-translation PTransform (same rigid body, same rotation), so V_vb/A_vb's angular
    // parts already equal V_p_/A_p_'s.
    Eigen::Matrix3<T> I_com_T = I_com_.cast<T>();
    Eigen::Vector3<T> omega_vb = V_vb.angular();
    Eigen::Vector3<T> alpha_vb = A_vb.angular();
    Eigen::Vector3<T> tau_vb_classical = I_com_T * alpha_vb + omega_vb.cross(I_com_T * omega_vb);
    Eigen::Vector3<T> tau_vb_classical_ds = X_p_ds.rotation() * tau_vb_classical;

    // The sensor measures the reaction to what the payload needs to accelerate: reaction =
    // gravity_load - mass * a_true - (I*alpha + omega x I*omega) (reduces to the static
    // gravity-only case when a_true = 0 and omega = alpha = 0, and correctly predicts zero
    // reading in free fall when a_true = gravity).
    sva::ForceVec<T> inertial_ds(tau_vb_classical_ds, mass[0] * a_vb_classical_ds);
    vb_f = vb_f - inertial_ds;

    sva::ForceVec<T> diff = measure_.cast<T>() - vb_f;
    residual[0] = diff.couple().x();
    residual[1] = diff.couple().y();
    residual[2] = diff.couple().z();
    residual[3] = diff.force().x();
    residual[4] = diff.force().y();
    residual[5] = diff.force().z();
    return true;
  }

private:
  const sva::PTransformd pos_;
  const sva::ForceVecd measure_;
  const sva::PTransformd X_p_f_;
  const sva::MotionVecd V_p_;
  const sva::MotionVecd A_p_;
  const Eigen::Matrix3d I_com_;
};
