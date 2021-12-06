#ifndef __quadruped_walkgen_quadruped_augmented_hxx__
#define __quadruped_walkgen_quadruped_augmented_hxx__

#include "crocoddyl/core/utils/exception.hpp"

namespace quadruped_walkgen {
template <typename Scalar>
ActionModelQuadrupedAugmentedTpl<Scalar>::ActionModelQuadrupedAugmentedTpl(
    typename Eigen::Matrix<Scalar, 3, 1> offset_CoM)
    : crocoddyl::ActionModelAbstractTpl<Scalar>(boost::make_shared<crocoddyl::StateVectorTpl<Scalar> >(20), 12, 32) {
  // Relative forces to compute the norm mof the command
  relative_forces = true;
  uref_.setZero();

  // Model parameters
  mu = Scalar(1);
  dt_ = Scalar(0.02);
  mass = Scalar(2.50000279);
  min_fz_in_contact = Scalar(0.0);
  max_fz_in_contact = Scalar(25.0);

  // Matrix model initialization
  g.setZero();
  g[8] = Scalar(-9.81) * dt_;
  gI.setZero();
  gI.diagonal() << Scalar(0.00578574), Scalar(0.01938108), Scalar(0.02476124);
  A.setIdentity();
  A.topRightCorner(6, 6) << Eigen::Matrix<Scalar, 6, 6>::Identity() * dt_;
  B.setZero();
  lever_arms.setZero();
  R.setZero();

  // Weight vectors initialization
  force_weights_.setConstant(Scalar(0.2));
  state_weights_ << Scalar(1.), Scalar(1.), Scalar(150.), Scalar(35.), Scalar(30.), Scalar(8.), Scalar(20.),
      Scalar(20.), Scalar(15.), Scalar(4.), Scalar(4.), Scalar(8.);
  friction_weight_ = Scalar(10);
  heuristic_weights_.setConstant(Scalar(1));
  stop_weights_.setConstant(Scalar(1));
  // pshoulder_ << Scalar(0.1946), Scalar(0.15005), Scalar(0.1946), Scalar(-0.15005), Scalar(-0.1946), Scalar(0.15005),
  //     Scalar(-0.1946), Scalar(-0.15005);
  pshoulder_0 << Scalar(0.18), Scalar(0.18), Scalar(-0.21), Scalar(-0.21), Scalar(0.14695), Scalar(-0.14695),
      Scalar(0.14695), Scalar(-0.14695);
  // pshoulder_tmp.setZero();
  // pcentrifugal_tmp_1.setZero();
  // pcentrifugal_tmp_2.setZero();
  // pcentrifugal_tmp.setZero();
  // UpperBound vector
  ub.setZero();
  for (int i = 0; i < 4; i = i + 1) {
    ub(6 * i + 5) = max_fz_in_contact;
  }

  // Temporary vector used
  Fa_x_u.setZero();
  rub_.setZero();
  rub_max_.setZero();
  Arr.setZero();
  r.setZero();
  lever_tmp.setZero();
  R_tmp.setZero();
  gait.setZero();
  base_vector_x << Scalar(1.), Scalar(0.), Scalar(0.);
  base_vector_y << Scalar(0.), Scalar(1.), Scalar(0.);
  base_vector_z << Scalar(0.), Scalar(0.), Scalar(1.);
  forces_3d.setZero();
  gait_double.setZero();

  // bool to add heuristic for foot position
  centrifugal_term = true;
  symmetry_term = true;
  T_gait = Scalar(0.48);

  // // Used for shoulder height weight
  // pshoulder_0 <<  Scalar(0.1946) ,   Scalar(0.1946) ,   Scalar(-0.1946),  Scalar(-0.1946) ,
  //                 Scalar(0.15005) ,  Scalar(-0.15005)  , Scalar(0.15005)  ,  Scalar(-0.15005) ;
  sh_hlim = Scalar(0.27);
  sh_weight.setConstant(Scalar(1.));
  sh_ub_max_.setZero();
  psh.setZero();
  pheuristic_.setZero();
  offset_com = offset_CoM;  // x, y, z offset

  shoulder_reference_position = false;  // Using predicted trajectory of the CoM
}

template <typename Scalar>
ActionModelQuadrupedAugmentedTpl<Scalar>::~ActionModelQuadrupedAugmentedTpl() {}

template <typename Scalar>
void ActionModelQuadrupedAugmentedTpl<Scalar>::calc(
    const boost::shared_ptr<crocoddyl::ActionDataAbstractTpl<Scalar> >& data,
    const Eigen::Ref<const typename MathBase::VectorXs>& x, const Eigen::Ref<const typename MathBase::VectorXs>& u) {
  if (static_cast<std::size_t>(x.size()) != state_->get_nx()) {
    throw_pretty("Invalid argument: "
                 << "x has wrong dimension (it should be " + std::to_string(state_->get_nx()) + ")");
  }
  if (static_cast<std::size_t>(u.size()) != nu_) {
    throw_pretty("Invalid argument: "
                 << "u has wrong dimension (it should be " + std::to_string(nu_) + ")");
  }

  ActionDataQuadrupedAugmentedTpl<Scalar>* d = static_cast<ActionDataQuadrupedAugmentedTpl<Scalar>*>(data.get());

  //  Update B :
  for (int i = 0; i < 4; i = i + 1) {
    lever_tmp.setZero();
    if (gait(i, 0) != 0) {
      lever_tmp.head(2) = x.block(12 + 2 * i, 0, 2, 1);
      lever_tmp += -x.block(0, 0, 3, 1);
      R_tmp << Scalar(0.0), -lever_tmp[2], lever_tmp[1], lever_tmp[2], Scalar(0.0), -lever_tmp[0], -lever_tmp[1],
          lever_tmp[0], Scalar(0.0);
      B.block(9, 3 * i, 3, 3) << dt_ * R * R_tmp;

      // Compute pdistance of the shoulder wrt contact point
      if (shoulder_reference_position) {
        // Ref vector as reference for the shoulder trajectory, roll and pitch at first
        psh.block(0, i, 3, 1) << xref_(0, 0) - offset_com(0, 0) + pshoulder_0(0, i) * cos(xref_(5, 0)) -
                                     pshoulder_0(1, i) * sin(xref_(5, 0)) - x(12 + 2 * i),
            xref_(1, 0) - offset_com(1, 0) + pshoulder_0(0, i) * sin(xref_(5, 0)) +
                pshoulder_0(1, i) * cos(xref_(5, 0)) - x(12 + 2 * i + 1),
            xref_(2, 0) - offset_com(2, 0);
      } else {
        // psh.block(0, i, 3, 1) << x(0) + pshoulder_0(0, i) - pshoulder_0(1, i) * x(5) - x(12 + 2 * i),
        //                          x(1) + pshoulder_0(1, i) + pshoulder_0(0, i) * x(5) - x(12 + 2 * i + 1),
        //                          x(2) - offset_com + pshoulder_0(1, i) * x(3) - pshoulder_0(0, i) * x(4);
        // Correction, no approximation for yaw
        psh.block(0, i, 3, 1) << x(0) - offset_com(0, 0) + pshoulder_0(0, i) * cos(x(5)) -
                                     pshoulder_0(1, i) * sin(x(5)) - x(12 + 2 * i),
            x(1) - offset_com(1, 0) + pshoulder_0(0, i) * sin(x(5)) + pshoulder_0(1, i) * cos(x(5)) -
                x(12 + 2 * i + 1),
            x(2) - offset_com(2, 0) + pshoulder_0(1, i) * x(3) - pshoulder_0(0, i) * x(4);
      }

    } else {
      // Compute pdistance of the shoulder wrt contact point
      psh.block(0, i, 3, 1).setZero();
    }
  };

  // Discrete dynamic : A*x + B*u + g
  d->xnext.template head<12>() = A.diagonal().cwiseProduct(x.block(0, 0, 12, 1)) + g;
  d->xnext.template head<6>() =
      d->xnext.template head<6>() + A.topRightCorner(6, 6).diagonal().cwiseProduct(x.block(6, 0, 6, 1));
  d->xnext.template segment<6>(6) = d->xnext.template segment<6>(6) + B.block(6, 0, 6, 12) * u;
  d->xnext.template tail<8>() = x.tail(8);

  // Residual cost on the state and force norm
  d->r.template head<12>() = state_weights_.cwiseProduct(x.head(12) - xref_);
  d->r.template segment<8>(12) =
      ((heuristic_weights_.cwiseProduct(x.tail(8) - pheuristic_)).array() * gait_double.array()).matrix();
  d->r.template tail<12>() = force_weights_.cwiseProduct(u - uref_);

  // Friction cone
  for (int i = 0; i < 4; i = i + 1) {
    Fa_x_u.segment(6 * i, 6) << u(3 * i) - mu * u(3 * i + 2), -u(3 * i) - mu * u(3 * i + 2),
        u(3 * i + 1) - mu * u(3 * i + 2), -u(3 * i + 1) - mu * u(3 * i + 2), -u(3 * i + 2), u(3 * i + 2);
  }
  rub_max_ = (Fa_x_u - ub).cwiseMax(Scalar(0.));

  // Shoulder height weight
  sh_ub_max_ << Scalar(0.5) * sh_weight(0) * (psh.block(0, 0, 3, 1).squaredNorm() - sh_hlim * sh_hlim),
      Scalar(0.5) * sh_weight(1) * (psh.block(0, 1, 3, 1).squaredNorm() - sh_hlim * sh_hlim),
      Scalar(0.5) * sh_weight(2) * (psh.block(0, 2, 3, 1).squaredNorm() - sh_hlim * sh_hlim),
      Scalar(0.5) * sh_weight(3) * (psh.block(0, 3, 3, 1).squaredNorm() - sh_hlim * sh_hlim);

  sh_ub_max_ = sh_ub_max_.cwiseMax(Scalar(0.));

  // Cost computation
  // d->cost = Scalar(0.5) * d->r.transpose() * d->r     + friction_weight_ * Scalar(0.5) * rub_max_.squaredNorm()
  //         + Scalar(0.5)*( (stop_weights_.cwiseProduct(x.tail(8) - pref_) ).array() * gait_double.array()
  //         ).matrix().squaredNorm()  ;

  d->cost =
      Scalar(0.5) * d->r.transpose() * d->r + friction_weight_ * Scalar(0.5) * rub_max_.squaredNorm() +
      Scalar(0.5) *
          ((stop_weights_.cwiseProduct(x.tail(8) - pstop_)).array() * gait_double.array()).matrix().squaredNorm() +
      sh_ub_max_.sum();
}

template <typename Scalar>
void ActionModelQuadrupedAugmentedTpl<Scalar>::calcDiff(
    const boost::shared_ptr<crocoddyl::ActionDataAbstractTpl<Scalar> >& data,
    const Eigen::Ref<const typename MathBase::VectorXs>& x, const Eigen::Ref<const typename MathBase::VectorXs>& u) {
  if (static_cast<std::size_t>(x.size()) != state_->get_nx()) {
    throw_pretty("Invalid argument: "
                 << "x has wrong dimension (it should be " + std::to_string(state_->get_nx()) + ")");
  }
  if (static_cast<std::size_t>(u.size()) != nu_) {
    throw_pretty("Invalid argument: "
                 << "u has wrong dimension (it should be " + std::to_string(nu_) + ")");
  }

  ActionDataQuadrupedAugmentedTpl<Scalar>* d = static_cast<ActionDataQuadrupedAugmentedTpl<Scalar>*>(data.get());

  // Cost derivatives : Lx
  d->Lx.setZero();
  d->Lx.template head<12>() = (state_weights_.array() * d->r.template head<12>().array()).matrix();
  d->Lx.template tail<8>() =
      (heuristic_weights_.array() * d->r.template segment<8>(12).array()).matrix() +
      ((stop_weights_.cwiseProduct(x.tail(8) - pstop_)).array() * gait_double.array() * stop_weights_.array())
          .matrix();

  // Cost derivative : Lu
  for (int i = 0; i < 4; i = i + 1) {
    r = friction_weight_ * rub_max_.segment(6 * i, 6);
    d->Lu.block(i * 3, 0, 3, 1) << r(0) - r(1), r(2) - r(3), -mu * (r(0) + r(1) + r(2) + r(3)) - r(4) + r(5);
  }
  d->Lu = d->Lu + (force_weights_.array() * d->r.template tail<12>().array()).matrix();

  // Hessian : Lxx
  // Hessian : Lxx
  d->Lxx.setZero();

  d->Lxx.diagonal().head(12) = (state_weights_.array() * state_weights_.array()).matrix();
  d->Lxx.diagonal().tail(8) = (gait_double.array() * heuristic_weights_.array() * heuristic_weights_.array()).matrix();
  d->Lxx.diagonal().tail(8) += (gait_double.array() * stop_weights_.array() * stop_weights_.array()).matrix();

  // Shoulder height derivative cost
  for (int j = 0; j < 4; j = j + 1) {
    if (sh_ub_max_[j] > Scalar(0.)) {
      if (shoulder_reference_position) {
        d->Lx(12 + 2 * j, 0) += -sh_weight(j) * psh(0, j);
        d->Lx(12 + 2 * j + 1, 0) += -sh_weight(j) * psh(1, j);

        d->Lxx(12 + 2 * j, 12 + 2 * j) += sh_weight(j);
        d->Lxx(12 + 2 * j + 1, 12 + 2 * j + 1) += sh_weight(j);

      } else {
        // Approximation of small even for yaw (wrong)
        // d->Lx(0, 0) += sh_weight(j) * psh(0, j);
        // d->Lx(1, 0) += sh_weight(j) * psh(1, j);
        // d->Lx(2, 0) += sh_weight(j) * psh(2, j);
        // d->Lx(3, 0) += sh_weight(j) * pshoulder_0(1, j) * psh(2, j);
        // d->Lx(4, 0) += -sh_weight(j) * pshoulder_0(0, j) * psh(2, j);
        // d->Lx(5, 0) += sh_weight(j) * (-pshoulder_0(1, j) * psh(0, j) + pshoulder_0(0, j) * psh(1, j));

        // d->Lx(12 + 2 * j, 0) += -sh_weight(j) * psh(0, j);
        // d->Lx(12 + 2 * j + 1, 0) += -sh_weight(j) * psh(1, j);

        // d->Lxx(0, 0) += sh_weight(j);
        // d->Lxx(1, 1) += sh_weight(j);
        // d->Lxx(2, 2) += sh_weight(j);
        // d->Lxx(3, 3) += sh_weight(j) * pshoulder_0(1, j) * pshoulder_0(1, j);
        // d->Lxx(3, 3) += sh_weight(j) * pshoulder_0(0, j) * pshoulder_0(0, j);
        // d->Lxx(5, 5) += sh_weight(j) * (pshoulder_0(1, j) * pshoulder_0(1, j) + pshoulder_0(0, j) * pshoulder_0(0,
        // j));

        // d->Lxx(12 + 2 * j, 12 + 2 * j) += sh_weight(j);
        // d->Lxx(12 + 2 * j + 1, 12 + 2 * j + 1) += sh_weight(j);

        // d->Lxx(0, 5) += -sh_weight(j) * pshoulder_0(1, j);
        // d->Lxx(5, 0) += -sh_weight(j) * pshoulder_0(1, j);

        // d->Lxx(1, 5) += sh_weight(j) * pshoulder_0(0, j);
        // d->Lxx(5, 1) += sh_weight(j) * pshoulder_0(0, j);

        // d->Lxx(2, 3) += sh_weight(j) * pshoulder_0(1, j);
        // d->Lxx(2, 4) += -sh_weight(j) * pshoulder_0(0, j);
        // d->Lxx(3, 2) += sh_weight(j) * pshoulder_0(1, j);
        // d->Lxx(4, 2) += -sh_weight(j) * pshoulder_0(0, j);

        // d->Lxx(3, 4) += -sh_weight(j) * pshoulder_0(1, j) * pshoulder_0(0, j);
        // d->Lxx(4, 3) += -sh_weight(j) * pshoulder_0(1, j) * pshoulder_0(0, j);

        // d->Lxx(0, 12 + 2 * j) += -sh_weight(j);
        // d->Lxx(12 + 2 * j, 0) += -sh_weight(j);

        // d->Lxx(5, 12 + 2 * j) += sh_weight(j) * pshoulder_0(1, j);
        // d->Lxx(12 + 2 * j, 5) += sh_weight(j) * pshoulder_0(1, j);

        // d->Lxx(1, 12 + 2 * j + 1) += -sh_weight(j);
        // d->Lxx(12 + 2 * j + 1, 1) += -sh_weight(j);

        // d->Lxx(5, 12 + 2 * j + 1) += -sh_weight(j) * pshoulder_0(0, j);
        // d->Lxx(12 + 2 * j + 1, 5) += -sh_weight(j) * pshoulder_0(0, j);
        d->Lx(0, 0) += sh_weight(j) * psh(0, j);
        d->Lx(1, 0) += sh_weight(j) * psh(1, j);
        d->Lx(2, 0) += sh_weight(j) * psh(2, j);
        d->Lx(3, 0) += sh_weight(j) * pshoulder_0(1, j) * psh(2, j);
        d->Lx(4, 0) += -sh_weight(j) * pshoulder_0(0, j) * psh(2, j);
        d->Lx(5, 0) += sh_weight(j) * ((-sin(x(5)) * pshoulder_0(0, j) - cos(x(5)) * pshoulder_0(1, j)) * psh(0, j) +
                                       (cos(x(5)) * pshoulder_0(0, j) - sin(x(5)) * pshoulder_0(1, j)) * psh(1, j));

        d->Lx(12 + 2 * j, 0) += -sh_weight(j) * psh(0, j);
        d->Lx(12 + 2 * j + 1, 0) += -sh_weight(j) * psh(1, j);

        d->Lxx(0, 0) += sh_weight(j);
        d->Lxx(1, 1) += sh_weight(j);
        d->Lxx(2, 2) += sh_weight(j);
        d->Lxx(3, 3) += sh_weight(j) * pshoulder_0(1, j) * pshoulder_0(1, j);
        d->Lxx(3, 3) += sh_weight(j) * pshoulder_0(0, j) * pshoulder_0(0, j);
        d->Lxx(5, 5) += sh_weight(j) * ((-cos(x(5)) * pshoulder_0(0, j) + sin(x(5)) * pshoulder_0(1, j)) * psh(0, j) +
                                        (-sin(x(5)) * pshoulder_0(0, j) - cos(x(5)) * pshoulder_0(1, j)) *
                                            (-sin(x(5)) * pshoulder_0(0, j) - cos(x(5)) * pshoulder_0(1, j)) +
                                        (-sin(x(5)) * pshoulder_0(0, j) - cos(x(5)) * pshoulder_0(1, j)) * psh(1, j) +
                                        (cos(x(5)) * pshoulder_0(0, j) - sin(x(5)) * pshoulder_0(1, j)) *
                                            (cos(x(5)) * pshoulder_0(0, j) - sin(x(5)) * pshoulder_0(1, j)));

        d->Lxx(12 + 2 * j, 12 + 2 * j) += sh_weight(j);
        d->Lxx(12 + 2 * j + 1, 12 + 2 * j + 1) += sh_weight(j);

        d->Lxx(0, 5) += sh_weight(j) * (-sin(x(5)) * pshoulder_0(0, j) - cos(x(5)) * pshoulder_0(1, j));
        d->Lxx(5, 0) += sh_weight(j) * (-sin(x(5)) * pshoulder_0(0, j) - cos(x(5)) * pshoulder_0(1, j));

        d->Lxx(1, 5) += sh_weight(j) * (cos(x(5)) * pshoulder_0(0, j) - sin(x(5)) * pshoulder_0(1, j));
        d->Lxx(5, 1) += sh_weight(j) * (cos(x(5)) * pshoulder_0(0, j) - sin(x(5)) * pshoulder_0(1, j));

        d->Lxx(2, 3) += sh_weight(j) * pshoulder_0(1, j);
        d->Lxx(2, 4) += -sh_weight(j) * pshoulder_0(0, j);
        d->Lxx(3, 2) += sh_weight(j) * pshoulder_0(1, j);
        d->Lxx(4, 2) += -sh_weight(j) * pshoulder_0(0, j);

        d->Lxx(3, 4) += -sh_weight(j) * pshoulder_0(1, j) * pshoulder_0(0, j);
        d->Lxx(4, 3) += -sh_weight(j) * pshoulder_0(1, j) * pshoulder_0(0, j);

        d->Lxx(0, 12 + 2 * j) += -sh_weight(j);
        d->Lxx(12 + 2 * j, 0) += -sh_weight(j);

        d->Lxx(5, 12 + 2 * j) += -sh_weight(j) * (-sin(x(5)) * pshoulder_0(0, j) - cos(x(5)) * pshoulder_0(1, j));
        d->Lxx(12 + 2 * j, 5) += -sh_weight(j) * (-sin(x(5)) * pshoulder_0(0, j) - cos(x(5)) * pshoulder_0(1, j));

        d->Lxx(1, 12 + 2 * j + 1) += -sh_weight(j);
        d->Lxx(12 + 2 * j + 1, 1) += -sh_weight(j);

        d->Lxx(5, 12 + 2 * j + 1) += -sh_weight(j) * (cos(x(5)) * pshoulder_0(0, j) - sin(x(5)) * pshoulder_0(1, j));
        d->Lxx(12 + 2 * j + 1, 5) += -sh_weight(j) * (cos(x(5)) * pshoulder_0(0, j) - sin(x(5)) * pshoulder_0(1, j));
      }
    }
  }

  // Hessian : Luu
  // Matrix friction cone hessian (20x12)
  Arr.diagonal() = ((Fa_x_u - ub).array() >= 0.).matrix().template cast<Scalar>();
  for (int i = 0; i < 4; i = i + 1) {
    r = friction_weight_ * Arr.diagonal().segment(6 * i, 6);
    d->Luu.block(3 * i, 3 * i, 3, 3) << r(0) + r(1), 0.0, mu * (r(1) - r(0)), 0.0, r(2) + r(3), mu * (r(3) - r(2)),
        mu * (r(1) - r(0)), mu * (r(3) - r(2)), mu * mu * (r(0) + r(1) + r(2) + r(3)) + r(4) + r(5);
  }
  d->Luu.diagonal() = d->Luu.diagonal() + (force_weights_.array() * force_weights_.array()).matrix();

  // Dynamic derivatives
  d->Fx.setZero();
  d->Fx.block(0, 0, 12, 12) << A;
  d->Fx.block(12, 12, 8, 8) << Eigen::Matrix<Scalar, 8, 8>::Identity();

  for (int i = 0; i < 4; i = i + 1) {
    if (gait(i, 0) != 0) {
      forces_3d = u.block(3 * i, 0, 3, 1);
      d->Fx.block(9, 0, 3, 1) += -dt_ * R * (base_vector_x.cross(forces_3d));
      d->Fx.block(9, 1, 3, 1) += -dt_ * R * (base_vector_y.cross(forces_3d));
      d->Fx.block(9, 2, 3, 1) += -dt_ * R * (base_vector_z.cross(forces_3d));

      d->Fx.block(9, 12 + 2 * i, 3, 1) += dt_ * R * (base_vector_x.cross(forces_3d));
      d->Fx.block(9, 12 + 2 * i + 1, 3, 1) += dt_ * R * (base_vector_y.cross(forces_3d));
    }
  }
  // d->Fu << Eigen::Matrix<Scalar, 20, 12>::Zero() ;
  d->Fu.block(0, 0, 12, 12) << B;
}

template <typename Scalar>
boost::shared_ptr<crocoddyl::ActionDataAbstractTpl<Scalar> > ActionModelQuadrupedAugmentedTpl<Scalar>::createData() {
  return boost::make_shared<ActionDataQuadrupedAugmentedTpl<Scalar> >(this);
}

////////////////////////////////
// get & set parameters ////////
////////////////////////////////

template <typename Scalar>
const typename Eigen::Matrix<Scalar, 12, 1>& ActionModelQuadrupedAugmentedTpl<Scalar>::get_force_weights() const {
  return force_weights_;
}
template <typename Scalar>
void ActionModelQuadrupedAugmentedTpl<Scalar>::set_force_weights(const typename MathBase::VectorXs& weights) {
  if (static_cast<std::size_t>(weights.size()) != 12) {
    throw_pretty("Invalid argument: "
                 << "Weights vector has wrong dimension (it should be 12)");
  }
  force_weights_ = weights;
}

template <typename Scalar>
const typename Eigen::Matrix<Scalar, 12, 1>& ActionModelQuadrupedAugmentedTpl<Scalar>::get_state_weights() const {
  return state_weights_;
}
template <typename Scalar>
void ActionModelQuadrupedAugmentedTpl<Scalar>::set_state_weights(const typename MathBase::VectorXs& weights) {
  if (static_cast<std::size_t>(weights.size()) != 12) {
    throw_pretty("Invalid argument: "
                 << "Weights vector has wrong dimension (it should be 12)");
  }
  state_weights_ = weights;
}

template <typename Scalar>
const typename Eigen::Matrix<Scalar, 8, 1>& ActionModelQuadrupedAugmentedTpl<Scalar>::get_heuristic_weights() const {
  return heuristic_weights_;
}
template <typename Scalar>
void ActionModelQuadrupedAugmentedTpl<Scalar>::set_heuristic_weights(const typename MathBase::VectorXs& weights) {
  if (static_cast<std::size_t>(weights.size()) != 8) {
    throw_pretty("Invalid argument: "
                 << "Weights vector has wrong dimension (it should be 8)");
  }
  heuristic_weights_ = weights;
}

template <typename Scalar>
const typename Eigen::Matrix<Scalar, 8, 1>& ActionModelQuadrupedAugmentedTpl<Scalar>::get_stop_weights() const {
  return stop_weights_;
}
template <typename Scalar>
void ActionModelQuadrupedAugmentedTpl<Scalar>::set_stop_weights(const typename MathBase::VectorXs& weights) {
  if (static_cast<std::size_t>(weights.size()) != 8) {
    throw_pretty("Invalid argument: "
                 << "Weights vector has wrong dimension (it should be 8)");
  }
  stop_weights_ = weights;
}

template <typename Scalar>
const Scalar& ActionModelQuadrupedAugmentedTpl<Scalar>::get_friction_weight() const {
  return friction_weight_;
}
template <typename Scalar>
void ActionModelQuadrupedAugmentedTpl<Scalar>::set_friction_weight(const Scalar& weight) {
  friction_weight_ = weight;
}

template <typename Scalar>
const Scalar& ActionModelQuadrupedAugmentedTpl<Scalar>::get_mu() const {
  return mu;
}
template <typename Scalar>
void ActionModelQuadrupedAugmentedTpl<Scalar>::set_mu(const Scalar& mu_coeff) {
  mu = mu_coeff;
}

template <typename Scalar>
const Scalar& ActionModelQuadrupedAugmentedTpl<Scalar>::get_mass() const {
  return mass;
}
template <typename Scalar>
void ActionModelQuadrupedAugmentedTpl<Scalar>::set_mass(const Scalar& m) {
  // The model need to be updated after this changed
  mass = m;
}

template <typename Scalar>
const Scalar& ActionModelQuadrupedAugmentedTpl<Scalar>::get_dt() const {
  return dt_;
}
template <typename Scalar>
void ActionModelQuadrupedAugmentedTpl<Scalar>::set_dt(const Scalar& dt) {
  // The model need to be updated after this changed
  dt_ = dt;
  g[8] = Scalar(-9.81) * dt_;
  A.topRightCorner(6, 6) << Eigen::Matrix<Scalar, 6, 6>::Identity() * dt_;
}

template <typename Scalar>
const typename Eigen::Matrix<Scalar, 3, 3>& ActionModelQuadrupedAugmentedTpl<Scalar>::get_gI() const {
  return gI;
}
template <typename Scalar>
void ActionModelQuadrupedAugmentedTpl<Scalar>::set_gI(const typename MathBase::Matrix3s& inertia_matrix) {
  // The model need to be updated after this changed
  if (static_cast<std::size_t>(inertia_matrix.size()) != 9) {
    throw_pretty("Invalid argument: "
                 << "gI has wrong dimension : 3x3");
  }
  gI = inertia_matrix;
}

template <typename Scalar>
const Scalar& ActionModelQuadrupedAugmentedTpl<Scalar>::get_min_fz_contact() const {
  // The model need to be updated after this changed
  return min_fz_in_contact;
}
template <typename Scalar>
void ActionModelQuadrupedAugmentedTpl<Scalar>::set_min_fz_contact(const Scalar& min_fz) {
  // The model need to be updated after this changed
  min_fz_in_contact = min_fz;
}

template <typename Scalar>
const Scalar& ActionModelQuadrupedAugmentedTpl<Scalar>::get_max_fz_contact() const {
  // The model need to be updated after this changed
  return max_fz_in_contact;
}
template <typename Scalar>
void ActionModelQuadrupedAugmentedTpl<Scalar>::set_max_fz_contact(const Scalar& max_fz) {
  // The model need to be updated after this changed
  max_fz_in_contact = max_fz;
  for (int i = 0; i < 4; i = i + 1) {
    ub(6 * i + 5) = max_fz_in_contact;
  }
}

template <typename Scalar>
const bool& ActionModelQuadrupedAugmentedTpl<Scalar>::get_symmetry_term() const {
  return symmetry_term;
}
template <typename Scalar>
void ActionModelQuadrupedAugmentedTpl<Scalar>::set_symmetry_term(const bool& sym_term) {
  // The model need to be updated after this changed
  symmetry_term = sym_term;
}

template <typename Scalar>
const bool& ActionModelQuadrupedAugmentedTpl<Scalar>::get_centrifugal_term() const {
  return centrifugal_term;
}
template <typename Scalar>
void ActionModelQuadrupedAugmentedTpl<Scalar>::set_centrifugal_term(const bool& cent_term) {
  // The model need to be updated after this changed
  centrifugal_term = cent_term;
}

template <typename Scalar>
const bool& ActionModelQuadrupedAugmentedTpl<Scalar>::get_shoulder_reference_position() const {
  return shoulder_reference_position;
}
template <typename Scalar>
void ActionModelQuadrupedAugmentedTpl<Scalar>::set_shoulder_reference_position(const bool& reference) {
  shoulder_reference_position = reference;
}

template <typename Scalar>
const Scalar& ActionModelQuadrupedAugmentedTpl<Scalar>::get_T_gait() const {
  // The model need to be updated after this changed
  return T_gait;
}
template <typename Scalar>
void ActionModelQuadrupedAugmentedTpl<Scalar>::set_T_gait(const Scalar& T_gait_) {
  // The model need to be updated after this changed
  T_gait = T_gait_;
}

template <typename Scalar>
const Scalar& ActionModelQuadrupedAugmentedTpl<Scalar>::get_shoulder_hlim() const {
  return sh_hlim;
}
template <typename Scalar>
void ActionModelQuadrupedAugmentedTpl<Scalar>::set_shoulder_hlim(const Scalar& hlim) {
  // The model need to be updated after this changed
  sh_hlim = hlim;
}

template <typename Scalar>
const typename Eigen::Matrix<Scalar, 4, 1>& ActionModelQuadrupedAugmentedTpl<Scalar>::get_shoulder_contact_weight()
    const {
  return sh_weight;
}
template <typename Scalar>
void ActionModelQuadrupedAugmentedTpl<Scalar>::set_shoulder_contact_weight(
    const typename Eigen::Matrix<Scalar, 4, 1>& weight) {
  // The model need to be updated after this changed
  sh_weight = weight;
}

///////////////////////////
//// get A & B matrix /////
///////////////////////////
template <typename Scalar>
const typename Eigen::Matrix<Scalar, 12, 12>& ActionModelQuadrupedAugmentedTpl<Scalar>::get_A() const {
  return A;
}
template <typename Scalar>
const typename Eigen::Matrix<Scalar, 12, 12>& ActionModelQuadrupedAugmentedTpl<Scalar>::get_B() const {
  return B;
}

// to modify the cost on the command : || fz - m*g/nb contact ||^2
// --> set to True
template <typename Scalar>
const bool& ActionModelQuadrupedAugmentedTpl<Scalar>::get_relative_forces() const {
  return relative_forces;
}
template <typename Scalar>
void ActionModelQuadrupedAugmentedTpl<Scalar>::set_relative_forces(const bool& rel_forces) {
  relative_forces = rel_forces;
  uref_.setZero();
  if (relative_forces) {
    for (int i = 0; i < 4; i = i + 1) {
      if (gait[i] == 1) {
        uref_[3 * i + 2] = (Scalar(9.81) * mass) / (gait.sum());
      }
    }
  }
}

////////////////////////
// Update current model
////////////////////////

template <typename Scalar>
void ActionModelQuadrupedAugmentedTpl<Scalar>::update_model(
    const Eigen::Ref<const typename MathBase::MatrixXs>& l_feet,
    const Eigen::Ref<const typename MathBase::MatrixXs>& l_stop,
    const Eigen::Ref<const typename MathBase::MatrixXs>& xref,
    const Eigen::Ref<const typename MathBase::MatrixXs>& S) {
  if (static_cast<std::size_t>(l_feet.size()) != 12) {
    throw_pretty("Invalid argument: "
                 << "l_feet matrix has wrong dimension (it should be : 3x4)");
  }
  if (static_cast<std::size_t>(xref.size()) != 12) {
    throw_pretty("Invalid argument: "
                 << "xref vector has wrong dimension (it should be 12 )");
  }
  if (static_cast<std::size_t>(S.size()) != 4) {
    throw_pretty("Invalid argument: "
                 << "S vector has wrong dimension (it should be 4x1)");
  }

  xref_ = xref;
  gait = S;

  uref_.setZero();
  if (relative_forces) {
    for (int i = 0; i < 4; i = i + 1) {
      if (gait[i] == 1) {
        uref_[3 * i + 2] = (Scalar(9.81) * mass) / (gait.sum());
      }
    }
  }

  for (int i = 0; i < 4; i = i + 1) {
    gait_double(2 * i, 0) = gait(i, 0);
    gait_double(2 * i + 1, 0) = gait(i, 0);

    pheuristic_.block(2 * i, 0, 2, 1) = l_feet.block(0, i, 2, 1);
    pstop_.block(2 * i, 0, 2, 1) = l_stop.block(0, i, 2, 1);
  }

  R_tmp << cos(xref(5, 0)), -sin(xref(5, 0)), Scalar(0), sin(xref(5, 0)), cos(xref(5, 0)), Scalar(0), Scalar(0),
      Scalar(0), Scalar(1.0);

  // Centrifual term
  // pcentrifugal_tmp_1 = xref.block(6, 0, 3, 1);
  // pcentrifugal_tmp_2 = xref.block(9, 0, 3, 1);
  // pcentrifugal_tmp = 0.5 * std::sqrt(xref(2, 0) / 9.81) * pcentrifugal_tmp_1.cross(pcentrifugal_tmp_2);

  // for (int i = 0; i < 4; i = i + 1) {
  //   pshoulder_tmp.block(0, i, 2, 1) =
  //       R_tmp.block(0, 0, 2, 2) *
  //       (pshoulder_0.block(0, i, 2, 1) + symmetry_term * 0.25 * T_gait * xref.block(6, 0, 2, 1) +
  //        centrifugal_term * pcentrifugal_tmp.block(0, 0, 2, 1));
  // }

  R = (R_tmp.transpose() * gI * R_tmp).inverse();  // I_inv

  for (int i = 0; i < 4; i = i + 1) {
    // pshoulder_[2 * i] = pshoulder_tmp(0, i) + xref(0, 0);
    // pshoulder_[2 * i + 1] = pshoulder_tmp(1, i) + xref(1, 0);

    if (S(i, 0) != 0) {
      // set limit for normal force, (foot in contact with the ground)
      ub[5 * i + 4] = -min_fz_in_contact;

      // B update
      B.block(6, 3 * i, 3, 3).diagonal() << dt_ / mass, dt_ / mass, dt_ / mass;

      //  Assuption 1 : levers arms not depends on the state, but on the predicted position (xfref)
      //  --> B will be updated with the update_B method for each calc function

      // lever_tmp = lever_arms.block(0,i,3,1) - xref.block(0,0,3,1) ;
      // R_tmp << 0.0, -lever_tmp[2], lever_tmp[1],
      // lever_tmp[2], 0.0, -lever_tmp[0], -lever_tmp[1], lever_tmp[0], 0.0 ;
      // B.block(9 , 3*i  , 3,3) << dt_ * R* R_tmp;
    } else {
      // set limit for normal force at 0.0
      ub[5 * i + 4] = Scalar(0.0);
      B.block(6, 3 * i, 3, 3).setZero();
      B.block(9, 3 * i, 3, 3).setZero();
    };
  };
}
}  // namespace quadruped_walkgen

#endif
