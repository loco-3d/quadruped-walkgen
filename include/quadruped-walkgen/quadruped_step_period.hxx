#ifndef __quadruped_walkgen_quadruped_step_period_hxx__
#define __quadruped_walkgen_quadruped_step_period_hxx__

#include "crocoddyl/core/utils/exception.hpp"

namespace quadruped_walkgen {
template <typename Scalar>
ActionModelQuadrupedStepPeriodTpl<Scalar>::ActionModelQuadrupedStepPeriodTpl()
    : crocoddyl::ActionModelAbstractTpl<Scalar>(boost::make_shared<crocoddyl::StateVectorTpl<Scalar> >(21), 5, 26) {
  B.setZero();
  dt_ref_.setConstant(Scalar(0.02));
  rub_max_.setZero();
  rub_max_bool.setZero();
  dt_min_.setConstant(Scalar(0.005));
  dt_max_.setConstant(Scalar(0.1));
  dt_weight_ = Scalar(1);
  dt_bound_weight = Scalar(10);
  state_weights_ << Scalar(1.), Scalar(1.), Scalar(150.), Scalar(35.), Scalar(30.), Scalar(8.), Scalar(20.),
      Scalar(20.), Scalar(15.), Scalar(4.), Scalar(4.), Scalar(8.);
  shoulder_weights_.setConstant(Scalar(1));
  pshoulder_ << Scalar(0.1946), Scalar(0.15005), Scalar(0.1946), Scalar(-0.15005), Scalar(-0.1946), Scalar(0.15005),
      Scalar(-0.1946), Scalar(-0.15005);
  pshoulder_0 << Scalar(0.1946), Scalar(0.1946), Scalar(-0.1946), Scalar(-0.1946), Scalar(0.15005), Scalar(-0.15005),
      Scalar(0.15005), Scalar(-0.15005);
  pshoulder_tmp.setZero();
  pcentrifugal_tmp_1.setZero();
  pcentrifugal_tmp_2.setZero();
  pcentrifugal_tmp.setZero();
  centrifugal_term = true;
  symmetry_term = true;
  T_gait = Scalar(0.64);

  step_weights_.setConstant(Scalar(1));

  // Optim dt
  nb_nodes = Scalar(15.);
  vlim = Scalar(2.);
  beta_lim = Scalar((64 * nb_nodes * nb_nodes * vlim * vlim) / 225);  // apparent speed used in the cost function
  speed_weight = Scalar(10.);
}

template <typename Scalar>
ActionModelQuadrupedStepPeriodTpl<Scalar>::~ActionModelQuadrupedStepPeriodTpl() {}

template <typename Scalar>
void ActionModelQuadrupedStepPeriodTpl<Scalar>::calc(
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

  ActionDataQuadrupedStepPeriodTpl<Scalar>* d = static_cast<ActionDataQuadrupedStepPeriodTpl<Scalar>*>(data.get());

  d->xnext.template head<12>() = x.head(12);
  d->xnext.template segment<8>(12) = x.segment(12, 8) + B * u.head(4);
  d->xnext.template tail<1>() = u.tail(1);

  // Residual cost on the state and force norm
  d->r.template head<12>() = state_weights_.cwiseProduct(x.head(12) - xref_);
  d->r.template segment<8>(12) = shoulder_weights_.cwiseProduct(x.segment(12, 8) - pshoulder_);
  d->r.template segment<1>(20) = dt_weight_ * (x.tail(1) - dt_ref_);
  d->r.template segment<4>(21) = step_weights_.cwiseProduct(u.head(4));
  d->r.template tail<1>() = dt_weight_ * (u.tail(1) - dt_ref_);

  rub_max_ << dt_min_ - x.tail(1), x.tail(1) - dt_max_, u[0] * u[0] + u[1] * u[1] - beta_lim * x[20] * x[20],
      u[2] * u[2] + u[3] * u[3] - beta_lim * x[20] * x[20];

  rub_max_bool = (rub_max_.array() >= Scalar(0.)).matrix().template cast<Scalar>();
  rub_max_ = rub_max_.cwiseMax(Scalar(0.));

  // d->cost = Scalar(0.5) * d->r.transpose() * d->r   + dt_bound_weight * Scalar(0.5) * rub_max_.head(2).squaredNorm()
  // + speed_weight * Scalar(0.5) * rub_max_.tail(2).squaredNorm();
  d->cost = Scalar(0.5) * d->r.transpose() * d->r + dt_bound_weight * Scalar(0.5) * rub_max_.head(2).squaredNorm() +
            speed_weight * Scalar(0.5) * rub_max_.tail(2).sum();
}

template <typename Scalar>
void ActionModelQuadrupedStepPeriodTpl<Scalar>::calcDiff(
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

  ActionDataQuadrupedStepPeriodTpl<Scalar>* d = static_cast<ActionDataQuadrupedStepPeriodTpl<Scalar>*>(data.get());

  // Cost derivatives : Lx
  d->Lx.template head<12>() = (state_weights_.array() * d->r.template head<12>().array()).matrix();
  d->Lx.template segment<8>(12) = (shoulder_weights_.array() * d->r.template segment<8>(12).array()).matrix();
  d->Lx.template tail<1>() << dt_bound_weight * (-rub_max_[0] + rub_max_[1]) -
                                  beta_lim * speed_weight * x(20) * rub_max_bool[2] -
                                  beta_lim * speed_weight * x(20) * rub_max_bool[3];
  d->Lx.template tail<1>() += dt_weight_ * d->r.template segment<1>(20);

  // cost period <--> distance
  // d->Lu << speed_weight*Scalar(2)*u[0]*rub_max_[2] , speed_weight*Scalar(2)*u[1]*rub_max_[2] ,
  //          speed_weight*Scalar(2)*u[2]*rub_max_[3] , speed_weight*Scalar(2)*u[3]*rub_max_[3] ,
  //          - speed_weight*Scalar(2)*beta_lim*u[4]*(rub_max_[2] + rub_max_[3]);
  d->Lu << speed_weight * u[0] * rub_max_bool[2], speed_weight * u[1] * rub_max_bool[2],
      speed_weight * u[2] * rub_max_bool[3], speed_weight * u[3] * rub_max_bool[3], Scalar(0.);

  d->Lu.template head<4>() += (step_weights_.array() * d->r.template segment<4>(21).array()).matrix();
  d->Lu.template tail<1>() += dt_weight_ * d->r.template tail<1>();

  // Hessian : Lxx
  d->Lxx.diagonal().head(12) = (state_weights_.array() * state_weights_.array()).matrix();
  d->Lxx.diagonal().segment(12, 8) = (shoulder_weights_.array() * shoulder_weights_.array()).matrix();
  d->Lxx.diagonal().tail(1) << dt_weight_ * dt_weight_ + dt_bound_weight * rub_max_bool[0] +
                                   dt_bound_weight * rub_max_bool[1];

  d->Lxx(20, 20) += -beta_lim * speed_weight * rub_max_bool[2] - beta_lim * speed_weight * rub_max_bool[3];

  d->Luu.diagonal() << speed_weight * rub_max_bool[2], speed_weight * rub_max_bool[2], speed_weight * rub_max_bool[3],
      speed_weight * rub_max_bool[3], Scalar(0.);

  d->Luu.diagonal().head(4) += (step_weights_.array() * step_weights_.array()).matrix();

  // d->Luu.diagonal() << speed_weight*Scalar(2)*(rub_max_[2] + 2*u[0]*u[0] )*rub_max_bool[2] ,
  //                      speed_weight*Scalar(2)*(rub_max_[2] + 2*u[1]*u[1] )*rub_max_bool[2] ,
  //                      speed_weight*Scalar(2)*(rub_max_[3] + 2*u[2]*u[2] )*rub_max_bool[3] ,
  //                      speed_weight*Scalar(2)*(rub_max_[3] + 2*u[3]*u[3] )*rub_max_bool[3] ,
  //                      dt_weight_*dt_weight_ - speed_weight*Scalar(2)*beta_lim*(rub_max_[2] -
  //                      2*beta_lim*u[4]*u[4])*rub_max_bool[2] -speed_weight*Scalar(2)*beta_lim*(rub_max_[3] -
  //                      2*beta_lim*u[4]*u[4])*rub_max_bool[3]  ;
  // d->Luu.diagonal().head(4) += (step_weights_.array() * step_weights_.array()).matrix() ;
  // d->Luu.block(0,1,1,1) << speed_weight*Scalar(4)*u[0]*u[1]*rub_max_bool[2] ;
  // d->Luu.block(1,0,1,1) << speed_weight*Scalar(4)*u[0]*u[1]*rub_max_bool[2] ;
  // d->Luu.block(0,4,1,1) << -speed_weight*Scalar(4)*beta_lim*u[0]*u[4]*rub_max_bool[2] ;
  // d->Luu.block(4,0,1,1) << -speed_weight*Scalar(4)*beta_lim*u[0]*u[4]*rub_max_bool[2] ;
  // d->Luu.block(1,4,1,1) << -speed_weight*Scalar(4)*beta_lim*u[1]*u[4]*rub_max_bool[2] ;
  // d->Luu.block(4,1,1,1) << -speed_weight*Scalar(4)*beta_lim*u[1]*u[4]*rub_max_bool[2] ;
  // d->Luu.block(2,3,1,1) << speed_weight*Scalar(4)*u[2]*u[3]*rub_max_bool[3] ;
  // d->Luu.block(3,2,1,1) << speed_weight*Scalar(4)*u[2]*u[3]*rub_max_bool[3] ;
  // d->Luu.block(0,4,1,1) << -speed_weight*Scalar(4)*beta_lim*u[2]*u[4]*rub_max_bool[3] ;
  // d->Luu.block(4,0,1,1) << -speed_weight*Scalar(4)*beta_lim*u[2]*u[4]*rub_max_bool[3] ;
  // d->Luu.block(1,4,1,1) << -speed_weight*Scalar(4)*beta_lim*u[3]*u[4]*rub_max_bool[3] ;
  // d->Luu.block(4,1,1,1) << -speed_weight*Scalar(4)*beta_lim*u[3]*u[4]*rub_max_bool[3] ;

  // Dynamic derivatives
  d->Fx.setIdentity();
  d->Fx.block(20, 20, 1, 1) << 0;
  d->Fu.block(12, 0, 8, 4) = B;
  d->Fu.block(20, 4, 1, 1) << 1;
}

template <typename Scalar>
boost::shared_ptr<crocoddyl::ActionDataAbstractTpl<Scalar> > ActionModelQuadrupedStepPeriodTpl<Scalar>::createData() {
  return boost::make_shared<ActionDataQuadrupedStepPeriodTpl<Scalar> >(this);
}

////////////////////////////////
// get & set parameters ////////
////////////////////////////////

template <typename Scalar>
const typename Eigen::Matrix<Scalar, 12, 1>& ActionModelQuadrupedStepPeriodTpl<Scalar>::get_state_weights() const {
  return state_weights_;
}
template <typename Scalar>
void ActionModelQuadrupedStepPeriodTpl<Scalar>::set_state_weights(const typename MathBase::VectorXs& weights) {
  if (static_cast<std::size_t>(weights.size()) != 12) {
    throw_pretty("Invalid argument: "
                 << "Weights vector has wrong dimension (it should be 12)");
  }
  state_weights_ = weights;
}

template <typename Scalar>
const typename Eigen::Matrix<Scalar, 4, 1>& ActionModelQuadrupedStepPeriodTpl<Scalar>::get_step_weights() const {
  return step_weights_;
}
template <typename Scalar>
void ActionModelQuadrupedStepPeriodTpl<Scalar>::set_step_weights(const typename MathBase::VectorXs& weights) {
  if (static_cast<std::size_t>(weights.size()) != 4) {
    throw_pretty("Invalid argument: "
                 << "Weights vector has wrong dimension (it should be 4)");
  }
  step_weights_ = weights;
}

template <typename Scalar>
const typename Eigen::Matrix<Scalar, 8, 1>& ActionModelQuadrupedStepPeriodTpl<Scalar>::get_shoulder_weights() const {
  return shoulder_weights_;
}
template <typename Scalar>
void ActionModelQuadrupedStepPeriodTpl<Scalar>::set_shoulder_weights(const typename MathBase::VectorXs& weights) {
  if (static_cast<std::size_t>(weights.size()) != 8) {
    throw_pretty("Invalid argument: "
                 << "Weights vector has wrong dimension (it should be 8)");
  }
  shoulder_weights_ = weights;
}

template <typename Scalar>
const typename Eigen::Matrix<Scalar, 8, 1>& ActionModelQuadrupedStepPeriodTpl<Scalar>::get_shoulder_position() const {
  return pshoulder_;
}
template <typename Scalar>
void ActionModelQuadrupedStepPeriodTpl<Scalar>::set_shoulder_position(const typename MathBase::VectorXs& pos) {
  if (static_cast<std::size_t>(pos.size()) != 8) {
    throw_pretty("Invalid argument: "
                 << "Weights vector has wrong dimension (it should be 8)");
  }
  pshoulder_ = pos;
}

template <typename Scalar>
const bool& ActionModelQuadrupedStepPeriodTpl<Scalar>::get_symmetry_term() const {
  return symmetry_term;
}
template <typename Scalar>
void ActionModelQuadrupedStepPeriodTpl<Scalar>::set_symmetry_term(const bool& sym_term) {
  // The model need to be updated after this changed
  symmetry_term = sym_term;
}

template <typename Scalar>
const bool& ActionModelQuadrupedStepPeriodTpl<Scalar>::get_centrifugal_term() const {
  return centrifugal_term;
}
template <typename Scalar>
void ActionModelQuadrupedStepPeriodTpl<Scalar>::set_centrifugal_term(const bool& cent_term) {
  // The model need to be updated after this changed
  centrifugal_term = cent_term;
}

template <typename Scalar>
const Scalar& ActionModelQuadrupedStepPeriodTpl<Scalar>::get_T_gait() const {
  // The model need to be updated after this changed
  return T_gait;
}
template <typename Scalar>
void ActionModelQuadrupedStepPeriodTpl<Scalar>::set_T_gait(const Scalar& T_gait_) {
  // The model need to be updated after this changed
  T_gait = T_gait_;
}

template <typename Scalar>
const Scalar& ActionModelQuadrupedStepPeriodTpl<Scalar>::get_dt_weight() const {
  // The model need to be updated after this changed
  return dt_weight_;
}
template <typename Scalar>
void ActionModelQuadrupedStepPeriodTpl<Scalar>::set_dt_weight(const Scalar& weight_) {
  // The model need to be updated after this changed
  dt_weight_ = weight_;
}

template <typename Scalar>
const Scalar& ActionModelQuadrupedStepPeriodTpl<Scalar>::get_speed_weight() const {
  // The model need to be updated after this changed
  return speed_weight;
}
template <typename Scalar>
void ActionModelQuadrupedStepPeriodTpl<Scalar>::set_speed_weight(const Scalar& weight_) {
  // The model need to be updated after this changed
  speed_weight = weight_;
}

template <typename Scalar>
const Scalar& ActionModelQuadrupedStepPeriodTpl<Scalar>::get_dt_bound_weight() const {
  // The model need to be updated after this changed
  return dt_bound_weight;
}
template <typename Scalar>
void ActionModelQuadrupedStepPeriodTpl<Scalar>::set_dt_bound_weight(const Scalar& weight_) {
  // The model need to be updated after this changed
  dt_bound_weight = weight_;
}

template <typename Scalar>
const Scalar& ActionModelQuadrupedStepPeriodTpl<Scalar>::get_nb_nodes() const {
  // The model need to be updated after this changed
  return nb_nodes;
}
template <typename Scalar>
void ActionModelQuadrupedStepPeriodTpl<Scalar>::set_nb_nodes(const Scalar& nodes_) {
  // The model need to be updated after this changed
  nb_nodes = nodes_;
  beta_lim = Scalar((64 * nb_nodes * nb_nodes * vlim * vlim) / 225);
  ;
}

template <typename Scalar>
const Scalar& ActionModelQuadrupedStepPeriodTpl<Scalar>::get_vlim() const {
  // The model need to be updated after this changed
  return vlim;
}
template <typename Scalar>
void ActionModelQuadrupedStepPeriodTpl<Scalar>::set_vlim(const Scalar& vlim_) {
  // The model need to be updated after this changed
  vlim = vlim_;
  beta_lim = Scalar((64 * nb_nodes * nb_nodes * vlim * vlim) / 225);
  ;
}

template <typename Scalar>
const Scalar& ActionModelQuadrupedStepPeriodTpl<Scalar>::get_dt_ref() const {
  return dt_ref_[0];
}
template <typename Scalar>
void ActionModelQuadrupedStepPeriodTpl<Scalar>::set_dt_ref(const Scalar& dt) {
  // The model need to be updated after this changed
  dt_ref_[0] = dt;
}

template <typename Scalar>
const Scalar& ActionModelQuadrupedStepPeriodTpl<Scalar>::get_dt_min() const {
  return dt_min_[0];
}
template <typename Scalar>
void ActionModelQuadrupedStepPeriodTpl<Scalar>::set_dt_min(const Scalar& dt) {
  // The model need to be updated after this changed
  dt_min_[0] = dt;
}

template <typename Scalar>
const Scalar& ActionModelQuadrupedStepPeriodTpl<Scalar>::get_dt_max() const {
  return dt_max_[0];
}
template <typename Scalar>
void ActionModelQuadrupedStepPeriodTpl<Scalar>::set_dt_max(const Scalar& dt) {
  // The model need to be updated after this changed
  dt_max_[0] = dt;
}

////////////////////////
// Update current model
////////////////////////

template <typename Scalar>
void ActionModelQuadrupedStepPeriodTpl<Scalar>::update_model(
    const Eigen::Ref<const typename MathBase::MatrixXs>& l_feet,
    const Eigen::Ref<const typename MathBase::MatrixXs>& xref,
    const Eigen::Ref<const typename MathBase::MatrixXs>& S) {
  if (static_cast<std::size_t>(l_feet.size()) != 12) {
    throw_pretty("Invalid argument: "
                 << "l_feet matrix has wrong dimension (it should be : 3x4)");
  }
  if (static_cast<std::size_t>(xref.size()) != 12) {
    throw_pretty("Invalid argument: "
                 << "Weights vector has wrong dimension (it should be " + std::to_string(state_->get_nx()) + ")");
  }
  if (static_cast<std::size_t>(S.size()) != 4) {
    throw_pretty("Invalid argument: "
                 << "S vector has wrong dimension (it should be 4x1)");
  }

  xref_ = xref;

  R_tmp << cos(xref(5, 0)), -sin(xref(5, 0)), Scalar(0), sin(xref(5, 0)), cos(xref(5, 0)), Scalar(0), Scalar(0),
      Scalar(0), Scalar(1.0);

  // Centrifual term
  pcentrifugal_tmp_1 = xref.block(6, 0, 3, 1);
  pcentrifugal_tmp_2 = xref.block(9, 0, 3, 1);
  pcentrifugal_tmp = 0.5 * std::sqrt(xref(2, 0) / 9.81) * pcentrifugal_tmp_1.cross(pcentrifugal_tmp_2);

  for (int i = 0; i < 4; i = i + 1) {
    pshoulder_tmp.block(0, i, 2, 1) =
        R_tmp.block(0, 0, 2, 2) *
        (pshoulder_0.block(0, i, 2, 1) + symmetry_term * 0.25 * T_gait * xref.block(6, 0, 2, 1) +
         centrifugal_term * pcentrifugal_tmp.block(0, 0, 2, 1));
    pshoulder_[2 * i] = pshoulder_tmp(0, i) + xref(0, 0);
    pshoulder_[2 * i + 1] = pshoulder_tmp(1, i) + xref(1, 0);
  }

  B.setZero();

  if (S(0, 0) == Scalar(1)) {
    B.block(0, 0, 2, 2).setIdentity();
    B.block(6, 2, 2, 2).setIdentity();
  } else {
    B.block(2, 0, 2, 2).setIdentity();
    B.block(4, 2, 2, 2).setIdentity();
  }
}
}  // namespace quadruped_walkgen

#endif
