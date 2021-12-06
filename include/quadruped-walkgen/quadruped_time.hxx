#ifndef __quadruped_walkgen_quadruped_time_hxx__
#define __quadruped_walkgen_quadruped_time_hxx__

#include "crocoddyl/core/utils/exception.hpp"

namespace quadruped_walkgen {
template <typename Scalar>
ActionModelQuadrupedTimeTpl<Scalar>::ActionModelQuadrupedTimeTpl()
    : crocoddyl::ActionModelAbstractTpl<Scalar>(boost::make_shared<crocoddyl::StateVectorTpl<Scalar> >(21), 1, 22) {
  // 4 costs
  state_weights_ << Scalar(1.), Scalar(1.), Scalar(150.), Scalar(35.), Scalar(30.), Scalar(8.), Scalar(20.),
      Scalar(20.), Scalar(15.), Scalar(4.), Scalar(4.), Scalar(8.);
  heuristic_weights_.setConstant(Scalar(1));

  dt_bound_weight_cmd = Scalar(100.);
  dt_weight_cmd = Scalar(0.);

  pheuristic_.setZero();
  gait_double_.setZero();

  // Shoulder heuristic position
  // pshoulder_ <<  Scalar(0.1946) ,  Scalar(0.15005),  Scalar(0.1946) ,  Scalar(-0.15005) ,
  //                Scalar(-0.1946),  Scalar(0.15005) , Scalar(-0.1946),  Scalar(-0.15005) ;
  // pshoulder_0 <<  Scalar(0.1946) ,   Scalar(0.1946) ,   Scalar(-0.1946),  Scalar(-0.1946) ,
  //                 Scalar(0.15005) ,  Scalar(-0.15005)  , Scalar(0.15005)  ,  Scalar(-0.15005) ;
  // pshoulder_tmp.setZero() ;
  // pcentrifugal_tmp_1.setZero() ;
  // pcentrifugal_tmp_2.setZero() ;
  // pcentrifugal_tmp.setZero() ;
  centrifugal_term = true;
  symmetry_term = true;
  T_gait = Scalar(0.64);
  // R_tmp.setZero() ;

  // Cost relative to the command
  rub_max_.setZero();
  rub_max_bool.setZero();
  dt_ref_.setConstant(Scalar(0.02));
  dt_min_.setConstant(Scalar(0.005));
  dt_max_.setConstant(Scalar(0.1));

  // Log cost
  cost_.setZero();
  log_cost = true;
}

template <typename Scalar>
ActionModelQuadrupedTimeTpl<Scalar>::~ActionModelQuadrupedTimeTpl() {}

template <typename Scalar>
void ActionModelQuadrupedTimeTpl<Scalar>::calc(
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

  ActionDataQuadrupedTimeTpl<Scalar>* d = static_cast<ActionDataQuadrupedTimeTpl<Scalar>*>(data.get());

  d->xnext.template head<12>() = x.head(12);
  d->xnext.template segment<8>(12) = x.segment(12, 8);
  d->xnext.template tail<1>() = u.cwiseAbs();

  // Residual cost on the state and force norm
  // State : delta*||X-Xref||
  d->r.template head<12>() = state_weights_.cwiseProduct(x.head(12) - xref_);
  // Feet placement : delta*||P-Pref||
  d->r.template segment<8>(12) =
      ((heuristic_weights_.cwiseProduct(x.segment(12, 8) - pheuristic_)).array() * gait_double_.array()).matrix();
  // Dt command, used to fix the optimisation to dt_ref_
  d->r.template tail<1>() << dt_weight_cmd * (u.cwiseAbs() - dt_ref_);

  // Penalisation if dt out of lower/upper bound
  rub_max_ << dt_min_ - u.cwiseAbs(), u.cwiseAbs() - dt_max_;
  rub_max_bool = (rub_max_.array() >= Scalar(0.)).matrix().template cast<Scalar>();
  rub_max_ = rub_max_.cwiseMax(Scalar(0.));

  d->cost = Scalar(0.5) * d->r.transpose() * d->r + dt_bound_weight_cmd * Scalar(0.5) * rub_max_.squaredNorm();

  if (log_cost) {
    // Length of the cost similar to the Augmented cost...
    cost_[0] = Scalar(0.5) * d->r.head(12).transpose() * d->r.head(12);              // state
    cost_[1] = Scalar(0.5) * d->r.segment(12, 8).transpose() * d->r.segment(12, 8);  // feet placement
    cost_[2] = Scalar(0.5) * d->r.tail(1).transpose() * d->r.tail(1);                // 0.5*\delta*|| U - dt_ref ||^2
    cost_[3] = dt_bound_weight_cmd * Scalar(0.5) * rub_max_.squaredNorm();           // u/l bound weight
    cost_[4] = Scalar(0);
    cost_[5] = Scalar(0);
    cost_[6] = Scalar(0);
  }
}

template <typename Scalar>
void ActionModelQuadrupedTimeTpl<Scalar>::calcDiff(
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

  ActionDataQuadrupedTimeTpl<Scalar>* d = static_cast<ActionDataQuadrupedTimeTpl<Scalar>*>(data.get());

  // Cost derivatives : Lx
  d->Lx.template head<12>() = (state_weights_.array() * d->r.template head<12>().array()).matrix();
  d->Lx.template segment<8>(12) =
      (heuristic_weights_.array() * d->r.template segment<8>(12).array()).matrix();  // * gait_double in d->r

  d->Lu << dt_bound_weight_cmd * std::copysign(1., u(0)) * (-rub_max_[0] + rub_max_[1]);
  d->Lu += dt_weight_cmd * std::copysign(1., u(0)) * d->r.template tail<1>();

  // Hessian : Lxx
  d->Lxx.diagonal().head(12) = (state_weights_.array() * state_weights_.array()).matrix();
  d->Lxx.diagonal().segment(12, 8) =
      (gait_double_.array() * heuristic_weights_.array() * heuristic_weights_.array()).matrix();

  d->Luu.diagonal() << dt_weight_cmd * dt_weight_cmd + dt_bound_weight_cmd * rub_max_bool[0] +
                           dt_bound_weight_cmd * rub_max_bool[1];

  // Dynamic derivatives
  d->Fx.setIdentity();
  d->Fx(20, 20) = Scalar(0.);
  d->Fu.block(20, 0, 1, 1) << std::copysign(1., u(0));
}

template <typename Scalar>
boost::shared_ptr<crocoddyl::ActionDataAbstractTpl<Scalar> > ActionModelQuadrupedTimeTpl<Scalar>::createData() {
  return boost::make_shared<ActionDataQuadrupedTimeTpl<Scalar> >(this);
}

////////////////////////////////
// get & set state weight //////
////////////////////////////////

template <typename Scalar>
const typename Eigen::Matrix<Scalar, 12, 1>& ActionModelQuadrupedTimeTpl<Scalar>::get_state_weights() const {
  return state_weights_;
}
template <typename Scalar>
void ActionModelQuadrupedTimeTpl<Scalar>::set_state_weights(const typename MathBase::VectorXs& weights) {
  if (static_cast<std::size_t>(weights.size()) != 12) {
    throw_pretty("Invalid argument: "
                 << "Weights vector has wrong dimension (it should be 12)");
  }
  state_weights_ = weights;
}

template <typename Scalar>
const typename Eigen::Matrix<Scalar, 8, 1>& ActionModelQuadrupedTimeTpl<Scalar>::get_heuristic_weights() const {
  return heuristic_weights_;
}
template <typename Scalar>
void ActionModelQuadrupedTimeTpl<Scalar>::set_heuristic_weights(const typename MathBase::VectorXs& weights) {
  if (static_cast<std::size_t>(weights.size()) != 8) {
    throw_pretty("Invalid argument: "
                 << "Weights vector has wrong dimension (it should be 8)");
  }
  heuristic_weights_ = weights;
}

//////////////////////////////////////////////////////////
//   Param relative to the heuristic position  //
//////////////////////////////////////////////////////////

template <typename Scalar>
const bool& ActionModelQuadrupedTimeTpl<Scalar>::get_symmetry_term() const {
  return symmetry_term;
}
template <typename Scalar>
void ActionModelQuadrupedTimeTpl<Scalar>::set_symmetry_term(const bool& sym_term) {
  // The model need to be updated after this changed
  symmetry_term = sym_term;
}

template <typename Scalar>
const bool& ActionModelQuadrupedTimeTpl<Scalar>::get_centrifugal_term() const {
  return centrifugal_term;
}
template <typename Scalar>
void ActionModelQuadrupedTimeTpl<Scalar>::set_centrifugal_term(const bool& cent_term) {
  // The model need to be updated after this changed
  centrifugal_term = cent_term;
}

template <typename Scalar>
const Scalar& ActionModelQuadrupedTimeTpl<Scalar>::get_T_gait() const {
  // The model need to be updated after this change
  return T_gait;
}
template <typename Scalar>
void ActionModelQuadrupedTimeTpl<Scalar>::set_T_gait(const Scalar& T_gait_) {
  // The model need to be updated after this change
  T_gait = T_gait_;
}

// dt_ref to fix the dt optimised
template <typename Scalar>
const Scalar& ActionModelQuadrupedTimeTpl<Scalar>::get_dt_ref() const {
  return dt_ref_[0];
}
template <typename Scalar>
void ActionModelQuadrupedTimeTpl<Scalar>::set_dt_ref(const Scalar& dt) {
  dt_ref_[0] = dt;
}

//////////////////////////////////////////////////////////
//   Upper and lower bound dt                           //
//////////////////////////////////////////////////////////
template <typename Scalar>
const Scalar& ActionModelQuadrupedTimeTpl<Scalar>::get_dt_min() const {
  return dt_min_[0];
}
template <typename Scalar>
void ActionModelQuadrupedTimeTpl<Scalar>::set_dt_min(const Scalar& dt) {
  dt_min_[0] = dt;
}
template <typename Scalar>
const Scalar& ActionModelQuadrupedTimeTpl<Scalar>::get_dt_max() const {
  return dt_max_[0];
}
template <typename Scalar>
void ActionModelQuadrupedTimeTpl<Scalar>::set_dt_max(const Scalar& dt) {
  dt_max_[0] = dt;
}
/////////////////////////////////////////////////////////////
// Get access and modify to the upper/lower bound for dt   //
/////////////////////////////////////////////////////////////
template <typename Scalar>
const Scalar& ActionModelQuadrupedTimeTpl<Scalar>::get_dt_bound_weight_cmd() const {
  return dt_bound_weight_cmd;
}
template <typename Scalar>
void ActionModelQuadrupedTimeTpl<Scalar>::set_dt_bound_weight_cmd(const Scalar& weight_) {
  dt_bound_weight_cmd = weight_;
}

///////////////////////////////////////////////////////////////
// Get access and modify to the \weight*||U-dt_ref|| weight  //
///////////////////////////////////////////////////////////////
template <typename Scalar>
const Scalar& ActionModelQuadrupedTimeTpl<Scalar>::get_dt_weight_cmd() const {
  return dt_weight_cmd;
}
template <typename Scalar>
void ActionModelQuadrupedTimeTpl<Scalar>::set_dt_weight_cmd(const Scalar& weight_) {
  dt_weight_cmd = weight_;
}

///////////////////////
// Logging cost      //
///////////////////////
template <typename Scalar>
const typename Eigen::Matrix<Scalar, 7, 1>& ActionModelQuadrupedTimeTpl<Scalar>::get_cost() const {
  return cost_;
}

////////////////////////
// Update current model
////////////////////////

template <typename Scalar>
void ActionModelQuadrupedTimeTpl<Scalar>::update_model(const Eigen::Ref<const typename MathBase::MatrixXs>& l_feet,
                                                       const Eigen::Ref<const typename MathBase::MatrixXs>& xref,
                                                       const Eigen::Ref<const typename MathBase::VectorXs>& S) {
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
  // l_feet and S useless, kept for now, to be consistent with others models

  xref_ = xref;
  for (int i = 0; i < 4; i = i + 1) {
    gait_double_[2 * i] = S[i];
    gait_double_[2 * i + 1] = S[i];
    pheuristic_[2 * i] = l_feet(0, i);
    pheuristic_[2 * i + 1] = l_feet(1, i);
  }

  // To compute heuristic here
  // R_tmp << cos(xref(5,0)) ,-sin(xref(5,0)) , Scalar(0),
  //     sin(xref(5,0)), cos(xref(5,0)), Scalar(0),
  //     Scalar(0),Scalar(0),Scalar(1.0) ;

  //  // Centrifual term
  // pcentrifugal_tmp_1 = xref.block(6,0,3,1) ;
  // pcentrifugal_tmp_2 = xref.block(9,0,3,1) ;
  // pcentrifugal_tmp = 0.5*std::sqrt(xref(2,0)/9.81) * pcentrifugal_tmp_1.cross(pcentrifugal_tmp_2) ;

  // for (int i=0; i<4; i=i+1){
  //   pshoulder_tmp.block(0,i,2,1) =  R_tmp.block(0,0,2,2)*(pshoulder_0.block(0,i,2,1) +   symmetry_term *
  //   0.25*T_gait*xref.block(6,0,2,1) + centrifugal_term * pcentrifugal_tmp.block(0,0,2,1) ); pshoulder_[2*i] =
  //   pshoulder_tmp(0,i) +   xref(0,0); pshoulder_[2*i+1] = pshoulder_tmp(1,i) +  xref(1,0);
  // }
}
}  // namespace quadruped_walkgen

#endif
