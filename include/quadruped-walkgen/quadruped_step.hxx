#ifndef __quadruped_walkgen_quadruped_step_hxx__
#define __quadruped_walkgen_quadruped_step_hxx__

#include "crocoddyl/core/utils/exception.hpp"

namespace quadruped_walkgen {
template <typename Scalar>
ActionModelQuadrupedStepTpl<Scalar>::ActionModelQuadrupedStepTpl()
    : crocoddyl::ActionModelAbstractTpl<Scalar>(boost::make_shared<crocoddyl::StateVectorTpl<Scalar> >(20), 8, 28) {
  B.setZero();
  state_weights_ << Scalar(1.), Scalar(1.), Scalar(150.), Scalar(35.), Scalar(30.), Scalar(8.), Scalar(20.),
      Scalar(20.), Scalar(15.), Scalar(4.), Scalar(4.), Scalar(8.);

  pheuristic_ << Scalar(0.18), Scalar(0.15005), Scalar(0.18), Scalar(-0.15005), Scalar(-0.21), Scalar(0.15005),
      Scalar(-0.21), Scalar(-0.15005);

  centrifugal_term = true;
  symmetry_term = true;
  T_gait = Scalar(0.48);

  step_weights_.setConstant(Scalar(1));
  heuristic_weights_.setConstant(Scalar(1));

  // Compute heuristic inside
  // pshoulder_0 << Scalar(0.1946), Scalar(0.1946), Scalar(-0.1946), Scalar(-0.1946), Scalar(0.15005),
  // Scalar(-0.15005),
  //     Scalar(0.15005), Scalar(-0.15005);
  // pshoulder_tmp.setZero();
  // pcentrifugal_tmp_1.setZero();
  // pcentrifugal_tmp_2.setZero();
  // pcentrifugal_tmp.setZero();

  N_sampling = 5;       // Number of point to sample the polynomial curve of the feet trajectory
  S_.setZero();         // Usefull to compute only the trajectory for moving feet
  position_.setZero();  // Xk+1 = Xk + Uk, Xk does not correspond to the current position of the flying feet, Delta_x
                        // is not straightforward

  // Cost on the acceleration of the feet :
  is_acc_activated_ = true;
  acc_weight_ = Scalar(1.);
  acc_lim_.setConstant(Scalar(50.));

  delta_ = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(N_sampling - 1, 4);
  gamma_ = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(N_sampling - 1, 3);
  for (int k = 1; k < N_sampling; k++) {
    delta_(k - 1, 0) = (float)k / (float)N_sampling;  // [1/N, 2/N, ... , (N-1)/N]
  }
  delta_.col(1) << delta_.col(0).pow(2);
  delta_.col(2) << delta_.col(0).pow(3);
  delta_.col(3) << delta_.col(0).pow(4);  // Only used for speed cost
  gamma_.col(0) = 60 * delta_.col(0) - 180 * delta_.col(1) + 120 * delta_.col(2);
  gamma_.col(1) = -36 * delta_.col(0) + 96 * delta_.col(1) - 60 * delta_.col(2);
  gamma_.col(2) = -9 * delta_.col(0) + 18 * delta_.col(1) - 10 * delta_.col(2);

  alpha_ = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(N_sampling - 1, 1);  // Common for 4 feet
  beta_x_ =
      Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(N_sampling - 1, 4);  // Depends on a0_x, v0_x of feet
  beta_y_ =
      Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(N_sampling - 1, 4);  // Depends on a0_y, v0_y of feet
  tmp_ones_ = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Ones(N_sampling - 1, 1);

  rb_accx_max_ = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(N_sampling - 1, 8);
  rb_accy_max_ = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(N_sampling - 1, 8);
  rb_accx_max_bool_ = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(N_sampling - 1, 8);
  rb_accy_max_bool_ = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(N_sampling - 1, 8);

  // Cost on the velocity of the feet :
  is_vel_activated_ = true;
  vel_weight_ = Scalar(1.);
  vel_lim_.setConstant(Scalar(3.));

  gamma_v = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(N_sampling - 1, 4);
  gamma_v.col(0) = 30 * delta_.col(1) - 60 * delta_.col(2) + 30 * delta_.col(3);
  gamma_v.col(1) = delta_.col(0);
  gamma_v.col(2) = -18 * delta_.col(1) + 32 * delta_.col(2) - 15 * delta_.col(3);
  gamma_v.col(3) = -4.5 * delta_.col(1) + 6 * delta_.col(2) - 2.5 * delta_.col(3);

  alpha_v = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(N_sampling - 1, 1);  // Common for 4 feet
  beta_x_v =
      Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(N_sampling - 1, 4);  // Depends on a0_x, v0_x of feet
  beta_y_v =
      Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(N_sampling - 1, 4);  // Depends on a0_y, v0_y of feet

  rb_velx_max_ = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(N_sampling - 1, 8);
  rb_vely_max_ = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(N_sampling - 1, 8);
  rb_velx_max_bool_ = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(N_sampling - 1, 8);
  rb_vely_max_bool_ = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(N_sampling - 1, 8);

  // Cost on the jerk at t=0
  is_jerk_activated_ = true;
  jerk_weight_ = Scalar(1.);
  alpha_j = Scalar(0.);                                                       // Common for 4 feet
  beta_j = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(2, 4);  // Depends on a0_x, v0_x of feet
  rb_jerk_.setZero();
}

template <typename Scalar>
ActionModelQuadrupedStepTpl<Scalar>::~ActionModelQuadrupedStepTpl() {}

template <typename Scalar>
void ActionModelQuadrupedStepTpl<Scalar>::calc(
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

  ActionDataQuadrupedStepTpl<Scalar>* d = static_cast<ActionDataQuadrupedStepTpl<Scalar>*>(data.get());

  d->xnext.template head<12>() = x.head(12);
  d->xnext.template tail<8>() = x.tail(8) + B * u;

  // Residual cost on the state and force norm
  d->r.template head<12>() = state_weights_.cwiseProduct(x.head(12) - xref_);
  d->r.template segment<8>(12) = heuristic_weights_.cwiseProduct(x.tail(8) - pheuristic_);
  d->r.template tail<8>() = step_weights_.cwiseProduct(u);

  d->cost = Scalar(0.5) * d->r.transpose() * d->r;

  // Weight on the feet acceleration :
  if (is_acc_activated_) {
    for (int i = 0; i < 4; i++) {
      if (S_(i) == Scalar(1.)) {
        // position_ expressed in base frame
        // rb_accx_max_.col(2*i) = (x(12+ 2*i) + u(2*i) - position_(0,i))*alpha_ + beta_x_.col(i) -
        // acc_lim_(0)*tmp_ones_; rb_accx_max_.col(2*i+1) = -( x(12+ 2*i) + u(2*i) - position_(0,i) )*alpha_ -
        // beta_x_.col(i) - acc_lim_(0)*tmp_ones_;

        // rb_accy_max_.col(2*i) = ( x(12+ 2*i+1) + u(2*i+1) - position_(1,i) )*alpha_ + beta_y_.col(i) -
        // acc_lim_(1)*tmp_ones_; rb_accy_max_.col(2*i+1) = -( x(12+ 2*i+1) + u(2*i+1) - position_(1,i) )*alpha_ -
        // beta_y_.col(i) - acc_lim_(1)*tmp_ones_;

        // position_ expressed in world frame
        rb_accx_max_.col(2 * i) = (oRh_(0, 0) * (x(12 + 2 * i) + u(2 * i)) +
                                   oRh_(0, 1) * (x(12 + 2 * i + 1) + u(2 * i + 1)) + oTh_(0) - position_(0, i)) *
                                      alpha_ +
                                  beta_x_.col(i) - acc_lim_(0) * tmp_ones_;
        rb_accx_max_.col(2 * i + 1) = -(oRh_(0, 0) * (x(12 + 2 * i) + u(2 * i)) +
                                        oRh_(0, 1) * (x(12 + 2 * i + 1) + u(2 * i + 1)) + oTh_(0) - position_(0, i)) *
                                          alpha_ +
                                      beta_x_.col(i) - acc_lim_(0) * tmp_ones_;
        ;

        rb_accy_max_.col(2 * i) = (oRh_(1, 0) * (x(12 + 2 * i) + u(2 * i)) +
                                   oRh_(1, 1) * (x(12 + 2 * i + 1) + u(2 * i + 1)) + oTh_(1) - position_(1, i)) *
                                      alpha_ +
                                  beta_y_.col(i) - acc_lim_(1) * tmp_ones_;
        rb_accy_max_.col(2 * i + 1) = -(oRh_(1, 0) * (x(12 + 2 * i) + u(2 * i)) +
                                        oRh_(1, 1) * (x(12 + 2 * i + 1) + u(2 * i + 1)) + oTh_(1) - position_(1, i)) *
                                          alpha_ -
                                      beta_y_.col(i) - acc_lim_(1) * tmp_ones_;
      } else {
        rb_accx_max_.col(2 * i).setZero();
        rb_accx_max_.col(2 * i + 1).setZero();
        rb_accy_max_.col(2 * i).setZero();
        rb_accy_max_.col(2 * i + 1).setZero();
      }
    }
    rb_accx_max_bool_ = (rb_accx_max_ > Scalar(0.)).template cast<Scalar>();  // Usefull to compute the derivatives
    rb_accy_max_bool_ = (rb_accy_max_ > Scalar(0.)).template cast<Scalar>();  // Usefull to compute the derivatives

    rb_accx_max_ = rb_accx_max_.cwiseMax(Scalar(0.));
    rb_accy_max_ = rb_accy_max_.cwiseMax(Scalar(0.));

    for (int foot = 0; foot < 4; foot++) {
      if (S_(foot) == Scalar(1.)) {
        for (int i = 0; i < (N_sampling - 1); i++) {
          if (rb_accx_max_bool_(i, 2 * foot)) {
            d->cost += Scalar(0.5) * acc_weight_ * pow(rb_accx_max_(i, 2 * foot), 2);
          }
          if (rb_accx_max_bool_(i, 2 * foot + 1)) {
            d->cost += Scalar(0.5) * acc_weight_ * pow(rb_accx_max_(i, 2 * foot + 1), 2);
          }
          if (rb_accy_max_bool_(i, 2 * foot)) {
            d->cost += Scalar(0.5) * acc_weight_ * pow(rb_accy_max_(i, 2 * foot), 2);
          }
          if (rb_accy_max_bool_(i, 2 * foot + 1)) {
            d->cost += Scalar(0.5) * acc_weight_ * pow(rb_accy_max_(i, 2 * foot + 1), 2);
          }
        }
      }
    }
  }

  // Weight on the feet velocity
  if (is_vel_activated_) {
    for (int i = 0; i < 4; i++) {
      if (S_(i) == Scalar(1.)) {
        // position_ expressed in local frame
        // rb_velx_max_.col(2*i) = (x(12+ 2*i) + u(2*i) - position_(0,i))*alpha_v + beta_x_v.col(i) -
        // vel_lim_(0)*tmp_ones_; rb_velx_max_.col(2*i+1) = -( x(12+ 2*i) + u(2*i) - position_(0,i) )*alpha_v -
        // beta_x_v.col(i) - vel_lim_(0)*tmp_ones_;

        // rb_vely_max_.col(2*i) = ( x(12+ 2*i+1) + u(2*i+1) - position_(1,i) )*alpha_v + beta_y_v.col(i) -
        // vel_lim_(1)*tmp_ones_; rb_vely_max_.col(2*i+1) = -( x(12+ 2*i+1) + u(2*i+1) - position_(1,i) )*alpha_v -
        // beta_y_v.col(i) - vel_lim_(1)*tmp_ones_;

        // position_ expressed in world frame
        rb_velx_max_.col(2 * i) = (oRh_(0, 0) * (x(12 + 2 * i) + u(2 * i)) +
                                   oRh_(0, 1) * (x(12 + 2 * i + 1) + u(2 * i + 1)) + oTh_(0) - position_(0, i)) *
                                      alpha_v +
                                  beta_x_v.col(i) - vel_lim_(0) * tmp_ones_;
        rb_velx_max_.col(2 * i + 1) = -(oRh_(0, 0) * (x(12 + 2 * i) + u(2 * i)) +
                                        oRh_(0, 1) * (x(12 + 2 * i + 1) + u(2 * i + 1)) + oTh_(0) - position_(0, i)) *
                                          alpha_v -
                                      beta_x_v.col(i) - vel_lim_(0) * tmp_ones_;

        rb_vely_max_.col(2 * i) = (oRh_(1, 0) * (x(12 + 2 * i) + u(2 * i)) +
                                   oRh_(1, 1) * (x(12 + 2 * i + 1) + u(2 * i + 1)) + oTh_(1) - position_(1, i)) *
                                      alpha_v +
                                  beta_y_v.col(i) - vel_lim_(1) * tmp_ones_;
        rb_vely_max_.col(2 * i + 1) = -(oRh_(1, 0) * (x(12 + 2 * i) + u(2 * i)) +
                                        oRh_(1, 1) * (x(12 + 2 * i + 1) + u(2 * i + 1)) + oTh_(1) - position_(1, i)) *
                                          alpha_v -
                                      beta_y_v.col(i) - vel_lim_(1) * tmp_ones_;
      } else {
        rb_velx_max_.col(2 * i).setZero();
        rb_velx_max_.col(2 * i + 1).setZero();
        rb_vely_max_.col(2 * i).setZero();
        rb_vely_max_.col(2 * i + 1).setZero();
      }
    }
    rb_velx_max_bool_ = (rb_velx_max_ > Scalar(0.)).template cast<Scalar>();  // Usefull to compute the derivatives
    rb_vely_max_bool_ = (rb_vely_max_ > Scalar(0.)).template cast<Scalar>();  // Usefull to compute the derivatives

    rb_velx_max_ = rb_velx_max_.cwiseMax(Scalar(0.));
    rb_vely_max_ = rb_vely_max_.cwiseMax(Scalar(0.));

    for (int foot = 0; foot < 4; foot++) {
      if (S_(foot) == Scalar(1.)) {
        for (int i = 0; i < (N_sampling - 1); i++) {
          if (rb_velx_max_bool_(i, 2 * foot)) {
            d->cost += Scalar(0.5) * vel_weight_ * pow(rb_velx_max_(i, 2 * foot), 2);
          }
          if (rb_velx_max_bool_(i, 2 * foot + 1)) {
            d->cost += Scalar(0.5) * vel_weight_ * pow(rb_velx_max_(i, 2 * foot + 1), 2);
          }
          if (rb_vely_max_bool_(i, 2 * foot)) {
            d->cost += Scalar(0.5) * vel_weight_ * pow(rb_vely_max_(i, 2 * foot), 2);
          }
          if (rb_vely_max_bool_(i, 2 * foot + 1)) {
            d->cost += Scalar(0.5) * vel_weight_ * pow(rb_vely_max_(i, 2 * foot + 1), 2);
          }
        }
      }
    }
  }

  // Weight on the feet velocity
  if (is_jerk_activated_) {
    for (int i = 0; i < 4; i++) {
      if (S_(i) == Scalar(1.)) {
        rb_jerk_(0, i) = (oRh_(0, 0) * (x(12 + 2 * i) + u(2 * i)) + oRh_(0, 1) * (x(12 + 2 * i + 1) + u(2 * i + 1)) +
                          oTh_(0) - position_(0, i)) *
                             alpha_j +
                         beta_j(0, i) - jerk_(0, i);
        rb_jerk_(1, i) = (oRh_(1, 0) * (x(12 + 2 * i) + u(2 * i)) + oRh_(1, 1) * (x(12 + 2 * i + 1) + u(2 * i + 1)) +
                          oTh_(1) - position_(1, i)) *
                             alpha_j +
                         beta_j(1, i) - jerk_(1, i);
        d->cost += Scalar(0.5) * jerk_weight_ * rb_jerk_.col(i).squaredNorm();
      } else {
        rb_jerk_.col(i).setZero();
      }
    }
  }
}

template <typename Scalar>
void ActionModelQuadrupedStepTpl<Scalar>::calcDiff(
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

  ActionDataQuadrupedStepTpl<Scalar>* d = static_cast<ActionDataQuadrupedStepTpl<Scalar>*>(data.get());

  d->Lxu.setZero();
  d->Luu.setZero();

  // Cost derivatives : Lx
  d->Lx.template head<12>() = (state_weights_.array() * d->r.template head<12>().array()).matrix();
  d->Lx.template tail<8>() = (heuristic_weights_.array() * d->r.template segment<8>(12).array()).matrix();

  d->Lu = (step_weights_.array() * d->r.template tail<8>().array()).matrix();

  // Hessian : Lxx
  d->Lxx.diagonal().head(12) = (state_weights_.array() * state_weights_.array()).matrix();
  d->Lxx.diagonal().tail(8) = (heuristic_weights_.array() * heuristic_weights_.array()).matrix();

  d->Luu.diagonal() = (step_weights_.array() * step_weights_.array()).matrix();

  if (is_acc_activated_) {
    for (int foot = 0; foot < 4; foot++) {
      if (S_[foot] == Scalar(1)) {
        for (int i = 0; i < (N_sampling - 1); i++) {
          // Position_ expressed in local frame
          // if (rb_accx_max_bool_(i,2*foot)){

          //   d->Lu(2*foot) += acc_weight_ * alpha_(i) * rb_accx_max_(i,2*foot);
          //   d->Luu(2*foot,2*foot) += acc_weight_ * pow(alpha_(i),2);

          //   d->Lx(12+2*foot) += acc_weight_ * alpha_(i) * rb_accx_max_(i,2*foot);
          //   d->Lxu(12+2*foot,2*foot) += acc_weight_ * pow(alpha_(i),2);
          //   d->Lxx(12+2*foot,12+2*foot) += acc_weight_ * pow(alpha_(i),2);
          // }
          // if (rb_accx_max_bool_(i,2*foot+1)){
          //   d->Lu(2*foot) += - acc_weight_ * alpha_(i) * rb_accx_max_(i,2*foot+1);
          //   d->Luu(2*foot,2*foot) += acc_weight_ * pow(alpha_(i),2);

          //   d->Lx(12+2*foot) += - acc_weight_ * alpha_(i) * rb_accx_max_(i,2*foot+1);
          //   d->Lxu(12+2*foot,2*foot) += acc_weight_ * pow(alpha_(i),2);
          //   d->Lxx(12+2*foot,12+2*foot) += acc_weight_ * pow(alpha_(i),2);
          // }
          // if (rb_accy_max_bool_(i,2*foot)){
          //   d->Lu(2*foot+1) += acc_weight_ * alpha_(i) * rb_accy_max_(i,2*foot);
          //   d->Luu(2*foot+1,2*foot+1) += acc_weight_ * pow(alpha_(i),2);

          //   d->Lx(12+2*foot+1) += acc_weight_ * alpha_(i) * rb_accy_max_(i,2*foot);
          //   d->Lxu(12+2*foot+1,2*foot+1) += acc_weight_ * pow(alpha_(i),2);
          //   d->Lxx(12+2*foot+1,12+2*foot+1) += acc_weight_ * pow(alpha_(i),2);
          // }
          // if (rb_accy_max_bool_(i,2*foot+1)){
          //   d->Lu(2*foot+1) += - acc_weight_ * alpha_(i) * rb_accy_max_(i,2*foot+1);
          //   d->Luu(2*foot+1,2*foot+1) += acc_weight_ * pow(alpha_(i),2);

          //   d->Lx(12+2*foot+1) += - acc_weight_ * alpha_(i) * rb_accy_max_(i,2*foot+1);
          //   d->Lxu(12+2*foot+1,2*foot+1) += acc_weight_ * pow(alpha_(i),2);
          //   d->Lxx(12+2*foot+1,12+2*foot+1) += acc_weight_ * pow(alpha_(i),2);
          // }

          // Position_ expressed in world frame
          if (rb_accx_max_bool_(i, 2 * foot)) {
            d->Lu(2 * foot) += acc_weight_ * oRh_(0, 0) * alpha_(i) * rb_accx_max_(i, 2 * foot);
            d->Lu(2 * foot + 1) += acc_weight_ * oRh_(0, 1) * alpha_(i) * rb_accx_max_(i, 2 * foot);

            d->Luu(2 * foot, 2 * foot) += acc_weight_ * pow(oRh_(0, 0), 2) * pow(alpha_(i), 2);
            d->Luu(2 * foot + 1, 2 * foot + 1) += acc_weight_ * pow(oRh_(0, 1), 2) * pow(alpha_(i), 2);

            d->Luu(2 * foot, 2 * foot + 1) += acc_weight_ * oRh_(0, 0) * oRh_(0, 1) * pow(alpha_(i), 2);
            d->Luu(2 * foot + 1, 2 * foot) += acc_weight_ * oRh_(0, 0) * oRh_(0, 1) * pow(alpha_(i), 2);

            d->Lx(12 + 2 * foot) += acc_weight_ * oRh_(0, 0) * alpha_(i) * rb_accx_max_(i, 2 * foot);
            d->Lx(12 + 2 * foot + 1) += acc_weight_ * oRh_(0, 1) * alpha_(i) * rb_accx_max_(i, 2 * foot);

            d->Lxx(12 + 2 * foot, 12 + 2 * foot) += acc_weight_ * pow(oRh_(0, 0), 2) * pow(alpha_(i), 2);
            d->Lxx(12 + 2 * foot + 1, 12 + 2 * foot + 1) += acc_weight_ * pow(oRh_(0, 1), 2) * pow(alpha_(i), 2);

            d->Lxx(12 + 2 * foot, 12 + 2 * foot + 1) += acc_weight_ * oRh_(0, 0) * oRh_(0, 1) * pow(alpha_(i), 2);
            d->Lxx(12 + 2 * foot + 1, 12 + 2 * foot) += acc_weight_ * oRh_(0, 0) * oRh_(0, 1) * pow(alpha_(i), 2);

            d->Lxu(12 + 2 * foot, 2 * foot) += acc_weight_ * pow(alpha_(i), 2) * pow(oRh_(0, 0), 2);
            d->Lxu(12 + 2 * foot, 2 * foot + 1) += acc_weight_ * pow(alpha_(i), 2) * oRh_(0, 0) * oRh_(0, 1);
            d->Lxu(12 + 2 * foot + 1, 2 * foot) += acc_weight_ * pow(alpha_(i), 2) * oRh_(0, 0) * oRh_(0, 1);
            d->Lxu(12 + 2 * foot + 1, 2 * foot + 1) += acc_weight_ * pow(alpha_(i), 2) * pow(oRh_(0, 1), 2);
          }
          if (rb_accx_max_bool_(i, 2 * foot + 1)) {
            d->Lu(2 * foot) += -acc_weight_ * oRh_(0, 0) * alpha_(i) * rb_accx_max_(i, 2 * foot + 1);
            d->Lu(2 * foot + 1) += -acc_weight_ * oRh_(0, 1) * alpha_(i) * rb_accx_max_(i, 2 * foot + 1);

            d->Luu(2 * foot, 2 * foot) += acc_weight_ * pow(oRh_(0, 0), 2) * pow(alpha_(i), 2);
            d->Luu(2 * foot + 1, 2 * foot + 1) += acc_weight_ * pow(oRh_(0, 1), 2) * pow(alpha_(i), 2);

            d->Luu(2 * foot, 2 * foot + 1) += acc_weight_ * oRh_(0, 0) * oRh_(0, 1) * pow(alpha_(i), 2);
            d->Luu(2 * foot + 1, 2 * foot) += acc_weight_ * oRh_(0, 0) * oRh_(0, 1) * pow(alpha_(i), 2);

            d->Lx(12 + 2 * foot) += -acc_weight_ * oRh_(0, 0) * alpha_(i) * rb_accx_max_(i, 2 * foot + 1);
            d->Lx(12 + 2 * foot + 1) += -acc_weight_ * oRh_(0, 1) * alpha_(i) * rb_accx_max_(i, 2 * foot + 1);

            d->Lxx(12 + 2 * foot, 12 + 2 * foot) += acc_weight_ * pow(oRh_(0, 0), 2) * pow(alpha_(i), 2);
            d->Lxx(12 + 2 * foot + 1, 12 + 2 * foot + 1) += acc_weight_ * pow(oRh_(0, 1), 2) * pow(alpha_(i), 2);

            d->Lxx(12 + 2 * foot, 12 + 2 * foot + 1) += acc_weight_ * oRh_(0, 0) * oRh_(0, 1) * pow(alpha_(i), 2);
            d->Lxx(12 + 2 * foot + 1, 12 + 2 * foot) += acc_weight_ * oRh_(0, 0) * oRh_(0, 1) * pow(alpha_(i), 2);

            d->Lxu(12 + 2 * foot, 2 * foot) += acc_weight_ * pow(alpha_(i), 2) * pow(oRh_(0, 0), 2);
            d->Lxu(12 + 2 * foot, 2 * foot + 1) += acc_weight_ * pow(alpha_(i), 2) * oRh_(0, 0) * oRh_(0, 1);
            d->Lxu(12 + 2 * foot + 1, 2 * foot) += acc_weight_ * pow(alpha_(i), 2) * oRh_(0, 0) * oRh_(0, 1);
            d->Lxu(12 + 2 * foot + 1, 2 * foot + 1) += acc_weight_ * pow(alpha_(i), 2) * pow(oRh_(0, 1), 2);
          }
          if (rb_accy_max_bool_(i, 2 * foot)) {
            d->Lu(2 * foot) += acc_weight_ * oRh_(1, 0) * alpha_(i) * rb_accy_max_(i, 2 * foot);
            d->Lu(2 * foot + 1) += acc_weight_ * oRh_(1, 1) * alpha_(i) * rb_accy_max_(i, 2 * foot);

            d->Luu(2 * foot, 2 * foot) += acc_weight_ * pow(oRh_(1, 0), 2) * pow(alpha_(i), 2);
            d->Luu(2 * foot + 1, 2 * foot + 1) += acc_weight_ * pow(oRh_(1, 1), 2) * pow(alpha_(i), 2);

            d->Luu(2 * foot, 2 * foot + 1) += acc_weight_ * oRh_(1, 0) * oRh_(1, 1) * pow(alpha_(i), 2);
            d->Luu(2 * foot + 1, 2 * foot) += acc_weight_ * oRh_(1, 0) * oRh_(1, 1) * pow(alpha_(i), 2);

            d->Lx(12 + 2 * foot) += acc_weight_ * oRh_(1, 0) * alpha_(i) * rb_accy_max_(i, 2 * foot);
            d->Lx(12 + 2 * foot + 1) += acc_weight_ * oRh_(1, 1) * alpha_(i) * rb_accy_max_(i, 2 * foot);

            d->Lxx(12 + 2 * foot, 12 + 2 * foot) += acc_weight_ * pow(oRh_(1, 0), 2) * pow(alpha_(i), 2);
            d->Lxx(12 + 2 * foot + 1, 12 + 2 * foot + 1) += acc_weight_ * pow(oRh_(1, 1), 2) * pow(alpha_(i), 2);

            d->Lxx(12 + 2 * foot, 12 + 2 * foot + 1) += acc_weight_ * oRh_(1, 0) * oRh_(1, 1) * pow(alpha_(i), 2);
            d->Lxx(12 + 2 * foot + 1, 12 + 2 * foot) += acc_weight_ * oRh_(1, 0) * oRh_(1, 1) * pow(alpha_(i), 2);

            d->Lxu(12 + 2 * foot, 2 * foot) += acc_weight_ * pow(alpha_(i), 2) * pow(oRh_(1, 0), 2);
            d->Lxu(12 + 2 * foot, 2 * foot + 1) += acc_weight_ * pow(alpha_(i), 2) * oRh_(1, 0) * oRh_(1, 1);
            d->Lxu(12 + 2 * foot + 1, 2 * foot) += acc_weight_ * pow(alpha_(i), 2) * oRh_(1, 0) * oRh_(1, 1);
            d->Lxu(12 + 2 * foot + 1, 2 * foot + 1) += acc_weight_ * pow(alpha_(i), 2) * pow(oRh_(1, 1), 2);
          }
          if (rb_accy_max_bool_(i, 2 * foot + 1)) {
            d->Lu(2 * foot) += -acc_weight_ * oRh_(1, 0) * alpha_(i) * rb_accy_max_(i, 2 * foot + 1);
            d->Lu(2 * foot + 1) += -acc_weight_ * oRh_(1, 1) * alpha_(i) * rb_accy_max_(i, 2 * foot + 1);

            d->Luu(2 * foot, 2 * foot) += acc_weight_ * pow(oRh_(1, 0), 2) * pow(alpha_(i), 2);
            d->Luu(2 * foot + 1, 2 * foot + 1) += acc_weight_ * pow(oRh_(1, 1), 2) * pow(alpha_(i), 2);

            d->Luu(2 * foot, 2 * foot + 1) += acc_weight_ * oRh_(1, 0) * oRh_(1, 1) * pow(alpha_(i), 2);
            d->Luu(2 * foot + 1, 2 * foot) += acc_weight_ * oRh_(1, 0) * oRh_(1, 1) * pow(alpha_(i), 2);

            d->Lx(12 + 2 * foot) += -acc_weight_ * oRh_(1, 0) * alpha_(i) * rb_accy_max_(i, 2 * foot + 1);
            d->Lx(12 + 2 * foot + 1) += -acc_weight_ * oRh_(1, 1) * alpha_(i) * rb_accy_max_(i, 2 * foot + 1);

            d->Lxx(12 + 2 * foot, 12 + 2 * foot) += acc_weight_ * pow(oRh_(1, 0), 2) * pow(alpha_(i), 2);
            d->Lxx(12 + 2 * foot + 1, 12 + 2 * foot + 1) += acc_weight_ * pow(oRh_(1, 1), 2) * pow(alpha_(i), 2);

            d->Lxx(12 + 2 * foot, 12 + 2 * foot + 1) += acc_weight_ * oRh_(1, 0) * oRh_(1, 1) * pow(alpha_(i), 2);
            d->Lxx(12 + 2 * foot + 1, 12 + 2 * foot) += acc_weight_ * oRh_(1, 0) * oRh_(1, 1) * pow(alpha_(i), 2);

            d->Lxu(12 + 2 * foot, 2 * foot) += acc_weight_ * pow(alpha_(i), 2) * pow(oRh_(1, 0), 2);
            d->Lxu(12 + 2 * foot, 2 * foot + 1) += acc_weight_ * pow(alpha_(i), 2) * oRh_(1, 0) * oRh_(1, 1);
            d->Lxu(12 + 2 * foot + 1, 2 * foot) += acc_weight_ * pow(alpha_(i), 2) * oRh_(1, 0) * oRh_(1, 1);
            d->Lxu(12 + 2 * foot + 1, 2 * foot + 1) += acc_weight_ * pow(alpha_(i), 2) * pow(oRh_(1, 1), 2);
          }
        }
      }
    }
  }

  if (is_vel_activated_) {
    for (int foot = 0; foot < 4; foot++) {
      if (S_[foot] == Scalar(1)) {
        for (int i = 0; i < (N_sampling - 1); i++) {
          // position_ in base frame
          // if (rb_velx_max_bool_(i,2*foot)){

          //   d->Lu(2*foot) += vel_weight_ * alpha_v(i) * rb_velx_max_(i,2*foot);
          //   d->Luu(2*foot,2*foot) += vel_weight_ * pow(alpha_v(i),2);

          //   d->Lx(12+2*foot) += vel_weight_ * alpha_v(i) * rb_velx_max_(i,2*foot);
          //   d->Lxu(12+2*foot,2*foot) += vel_weight_ * pow(alpha_v(i),2);
          //   d->Lxx(12+2*foot,12+2*foot) += vel_weight_ * pow(alpha_v(i),2);
          // }
          // if (rb_velx_max_bool_(i,2*foot+1)){
          //   d->Lu(2*foot) += - vel_weight_ * alpha_v(i) * rb_velx_max_(i,2*foot+1);
          //   d->Luu(2*foot,2*foot) += vel_weight_ * pow(alpha_v(i),2);

          //   d->Lx(12+2*foot) += - vel_weight_ * alpha_v(i) * rb_velx_max_(i,2*foot+1);
          //   d->Lxu(12+2*foot,2*foot) += vel_weight_ * pow(alpha_v(i),2);
          //   d->Lxx(12+2*foot,12+2*foot) += vel_weight_ * pow(alpha_v(i),2);
          // }
          // if (rb_vely_max_bool_(i,2*foot)){
          //   d->Lu(2*foot+1) += vel_weight_ * alpha_v(i) * rb_vely_max_(i,2*foot);
          //   d->Luu(2*foot+1,2*foot+1) += vel_weight_ * pow(alpha_v(i),2);

          //   d->Lx(12+2*foot+1) += vel_weight_ * alpha_v(i) * rb_vely_max_(i,2*foot);
          //   d->Lxu(12+2*foot+1,2*foot+1) += vel_weight_ * pow(alpha_v(i),2);
          //   d->Lxx(12+2*foot+1,12+2*foot+1) += vel_weight_ * pow(alpha_v(i),2);
          // }
          // if (rb_vely_max_bool_(i,2*foot+1)){
          //   d->Lu(2*foot+1) += - vel_weight_ * alpha_v(i) * rb_vely_max_(i,2*foot+1);
          //   d->Luu(2*foot+1,2*foot+1) += vel_weight_ * pow(alpha_v(i),2);

          //   d->Lx(12+2*foot+1) += - vel_weight_ * alpha_v(i) * rb_vely_max_(i,2*foot+1);
          //   d->Lxu(12+2*foot+1,2*foot+1) += vel_weight_ * pow(alpha_v(i),2);
          //   d->Lxx(12+2*foot+1,12+2*foot+1) += vel_weight_ * pow(alpha_v(i),2);
          // }

          // position_ in world frame
          if (rb_velx_max_bool_(i, 2 * foot)) {
            d->Lu(2 * foot) += vel_weight_ * oRh_(0, 0) * alpha_v(i) * rb_velx_max_(i, 2 * foot);
            d->Lu(2 * foot + 1) += vel_weight_ * oRh_(0, 1) * alpha_v(i) * rb_velx_max_(i, 2 * foot);

            d->Luu(2 * foot, 2 * foot) += vel_weight_ * pow(oRh_(0, 0), 2) * pow(alpha_v(i), 2);
            d->Luu(2 * foot + 1, 2 * foot + 1) += vel_weight_ * pow(oRh_(0, 1), 2) * pow(alpha_v(i), 2);

            d->Luu(2 * foot, 2 * foot + 1) += vel_weight_ * oRh_(0, 0) * oRh_(0, 1) * pow(alpha_v(i), 2);
            d->Luu(2 * foot + 1, 2 * foot) += vel_weight_ * oRh_(0, 0) * oRh_(0, 1) * pow(alpha_v(i), 2);

            d->Lx(12 + 2 * foot) += vel_weight_ * oRh_(0, 0) * alpha_v(i) * rb_velx_max_(i, 2 * foot);
            d->Lx(12 + 2 * foot + 1) += vel_weight_ * oRh_(0, 1) * alpha_v(i) * rb_velx_max_(i, 2 * foot);

            d->Lxx(12 + 2 * foot, 12 + 2 * foot) += vel_weight_ * pow(oRh_(0, 0), 2) * pow(alpha_v(i), 2);
            d->Lxx(12 + 2 * foot + 1, 12 + 2 * foot + 1) += vel_weight_ * pow(oRh_(0, 1), 2) * pow(alpha_v(i), 2);

            d->Lxx(12 + 2 * foot, 12 + 2 * foot + 1) += vel_weight_ * oRh_(0, 0) * oRh_(0, 1) * pow(alpha_v(i), 2);
            d->Lxx(12 + 2 * foot + 1, 12 + 2 * foot) += vel_weight_ * oRh_(0, 0) * oRh_(0, 1) * pow(alpha_v(i), 2);

            d->Lxu(12 + 2 * foot, 2 * foot) += vel_weight_ * pow(alpha_v(i), 2) * pow(oRh_(0, 0), 2);
            d->Lxu(12 + 2 * foot, 2 * foot + 1) += vel_weight_ * pow(alpha_v(i), 2) * oRh_(0, 0) * oRh_(0, 1);
            d->Lxu(12 + 2 * foot + 1, 2 * foot) += vel_weight_ * pow(alpha_v(i), 2) * oRh_(0, 0) * oRh_(0, 1);
            d->Lxu(12 + 2 * foot + 1, 2 * foot + 1) += vel_weight_ * pow(alpha_v(i), 2) * pow(oRh_(0, 1), 2);
          }
          if (rb_velx_max_bool_(i, 2 * foot + 1)) {
            d->Lu(2 * foot) += -vel_weight_ * oRh_(0, 0) * alpha_v(i) * rb_velx_max_(i, 2 * foot + 1);
            d->Lu(2 * foot + 1) += -vel_weight_ * oRh_(0, 1) * alpha_v(i) * rb_velx_max_(i, 2 * foot + 1);

            d->Luu(2 * foot, 2 * foot) += vel_weight_ * pow(oRh_(0, 0), 2) * pow(alpha_v(i), 2);
            d->Luu(2 * foot + 1, 2 * foot + 1) += vel_weight_ * pow(oRh_(0, 1), 2) * pow(alpha_v(i), 2);

            d->Luu(2 * foot, 2 * foot + 1) += vel_weight_ * oRh_(0, 0) * oRh_(0, 1) * pow(alpha_v(i), 2);
            d->Luu(2 * foot + 1, 2 * foot) += vel_weight_ * oRh_(0, 0) * oRh_(0, 1) * pow(alpha_v(i), 2);

            d->Lx(12 + 2 * foot) += -vel_weight_ * oRh_(0, 0) * alpha_v(i) * rb_velx_max_(i, 2 * foot + 1);
            d->Lx(12 + 2 * foot + 1) += -vel_weight_ * oRh_(0, 1) * alpha_v(i) * rb_velx_max_(i, 2 * foot + 1);

            d->Lxx(12 + 2 * foot, 12 + 2 * foot) += vel_weight_ * pow(oRh_(0, 0), 2) * pow(alpha_v(i), 2);
            d->Lxx(12 + 2 * foot + 1, 12 + 2 * foot + 1) += vel_weight_ * pow(oRh_(0, 1), 2) * pow(alpha_v(i), 2);

            d->Lxx(12 + 2 * foot, 12 + 2 * foot + 1) += vel_weight_ * oRh_(0, 0) * oRh_(0, 1) * pow(alpha_v(i), 2);
            d->Lxx(12 + 2 * foot + 1, 12 + 2 * foot) += vel_weight_ * oRh_(0, 0) * oRh_(0, 1) * pow(alpha_v(i), 2);

            d->Lxu(12 + 2 * foot, 2 * foot) += vel_weight_ * pow(alpha_v(i), 2) * pow(oRh_(0, 0), 2);
            d->Lxu(12 + 2 * foot, 2 * foot + 1) += vel_weight_ * pow(alpha_v(i), 2) * oRh_(0, 0) * oRh_(0, 1);
            d->Lxu(12 + 2 * foot + 1, 2 * foot) += vel_weight_ * pow(alpha_v(i), 2) * oRh_(0, 0) * oRh_(0, 1);
            d->Lxu(12 + 2 * foot + 1, 2 * foot + 1) += vel_weight_ * pow(alpha_v(i), 2) * pow(oRh_(0, 1), 2);
          }
          if (rb_vely_max_bool_(i, 2 * foot)) {
            d->Lu(2 * foot) += vel_weight_ * oRh_(1, 0) * alpha_v(i) * rb_vely_max_(i, 2 * foot);
            d->Lu(2 * foot + 1) += vel_weight_ * oRh_(1, 1) * alpha_v(i) * rb_vely_max_(i, 2 * foot);

            d->Luu(2 * foot, 2 * foot) += vel_weight_ * pow(oRh_(1, 0), 2) * pow(alpha_v(i), 2);
            d->Luu(2 * foot + 1, 2 * foot + 1) += vel_weight_ * pow(oRh_(1, 1), 2) * pow(alpha_v(i), 2);

            d->Luu(2 * foot, 2 * foot + 1) += vel_weight_ * oRh_(1, 0) * oRh_(1, 1) * pow(alpha_v(i), 2);
            d->Luu(2 * foot + 1, 2 * foot) += vel_weight_ * oRh_(1, 0) * oRh_(1, 1) * pow(alpha_v(i), 2);

            d->Lx(12 + 2 * foot) += vel_weight_ * oRh_(1, 0) * alpha_v(i) * rb_vely_max_(i, 2 * foot);
            d->Lx(12 + 2 * foot + 1) += vel_weight_ * oRh_(1, 1) * alpha_v(i) * rb_vely_max_(i, 2 * foot);

            d->Lxx(12 + 2 * foot, 12 + 2 * foot) += vel_weight_ * pow(oRh_(1, 0), 2) * pow(alpha_v(i), 2);
            d->Lxx(12 + 2 * foot + 1, 12 + 2 * foot + 1) += vel_weight_ * pow(oRh_(1, 1), 2) * pow(alpha_v(i), 2);

            d->Lxx(12 + 2 * foot, 12 + 2 * foot + 1) += vel_weight_ * oRh_(1, 0) * oRh_(1, 1) * pow(alpha_v(i), 2);
            d->Lxx(12 + 2 * foot + 1, 12 + 2 * foot) += vel_weight_ * oRh_(1, 0) * oRh_(1, 1) * pow(alpha_v(i), 2);

            d->Lxu(12 + 2 * foot, 2 * foot) += vel_weight_ * pow(alpha_v(i), 2) * pow(oRh_(1, 0), 2);
            d->Lxu(12 + 2 * foot, 2 * foot + 1) += vel_weight_ * pow(alpha_v(i), 2) * oRh_(1, 0) * oRh_(1, 1);
            d->Lxu(12 + 2 * foot + 1, 2 * foot) += vel_weight_ * pow(alpha_v(i), 2) * oRh_(1, 0) * oRh_(1, 1);
            d->Lxu(12 + 2 * foot + 1, 2 * foot + 1) += vel_weight_ * pow(alpha_v(i), 2) * pow(oRh_(1, 1), 2);
          }
          if (rb_vely_max_bool_(i, 2 * foot + 1)) {
            d->Lu(2 * foot) += -vel_weight_ * oRh_(1, 0) * alpha_v(i) * rb_vely_max_(i, 2 * foot + 1);
            d->Lu(2 * foot + 1) += -vel_weight_ * oRh_(1, 1) * alpha_v(i) * rb_vely_max_(i, 2 * foot + 1);

            d->Luu(2 * foot, 2 * foot) += vel_weight_ * pow(oRh_(1, 0), 2) * pow(alpha_v(i), 2);
            d->Luu(2 * foot + 1, 2 * foot + 1) += vel_weight_ * pow(oRh_(1, 1), 2) * pow(alpha_v(i), 2);

            d->Luu(2 * foot, 2 * foot + 1) += vel_weight_ * oRh_(1, 0) * oRh_(1, 1) * pow(alpha_v(i), 2);
            d->Luu(2 * foot + 1, 2 * foot) += vel_weight_ * oRh_(1, 0) * oRh_(1, 1) * pow(alpha_v(i), 2);

            d->Lx(12 + 2 * foot) += -vel_weight_ * oRh_(1, 0) * alpha_v(i) * rb_vely_max_(i, 2 * foot + 1);
            d->Lx(12 + 2 * foot + 1) += -vel_weight_ * oRh_(1, 1) * alpha_v(i) * rb_vely_max_(i, 2 * foot + 1);

            d->Lxx(12 + 2 * foot, 12 + 2 * foot) += vel_weight_ * pow(oRh_(1, 0), 2) * pow(alpha_v(i), 2);
            d->Lxx(12 + 2 * foot + 1, 12 + 2 * foot + 1) += vel_weight_ * pow(oRh_(1, 1), 2) * pow(alpha_v(i), 2);

            d->Lxx(12 + 2 * foot, 12 + 2 * foot + 1) += vel_weight_ * oRh_(1, 0) * oRh_(1, 1) * pow(alpha_v(i), 2);
            d->Lxx(12 + 2 * foot + 1, 12 + 2 * foot) += vel_weight_ * oRh_(1, 0) * oRh_(1, 1) * pow(alpha_v(i), 2);

            d->Lxu(12 + 2 * foot, 2 * foot) += vel_weight_ * pow(alpha_v(i), 2) * pow(oRh_(1, 0), 2);
            d->Lxu(12 + 2 * foot, 2 * foot + 1) += vel_weight_ * pow(alpha_v(i), 2) * oRh_(1, 0) * oRh_(1, 1);
            d->Lxu(12 + 2 * foot + 1, 2 * foot) += vel_weight_ * pow(alpha_v(i), 2) * oRh_(1, 0) * oRh_(1, 1);
            d->Lxu(12 + 2 * foot + 1, 2 * foot + 1) += vel_weight_ * pow(alpha_v(i), 2) * pow(oRh_(1, 1), 2);
          }
        }
      }
    }
  }

  if (is_jerk_activated_) {
    for (int foot = 0; foot < 4; foot++) {
      if (S_[foot] == Scalar(1)) {
        // derivatives relative to jerk on x axis
        d->Lu(2 * foot) += jerk_weight_ * oRh_(0, 0) * alpha_j * rb_jerk_(0, foot);
        d->Lu(2 * foot + 1) += jerk_weight_ * oRh_(0, 1) * alpha_j * rb_jerk_(0, foot);

        d->Luu(2 * foot, 2 * foot) += jerk_weight_ * pow(oRh_(0, 0), 2) * pow(alpha_j, 2);
        d->Luu(2 * foot + 1, 2 * foot + 1) += jerk_weight_ * pow(oRh_(0, 1), 2) * pow(alpha_j, 2);

        d->Luu(2 * foot, 2 * foot + 1) += jerk_weight_ * oRh_(0, 0) * oRh_(0, 1) * pow(alpha_j, 2);
        d->Luu(2 * foot + 1, 2 * foot) += jerk_weight_ * oRh_(0, 0) * oRh_(0, 1) * pow(alpha_j, 2);

        d->Lx(12 + 2 * foot) += jerk_weight_ * oRh_(0, 0) * alpha_j * rb_jerk_(0, foot);
        d->Lx(12 + 2 * foot + 1) += jerk_weight_ * oRh_(0, 1) * alpha_j * rb_jerk_(0, foot);

        d->Lxx(12 + 2 * foot, 12 + 2 * foot) += jerk_weight_ * pow(oRh_(0, 0), 2) * pow(alpha_j, 2);
        d->Lxx(12 + 2 * foot + 1, 12 + 2 * foot + 1) += jerk_weight_ * pow(oRh_(0, 1), 2) * pow(alpha_j, 2);

        d->Lxx(12 + 2 * foot, 12 + 2 * foot + 1) += jerk_weight_ * oRh_(0, 0) * oRh_(0, 1) * pow(alpha_j, 2);
        d->Lxx(12 + 2 * foot + 1, 12 + 2 * foot) += jerk_weight_ * oRh_(0, 0) * oRh_(0, 1) * pow(alpha_j, 2);

        d->Lxu(12 + 2 * foot, 2 * foot) += jerk_weight_ * pow(alpha_j, 2) * pow(oRh_(0, 0), 2);
        d->Lxu(12 + 2 * foot, 2 * foot + 1) += jerk_weight_ * pow(alpha_j, 2) * oRh_(0, 0) * oRh_(0, 1);
        d->Lxu(12 + 2 * foot + 1, 2 * foot) += jerk_weight_ * pow(alpha_j, 2) * oRh_(0, 0) * oRh_(0, 1);
        d->Lxu(12 + 2 * foot + 1, 2 * foot + 1) += jerk_weight_ * pow(alpha_j, 2) * pow(oRh_(0, 1), 2);

        // derivatives relative to jerk on y axis
        d->Lu(2 * foot) += jerk_weight_ * oRh_(1, 0) * alpha_j * rb_jerk_(1, foot);
        d->Lu(2 * foot + 1) += jerk_weight_ * oRh_(1, 1) * alpha_j * rb_jerk_(1, foot);

        d->Luu(2 * foot, 2 * foot) += jerk_weight_ * pow(oRh_(1, 0), 2) * pow(alpha_j, 2);
        d->Luu(2 * foot + 1, 2 * foot + 1) += jerk_weight_ * pow(oRh_(1, 1), 2) * pow(alpha_j, 2);

        d->Luu(2 * foot, 2 * foot + 1) += jerk_weight_ * oRh_(1, 0) * oRh_(1, 1) * pow(alpha_j, 2);
        d->Luu(2 * foot + 1, 2 * foot) += jerk_weight_ * oRh_(1, 0) * oRh_(1, 1) * pow(alpha_j, 2);

        d->Lx(12 + 2 * foot) += jerk_weight_ * oRh_(1, 0) * alpha_j * rb_jerk_(1, foot);
        d->Lx(12 + 2 * foot + 1) += jerk_weight_ * oRh_(1, 1) * alpha_j * rb_jerk_(1, foot);

        d->Lxx(12 + 2 * foot, 12 + 2 * foot) += jerk_weight_ * pow(oRh_(1, 0), 2) * pow(alpha_j, 2);
        d->Lxx(12 + 2 * foot + 1, 12 + 2 * foot + 1) += jerk_weight_ * pow(oRh_(1, 1), 2) * pow(alpha_j, 2);

        d->Lxx(12 + 2 * foot, 12 + 2 * foot + 1) += jerk_weight_ * oRh_(1, 0) * oRh_(1, 1) * pow(alpha_j, 2);
        d->Lxx(12 + 2 * foot + 1, 12 + 2 * foot) += jerk_weight_ * oRh_(1, 0) * oRh_(1, 1) * pow(alpha_j, 2);

        d->Lxu(12 + 2 * foot, 2 * foot) += jerk_weight_ * pow(alpha_j, 2) * pow(oRh_(1, 0), 2);
        d->Lxu(12 + 2 * foot, 2 * foot + 1) += jerk_weight_ * pow(alpha_j, 2) * oRh_(1, 0) * oRh_(1, 1);
        d->Lxu(12 + 2 * foot + 1, 2 * foot) += jerk_weight_ * pow(alpha_j, 2) * oRh_(1, 0) * oRh_(1, 1);
        d->Lxu(12 + 2 * foot + 1, 2 * foot + 1) += jerk_weight_ * pow(alpha_j, 2) * pow(oRh_(1, 1), 2);
      }
    }
  }

  // Dynamic derivatives
  d->Fx.setIdentity();
  d->Fu.block(12, 0, 8, 8) = B;
}

template <typename Scalar>
boost::shared_ptr<crocoddyl::ActionDataAbstractTpl<Scalar> > ActionModelQuadrupedStepTpl<Scalar>::createData() {
  return boost::make_shared<ActionDataQuadrupedStepTpl<Scalar> >(this);
}

////////////////////////////////
// get & set parameters ////////
////////////////////////////////

template <typename Scalar>
const typename Eigen::Matrix<Scalar, 12, 1>& ActionModelQuadrupedStepTpl<Scalar>::get_state_weights() const {
  return state_weights_;
}
template <typename Scalar>
void ActionModelQuadrupedStepTpl<Scalar>::set_state_weights(const typename MathBase::VectorXs& weights) {
  if (static_cast<std::size_t>(weights.size()) != 12) {
    throw_pretty("Invalid argument: "
                 << "Weights vector has wrong dimension (it should be 12)");
  }
  state_weights_ = weights;
}

template <typename Scalar>
const typename Eigen::Matrix<Scalar, 8, 1>& ActionModelQuadrupedStepTpl<Scalar>::get_step_weights() const {
  return step_weights_;
}
template <typename Scalar>
void ActionModelQuadrupedStepTpl<Scalar>::set_step_weights(const typename MathBase::VectorXs& weights) {
  if (static_cast<std::size_t>(weights.size()) != 8) {
    throw_pretty("Invalid argument: "
                 << "Weights vector has wrong dimension (it should be 4)");
  }
  step_weights_ = weights;
}

template <typename Scalar>
const typename Eigen::Matrix<Scalar, 8, 1>& ActionModelQuadrupedStepTpl<Scalar>::get_heuristic_weights() const {
  return heuristic_weights_;
}
template <typename Scalar>
void ActionModelQuadrupedStepTpl<Scalar>::set_heuristic_weights(const typename MathBase::VectorXs& weights) {
  if (static_cast<std::size_t>(weights.size()) != 8) {
    throw_pretty("Invalid argument: "
                 << "Weights vector has wrong dimension (it should be 8)");
  }
  heuristic_weights_ = weights;
}

template <typename Scalar>
const bool& ActionModelQuadrupedStepTpl<Scalar>::get_symmetry_term() const {
  return symmetry_term;
}
template <typename Scalar>
void ActionModelQuadrupedStepTpl<Scalar>::set_symmetry_term(const bool& sym_term) {
  // The model need to be updated after this changed
  symmetry_term = sym_term;
}

template <typename Scalar>
const bool& ActionModelQuadrupedStepTpl<Scalar>::get_centrifugal_term() const {
  return centrifugal_term;
}
template <typename Scalar>
void ActionModelQuadrupedStepTpl<Scalar>::set_centrifugal_term(const bool& cent_term) {
  // The model need to be updated after this changed
  centrifugal_term = cent_term;
}

template <typename Scalar>
const Scalar& ActionModelQuadrupedStepTpl<Scalar>::get_T_gait() const {
  // The model need to be updated after this changed
  return T_gait;
}
template <typename Scalar>
void ActionModelQuadrupedStepTpl<Scalar>::set_T_gait(const Scalar& T_gait_) {
  // The model need to be updated after this changed
  T_gait = T_gait_;
}

template <typename Scalar>
const bool& ActionModelQuadrupedStepTpl<Scalar>::get_acc_activated() const {
  return is_acc_activated_;
}
template <typename Scalar>
void ActionModelQuadrupedStepTpl<Scalar>::set_acc_activated(const bool& is_activated) {
  is_acc_activated_ = is_activated;
}

template <typename Scalar>
const typename Eigen::Matrix<Scalar, 2, 1>& ActionModelQuadrupedStepTpl<Scalar>::get_acc_lim() const {
  return acc_lim_;
}
template <typename Scalar>
void ActionModelQuadrupedStepTpl<Scalar>::set_acc_lim(const typename MathBase::VectorXs& acceleration_lim_) {
  if (static_cast<std::size_t>(acceleration_lim_.size()) != 2) {
    throw_pretty("Invalid argument: "
                 << "Acceleration limit vector [ax_max, ay_max] has wrong dimension (it should be 2)");
  }
  acc_lim_ = acceleration_lim_;
}

template <typename Scalar>
const Scalar& ActionModelQuadrupedStepTpl<Scalar>::get_acc_weight() const {
  return acc_weight_;
}
template <typename Scalar>
void ActionModelQuadrupedStepTpl<Scalar>::set_acc_weight(const Scalar& weight_) {
  acc_weight_ = weight_;
}

template <typename Scalar>
const bool& ActionModelQuadrupedStepTpl<Scalar>::get_vel_activated() const {
  return is_vel_activated_;
}
template <typename Scalar>
void ActionModelQuadrupedStepTpl<Scalar>::set_vel_activated(const bool& is_activated) {
  is_vel_activated_ = is_activated;
}

template <typename Scalar>
const typename Eigen::Matrix<Scalar, 2, 1>& ActionModelQuadrupedStepTpl<Scalar>::get_vel_lim() const {
  return vel_lim_;
}
template <typename Scalar>
void ActionModelQuadrupedStepTpl<Scalar>::set_vel_lim(const typename MathBase::VectorXs& velocity_lim_) {
  if (static_cast<std::size_t>(velocity_lim_.size()) != 2) {
    throw_pretty("Invalid argument: "
                 << "Velocity limit vector [vx_max, vy_max] has wrong dimension (it should be 2)");
  }
  vel_lim_ = velocity_lim_;
}

template <typename Scalar>
const Scalar& ActionModelQuadrupedStepTpl<Scalar>::get_vel_weight() const {
  return vel_weight_;
}
template <typename Scalar>
void ActionModelQuadrupedStepTpl<Scalar>::set_vel_weight(const Scalar& weight_) {
  vel_weight_ = weight_;
}

template <typename Scalar>
const Scalar& ActionModelQuadrupedStepTpl<Scalar>::get_jerk_weight() const {
  return jerk_weight_;
}
template <typename Scalar>
void ActionModelQuadrupedStepTpl<Scalar>::set_jerk_weight(const Scalar& weight_) {
  jerk_weight_ = weight_;
}

template <typename Scalar>
const bool& ActionModelQuadrupedStepTpl<Scalar>::get_jerk_activated() const {
  return is_jerk_activated_;
}
template <typename Scalar>
void ActionModelQuadrupedStepTpl<Scalar>::set_jerk_activated(const bool& is_activated) {
  is_jerk_activated_ = is_activated;
}

template <typename Scalar>
void ActionModelQuadrupedStepTpl<Scalar>::set_sample_feet_traj(const int& n_sample) {
  N_sampling = n_sample;

  // Acceleration cost
  delta_ = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(N_sampling - 1, 4);
  gamma_ = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(N_sampling - 1, 3);
  for (int k = 1; k < N_sampling; k++) {
    delta_(k - 1, 0) = (float)k / (float)N_sampling;  // [1/N, 2/N, ... , (N-1)/N]
  }
  delta_.col(1) << delta_.col(0).pow(2);
  delta_.col(2) << delta_.col(0).pow(3);
  delta_.col(3) << delta_.col(0).pow(4);

  gamma_.col(0) = 60 * delta_.col(0) - 180 * delta_.col(1) + 120 * delta_.col(2);
  gamma_.col(1) = -36 * delta_.col(0) + 96 * delta_.col(1) - 60 * delta_.col(2);
  gamma_.col(2) = -9 * delta_.col(0) + 18 * delta_.col(1) - 10 * delta_.col(2);

  alpha_ = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(N_sampling - 1, 1);   // Common for 4 feet
  beta_x_ = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(N_sampling - 1, 4);  // Depends on ao,vo of feet
  beta_y_ = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(N_sampling - 1, 4);  // Depends on ao,vo of feet
  tmp_ones_ = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Ones(N_sampling - 1, 1);

  rb_accx_max_ = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(N_sampling - 1, 8);
  rb_accy_max_ = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(N_sampling - 1, 8);
  rb_accx_max_bool_ = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(N_sampling - 1, 8);
  rb_accy_max_bool_ = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(N_sampling - 1, 8);

  // Velocity cost
  gamma_v = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(N_sampling - 1, 4);
  gamma_v.col(0) = 30 * delta_.col(1) - 60 * delta_.col(2) + 30 * delta_.col(3);
  gamma_v.col(1) = delta_.col(0);
  gamma_v.col(2) = -18 * delta_.col(1) + 32 * delta_.col(2) - 15 * delta_.col(3);
  gamma_v.col(3) = -4.5 * delta_.col(1) + 6 * delta_.col(2) - 2.5 * delta_.col(3);

  alpha_v = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(N_sampling - 1, 1);  // Common for 4 feet
  beta_x_v =
      Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(N_sampling - 1, 4);  // Depends on a0_x, v0_x of feet
  beta_y_v =
      Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(N_sampling - 1, 4);  // Depends on a0_y, v0_y of feet

  rb_velx_max_ = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(N_sampling - 1, 8);
  rb_vely_max_ = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(N_sampling - 1, 8);
  rb_velx_max_bool_ = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(N_sampling - 1, 8);
  rb_vely_max_bool_ = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(N_sampling - 1, 8);
}

////////////////////////
// Update current model
////////////////////////

template <typename Scalar>
void ActionModelQuadrupedStepTpl<Scalar>::update_model(
    const Eigen::Ref<const typename MathBase::MatrixXs>& l_feet,
    const Eigen::Ref<const typename MathBase::MatrixXs>& xref, const Eigen::Ref<const typename MathBase::VectorXs>& S,
    const Eigen::Ref<const typename MathBase::MatrixXs>& position,
    const Eigen::Ref<const typename MathBase::MatrixXs>& velocity,
    const Eigen::Ref<const typename MathBase::MatrixXs>& acceleration,
    const Eigen::Ref<const typename MathBase::MatrixXs>& jerk,
    const Eigen::Ref<const typename MathBase::MatrixXs>& oRh, const Eigen::Ref<const typename MathBase::MatrixXs>& oTh,
    const Scalar& delta_T) {
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

  // Velocity :  [[vx_0, vx_1, vx_2, vx_3],
  //              [vy_0, vy_1, vy_2, vy_3],
  //              [vz_0, vz_1, vz_2, vz_3]]

  xref_ = xref;
  S_ = S;
  position_ = position;
  oRh_ = oRh;
  oTh_ = oTh;
  jerk_ = jerk;

  /* R_tmp << cos(xref(5, 0)), -sin(xref(5, 0)), Scalar(0), sin(xref(5, 0)), cos(xref(5, 0)), Scalar(0), Scalar(0),
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
  } */

  for (int i = 0; i < 4; i = i + 1) {
    pheuristic_[2 * i] = l_feet(0, i);
    pheuristic_[2 * i + 1] = l_feet(1, i);
  }

  /* std::cout << pshoulder_ << std::endl; */

  B.setZero();

  if (S[0] == Scalar(1)) {
    B.block(0, 0, 2, 2).setIdentity();
  }
  if (S[1] == Scalar(1)) {
    B.block(2, 2, 2, 2).setIdentity();
  }
  if (S[2] == Scalar(1)) {
    B.block(4, 4, 2, 2).setIdentity();
  }
  if (S[3] == Scalar(1)) {
    B.block(6, 6, 2, 2).setIdentity();
  }

  alpha_ = (1 / pow(delta_T, 2)) * gamma_.col(0);
  alpha_j = (60 / pow(delta_T, 3));
  alpha_v = (1 / delta_T) * gamma_v.col(0);

  for (int i = 0; i < 4; i++) {
    if (S[i] == Scalar(1) && is_acc_activated_) {
      beta_x_.col(i) = acceleration(0, i) * tmp_ones_ + (velocity(0, i) / delta_T) * gamma_.col(1) +
                       acceleration(0, i) * gamma_.col(2);
      beta_y_.col(i) = acceleration(1, i) * tmp_ones_ + (velocity(1, i) / delta_T) * gamma_.col(1) +
                       acceleration(1, i) * gamma_.col(2);
    } else {
      beta_x_.col(i).setZero();
      beta_y_.col(i).setZero();
    }

    if (S[i] == Scalar(1) && is_vel_activated_) {
      beta_x_v.col(i) = velocity(0, i) * tmp_ones_ + (acceleration(0, i) * delta_T) * gamma_v.col(1) +
                        velocity(0, i) * gamma_.col(2) + (acceleration(0, i) * delta_T) * gamma_v.col(3);
      beta_y_v.col(i) = velocity(1, i) * tmp_ones_ + (acceleration(1, i) * delta_T) * gamma_v.col(1) +
                        velocity(1, i) * gamma_.col(2) + (acceleration(1, i) * delta_T) * gamma_v.col(3);
    } else {
      beta_x_v.col(i).setZero();
      beta_y_v.col(i).setZero();
    }

    if (S[i] == Scalar(1) && is_jerk_activated_) {
      beta_j(0, i) = -(36 * velocity(0, i)) / pow(delta_T, 2) - (9 * acceleration(0, i)) / delta_T;
      beta_j(1, i) = -(36 * velocity(1, i)) / pow(delta_T, 2) - (9 * acceleration(1, i)) / delta_T;
    } else {
      beta_j.col(i).setZero();
    }
  }
}
}  // namespace quadruped_walkgen

#endif
