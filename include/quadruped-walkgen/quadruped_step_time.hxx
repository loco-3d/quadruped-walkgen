#ifndef __quadruped_walkgen_quadruped_step_time_hxx__
#define __quadruped_walkgen_quadruped_step_time_hxx__

#include "crocoddyl/core/utils/exception.hpp"

namespace quadruped_walkgen {
template <typename Scalar>
ActionModelQuadrupedStepTimeTpl<Scalar>::ActionModelQuadrupedStepTimeTpl()
    : crocoddyl::ActionModelAbstractTpl<Scalar>(boost::make_shared<crocoddyl::StateVectorTpl<Scalar> >(21), 8, 29) {
  B.setZero();  // x_next = x + B * u
  rub_max_.setZero();
  rub_max_bool.setZero();

  state_weights_ << Scalar(1.), Scalar(1.), Scalar(150.), Scalar(35.), Scalar(30.), Scalar(8.), Scalar(20.),
      Scalar(20.), Scalar(15.), Scalar(4.), Scalar(4.), Scalar(8.);
  heuristicWeights.setConstant(Scalar(0.));
  step_weights_.setConstant(Scalar(1));
  pheuristic_.setZero();

  // Compute heuristic inside update Model
  // pshoulder_0 <<  Scalar(0.1946) ,   Scalar(0.1946) ,   Scalar(-0.1946),  Scalar(-0.1946) ,
  //                 Scalar(0.15005) ,  Scalar(-0.15005)  , Scalar(0.15005)  ,  Scalar(-0.15005) ;
  // pshoulder_tmp.setZero() ;
  // pcentrifugal_tmp_1.setZero() ;
  // pcentrifugal_tmp_2.setZero() ;
  // pcentrifugal_tmp.setZero() ;
  // T_gait = Scalar(0.64) ;
  centrifugal_term = true;
  symmetry_term = true;

  // Weight on the speed ot the feet
  nb_nodes = Scalar(15.);
  vlim = Scalar(2.);
  beta_lim = Scalar((64 * nb_nodes * nb_nodes * vlim * vlim) / 225);  // apparent speed used in the cost function
  speed_weight = Scalar(10.);

  // Logging cost
  cost_.setZero();
  log_cost = true;

  // indicates whether it t the 1st step, otherwise the cost function is much simpler (acc, speed = 0)
  first_step = false;

  // Coefficients for sample velocity of the feet
  nb_alpha_ = 4;
  alpha = MathBase::ArrayXs::Zero(nb_alpha_);
  alpha2 = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(nb_alpha_, 4);
  b_coeff = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(
      nb_alpha_, 3);  // Constant for all feet, avoid re-computing them

  // Cost = DT * b0(alpha) + DT**2 * b1(alpha) + DX * b2(alpha) for x velocity
  b_coeff_x0 = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(nb_alpha_, 4);  // col(i) --> foot i
  b_coeff_x1 = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(nb_alpha_, 4);
  b_coeff_x2 = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(nb_alpha_, 4);

  // Cost = DT * b0(alpha) + DT**2 * b1(alpha) + DX * b2(alpha) for y velocity
  b_coeff_y0 = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(nb_alpha_, 4);  // col(i) --> foot i
  b_coeff_y1 = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(nb_alpha_, 4);
  b_coeff_y2 = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(nb_alpha_, 4);

  rub_max_first_x = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(nb_alpha_, 4);
  rub_max_first_y = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(nb_alpha_, 4);
  rub_max_first_2 = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(nb_alpha_, 4);
  rub_max_first_bool = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(nb_alpha_, 4);

  alpha.setLinSpaced(nb_alpha_, Scalar(0.0), Scalar(1.0));
  alpha2.col(0) << alpha;
  alpha2.col(1) << alpha.pow(2);
  alpha2.col(2) << alpha.pow(3);
  alpha2.col(3) << alpha.pow(4);

  b_coeff.col(0) =
      Scalar(1.0) - Scalar(18.) * alpha2.col(1) + Scalar(32.) * alpha2.col(2) - Scalar(15.) * alpha2.col(3);
  b_coeff.col(1) =
      alpha2.col(0) - Scalar(4.5) * alpha2.col(1) + Scalar(6.) * alpha2.col(2) - Scalar(2.5) * alpha2.col(3);
  b_coeff.col(2) = Scalar(30.) * alpha2.col(1) - Scalar(60.) * alpha2.col(2) + Scalar(30.) * alpha2.col(3);

  lfeet.setZero();
}

template <typename Scalar>
ActionModelQuadrupedStepTimeTpl<Scalar>::~ActionModelQuadrupedStepTimeTpl() {}

template <typename Scalar>
void ActionModelQuadrupedStepTimeTpl<Scalar>::calc(
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

  ActionDataQuadrupedStepTimeTpl<Scalar>* d = static_cast<ActionDataQuadrupedStepTimeTpl<Scalar>*>(data.get());

  // Update position of the feet
  d->xnext.template head<12>() = x.head(12);
  d->xnext.template segment<8>(12) = x.segment(12, 8) + B * u;
  d->xnext.template tail<1>() = x.tail(1);

  // Residual cost on the state and force norm
  d->r.template head<12>() = state_weights_.cwiseProduct(x.head(12) - xref_);
  d->r.template segment<8>(12) = heuristicWeights.cwiseProduct(
      x.segment(12, 8) -
      pheuristic_);  // Not used, set to 0, S matrix is for moving feet and not feet already on the ground
  d->r.template tail<4>() = step_weights_.cwiseProduct(u);

  d->cost = Scalar(0.5) * d->r.transpose() * d->r;

  if (first_step) {
    for (int i = 0; i < 4; i++) {
      if (S_[i] == Scalar(1)) {
        rub_max_first_x.col(i) =
            x(20) * b_coeff_x0.col(i) + x(20) * x(20) * b_coeff_x1.col(i) + u(2 * i) * b_coeff_x2.col(i);
        rub_max_first_y.col(i) =
            x(20) * b_coeff_y0.col(i) + x(20) * x(20) * b_coeff_y1.col(i) + u(2 * i + 1) * b_coeff_y2.col(i);

        rub_max_first_2.col(i) = rub_max_first_x.col(i).pow(2) + rub_max_first_y.col(i).pow(2) -
                                 x(20) * x(20) * vlim * vlim * nb_nodes * nb_nodes;
      } else {
        rub_max_first_2.col(i).setZero();
      }
    }

    rub_max_first_bool = (rub_max_first_2 > Scalar(0.)).template cast<Scalar>();  // Usefull to compute the derivatives
    rub_max_first_2 = rub_max_first_2.cwiseMax(Scalar(0.));                       // Remove <0 terms

    for (int i = 0; i < nb_alpha_; i++) {
      d->cost += speed_weight * Scalar(0.5) * rub_max_first_2.row(i).sum();
    }
  } else {
    rub_max_ << u[0] * u[0] + u[1] * u[1] - beta_lim * x[20] * x[20],
        u[2] * u[2] + u[3] * u[3] - beta_lim * x[20] * x[20];

    rub_max_bool = (rub_max_.array() >= Scalar(0.)).matrix().template cast<Scalar>();
    rub_max_ = rub_max_.cwiseMax(Scalar(0.));

    d->cost += speed_weight * Scalar(0.5) * rub_max_.sum();
  }

  if (log_cost) {
    cost_[3] = 0;
    // Length to be consistent with others models
    cost_[0] = Scalar(0.5) * d->r.head(12).transpose() * d->r.head(12);              // State cost
    cost_[1] = Scalar(0.5) * d->r.segment(12, 8).transpose() * d->r.segment(12, 8);  // heuristic cost
    cost_[2] = Scalar(0.5) * d->r.tail(4).transpose() * d->r.tail(4);                // Delta feet cost

    if (first_step) {
      for (int i = 0; i < 3; i++) {
        cost_[3] += speed_weight * Scalar(0.5) * rub_max_first_2.row(i).sum();
      }
    } else {
      cost_[3] = speed_weight * Scalar(0.5) * rub_max_.sum();
    }
  }
}

template <typename Scalar>
void ActionModelQuadrupedStepTimeTpl<Scalar>::calcDiff(
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

  ActionDataQuadrupedStepTimeTpl<Scalar>* d = static_cast<ActionDataQuadrupedStepTimeTpl<Scalar>*>(data.get());

  d->Lx.setZero();
  d->Lu.setZero();
  d->Lxu.setZero();
  d->Lxx.setZero();
  d->Luu.setZero();
  // Cost derivatives : Lx
  d->Lx.template head<12>() = (state_weights_.array() * d->r.template head<12>().array()).matrix();
  d->Lx.template segment<8>(12) = (heuristicWeights.array() * d->r.template segment<8>(12).array()).matrix();

  if (first_step) {
    for (int foot = 0; foot < 4; foot++) {
      if (S_[foot] == Scalar(1)) {
        for (int i = 0; i < nb_alpha_; i++) {
          if (rub_max_first_bool(i, foot)) {
            d->Lx(20) += speed_weight * (b_coeff_x0(i, foot) + Scalar(2) * x(20) * b_coeff_x1(i, foot)) *
                             rub_max_first_x(i, foot) +
                         speed_weight * (b_coeff_y0(i, foot) + Scalar(2) * x(20) * b_coeff_y1(i, foot)) *
                             rub_max_first_y(i, foot) -
                         speed_weight * x(20) * vlim * vlim * nb_nodes * nb_nodes;
            d->Lu(2 * foot) += speed_weight * b_coeff_x2(i, foot) * rub_max_first_x(i, foot);
            d->Lu(2 * foot + 1) += speed_weight * b_coeff_y2(i, foot) * rub_max_first_y(i, foot);

            d->Luu(2 * foot, 2 * foot) += speed_weight * b_coeff_x2(i, foot) * b_coeff_x2(i, foot);
            d->Luu(2 * foot + 1, 2 * foot + 1) += speed_weight * b_coeff_y2(i, foot) * b_coeff_y2(i, foot);
            d->Lxu(20, 2 * foot) +=
                speed_weight * (b_coeff_x0(i, foot) + Scalar(2) * x(20) * b_coeff_x1(i, foot)) * b_coeff_x2(i, foot);
            d->Lxu(20, 2 * foot + 1) +=
                speed_weight * (b_coeff_y0(i, foot) + Scalar(2) * x(20) * b_coeff_y1(i, foot)) * b_coeff_y2(i, foot);
            d->Lxx(20, 20) +=
                speed_weight * std::pow(b_coeff_x0(i, foot) + Scalar(2) * x(20) * b_coeff_x1(i, foot), 2) +
                speed_weight * Scalar(2) * b_coeff_x1(i, foot) * rub_max_first_x(i, foot) +
                speed_weight * std::pow(b_coeff_y0(i, foot) + Scalar(2) * x(20) * b_coeff_x1(i, foot), 2) +
                speed_weight * Scalar(2) * b_coeff_y1(i, foot) * rub_max_first_y(i, foot) -
                speed_weight * vlim * vlim * nb_nodes * nb_nodes;
          }
        }
      }
    }

  }

  else {
    d->Lx.template tail<1>() << -beta_lim * speed_weight * x(20) * rub_max_bool[0] -
                                    beta_lim * speed_weight * x(20) * rub_max_bool[1];

    d->Lu << speed_weight * u[0] * rub_max_bool[0], speed_weight * u[1] * rub_max_bool[0],
        speed_weight * u[2] * rub_max_bool[1], speed_weight * u[3] * rub_max_bool[1];

    d->Lxx(20, 20) = -beta_lim * speed_weight * rub_max_bool[0] - beta_lim * speed_weight * rub_max_bool[1];

    d->Luu.diagonal() << speed_weight * rub_max_bool[0], speed_weight * rub_max_bool[0],
        speed_weight * rub_max_bool[1], speed_weight * rub_max_bool[1];
  }

  d->Lu += (step_weights_.array() * d->r.template tail<4>().array()).matrix();

  // Hessian : Lxx
  d->Lxx.diagonal().head(12) = (state_weights_.array() * state_weights_.array()).matrix();
  d->Lxx.diagonal().segment(12, 8) = (heuristicWeights.array() * heuristicWeights.array()).matrix();

  d->Luu.diagonal() += (step_weights_.array() * step_weights_.array()).matrix();

  // Dynamic derivatives
  d->Fx.setIdentity();
  d->Fu.block(12, 0, 8, 8) = B;
}

template <typename Scalar>
boost::shared_ptr<crocoddyl::ActionDataAbstractTpl<Scalar> > ActionModelQuadrupedStepTimeTpl<Scalar>::createData() {
  return boost::make_shared<ActionDataQuadrupedStepTimeTpl<Scalar> >(this);
}

////////////////////////////////
// get & set parameters ////////
////////////////////////////////

template <typename Scalar>
const typename Eigen::Matrix<Scalar, 12, 1>& ActionModelQuadrupedStepTimeTpl<Scalar>::get_state_weights() const {
  return state_weights_;
}
template <typename Scalar>
void ActionModelQuadrupedStepTimeTpl<Scalar>::set_state_weights(const typename MathBase::VectorXs& weights) {
  if (static_cast<std::size_t>(weights.size()) != 12) {
    throw_pretty("Invalid argument: "
                 << "Weights vector has wrong dimension (it should be 12)");
  }
  state_weights_ = weights;
}

template <typename Scalar>
const typename Eigen::Matrix<Scalar, 4, 1>& ActionModelQuadrupedStepTimeTpl<Scalar>::get_step_weights() const {
  return step_weights_;
}
template <typename Scalar>
void ActionModelQuadrupedStepTimeTpl<Scalar>::set_step_weights(const typename MathBase::VectorXs& weights) {
  if (static_cast<std::size_t>(weights.size()) != 8) {
    throw_pretty("Invalid argument: "
                 << "Weights vector has wrong dimension (it should be 8)");
  }
  step_weights_ = weights;
}

template <typename Scalar>
const typename Eigen::Matrix<Scalar, 8, 1>& ActionModelQuadrupedStepTimeTpl<Scalar>::get_heuristic_weights() const {
  return heuristicWeights;
}
template <typename Scalar>
void ActionModelQuadrupedStepTimeTpl<Scalar>::set_heuristic_weights(const typename MathBase::VectorXs& weights) {
  if (static_cast<std::size_t>(weights.size()) != 8) {
    throw_pretty("Invalid argument: "
                 << "Weights vector has wrong dimension (it should be 8)");
  }
  heuristicWeights = weights;
}

template <typename Scalar>
const bool& ActionModelQuadrupedStepTimeTpl<Scalar>::get_symmetry_term() const {
  return symmetry_term;
}
template <typename Scalar>
void ActionModelQuadrupedStepTimeTpl<Scalar>::set_symmetry_term(const bool& sym_term) {
  // The model need to be updated after this changed
  symmetry_term = sym_term;
}

template <typename Scalar>
const bool& ActionModelQuadrupedStepTimeTpl<Scalar>::get_centrifugal_term() const {
  return centrifugal_term;
}
template <typename Scalar>
void ActionModelQuadrupedStepTimeTpl<Scalar>::set_centrifugal_term(const bool& cent_term) {
  // The model need to be updated after this changed
  centrifugal_term = cent_term;
}

template <typename Scalar>
const Scalar& ActionModelQuadrupedStepTimeTpl<Scalar>::get_T_gait() const {
  // The model need to be updated after this changed
  return T_gait;
}
template <typename Scalar>
void ActionModelQuadrupedStepTimeTpl<Scalar>::set_T_gait(const Scalar& T_gait_) {
  // The model need to be updated after this changed
  T_gait = T_gait_;
}

/////////////////////////////////////////////
// Get and modify param in speed cost      //
/////////////////////////////////////////////
template <typename Scalar>
const Scalar& ActionModelQuadrupedStepTimeTpl<Scalar>::get_speed_weight() const {
  return speed_weight;
}
template <typename Scalar>
void ActionModelQuadrupedStepTimeTpl<Scalar>::set_speed_weight(const Scalar& weight_) {
  speed_weight = weight_;
}

template <typename Scalar>
const Scalar& ActionModelQuadrupedStepTimeTpl<Scalar>::get_nb_nodes() const {
  return nb_nodes;
}
template <typename Scalar>
void ActionModelQuadrupedStepTimeTpl<Scalar>::set_nb_nodes(const Scalar& nodes_) {
  nb_nodes = nodes_;
  beta_lim = Scalar((64 * nb_nodes * nb_nodes * vlim * vlim) / 225);
  ;
}

template <typename Scalar>
const Scalar& ActionModelQuadrupedStepTimeTpl<Scalar>::get_vlim() const {
  return vlim;
}
template <typename Scalar>
void ActionModelQuadrupedStepTimeTpl<Scalar>::set_vlim(const Scalar& vlim_) {
  vlim = vlim_;
  beta_lim = Scalar((64 * nb_nodes * nb_nodes * vlim * vlim) / 225);
  ;
}

///////////////
// Log cost  //
///////////////
template <typename Scalar>
const typename Eigen::Matrix<Scalar, 7, 1>& ActionModelQuadrupedStepTimeTpl<Scalar>::get_cost() const {
  return cost_;
}

// indicates whether it t the 1st step, otherwise the cost function is much simpler (acc, speed = 0)
template <typename Scalar>
const bool& ActionModelQuadrupedStepTimeTpl<Scalar>::get_first_step() const {
  return first_step;
}
template <typename Scalar>
void ActionModelQuadrupedStepTimeTpl<Scalar>::set_first_step(const bool& first) {
  // The model need to be updated after this changed
  first_step = first;
}

////////////////////////
// Update current model
////////////////////////

template <typename Scalar>
void ActionModelQuadrupedStepTimeTpl<Scalar>::update_model(
    const Eigen::Ref<const typename MathBase::MatrixXs>& l_feet,
    const Eigen::Ref<const typename MathBase::MatrixXs>& velocity,
    const Eigen::Ref<const typename MathBase::MatrixXs>& acceleration,
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
  // Velocity :  [[vx_0, vx_1, vx_2, vx_3],
  //              [vy_0, vy_1, vy_2, vy_3],
  //              [vz_0, vz_1, vz_2, vz_3]]

  for (int i = 0; i < 4; i = i + 1) {
    pheuristic_.block(2 * i, 0, 2, 1) = l_feet.block(0, i, 2, 1);
  }
  xref_ = xref;
  S_ = S;

  for (int i = 0; i < 4; i++) {
    if (S[i] == Scalar(1)) {
      // Coeff for x velocity
      b_coeff_x0.col(i) = nb_nodes * velocity(0, i) * b_coeff.col(0);
      b_coeff_x1.col(i) = nb_nodes * nb_nodes * acceleration(0, i) * b_coeff.col(1);
      b_coeff_x2.col(i) = b_coeff.col(2);

      // Coeff for y velocity
      b_coeff_y0.col(i) = nb_nodes * velocity(1, i) * b_coeff.col(0);
      b_coeff_y1.col(i) = nb_nodes * nb_nodes * acceleration(1, i) * b_coeff.col(1);
      b_coeff_y2.col(i) = b_coeff.col(2);
    } else {
      b_coeff_x0.col(i).setZero();
      b_coeff_x1.col(i).setZero();
      b_coeff_x2.col(i).setZero();
      b_coeff_y0.col(i).setZero();
      b_coeff_y1.col(i).setZero();
      b_coeff_y2.col(i).setZero();
    }
  }

  // Compute heuristic inside update_model
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

  B.setZero();
  // Set B matrix according to the moving feet : S = gait - gait_old
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
}
}  // namespace quadruped_walkgen

#endif
