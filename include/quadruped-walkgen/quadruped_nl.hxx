#ifndef __quadruped_walkgen_quadruped_nl_hxx__
#define __quadruped_walkgen_quadruped_nl_hxx__

#include "crocoddyl/core/utils/exception.hpp"

namespace quadruped_walkgen {
template <typename Scalar>
ActionModelQuadrupedNonLinearTpl<Scalar>::ActionModelQuadrupedNonLinearTpl(
    typename Eigen::Matrix<Scalar, 3, 1> offset_CoM)
    : crocoddyl::ActionModelAbstractTpl<Scalar>(boost::make_shared<crocoddyl::StateVectorTpl<Scalar> >(12), 12, 24) {
  // Relative forces to compute the norm mof the command
  relative_forces = false;
  uref_.setZero();

  // Model parameters
  mu = Scalar(1);
  dt_ = Scalar(0.02);
  mass = Scalar(2.50000279);
  min_fz_in_contact = Scalar(0.0);
  max_fz = Scalar(25.);

  // Matrix model initialization
  g.setZero();
  g[8] = Scalar(-9.81) * dt_;
  gI.setZero();
  gI.diagonal() << Scalar(3.09249e-2), Scalar(5.106100e-2), Scalar(6.939757e-2);
  A.setIdentity();
  A.topRightCorner(6, 6) << Eigen::Matrix<Scalar, 6, 6>::Identity() * dt_;
  B.setZero();
  lever_arms.setZero();
  I_inv.setZero();

  // Weight vectors initialization
  force_weights_.setConstant(0.2);
  state_weights_ << Scalar(1.), Scalar(1.), Scalar(150.), Scalar(35.), Scalar(30.), Scalar(8.), Scalar(20.),
      Scalar(20.), Scalar(15.), Scalar(4.), Scalar(4.), Scalar(8.);
  friction_weight_ = Scalar(10);

  // UpperBound vector
  ub.setZero();
  for (int i = 0; i < 4; i = i + 1) {
    ub(6 * i + 5) = max_fz;
  }

  // Temporary vector used
  Fa_x_u.setZero();
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

  // Used for shoulder height weight
  pshoulder_0 << Scalar(0.1946), Scalar(0.1946), Scalar(-0.1946), Scalar(-0.1946), Scalar(0.14695), Scalar(-0.14695),
      Scalar(0.14695), Scalar(-0.14695);
  sh_hlim = Scalar(0.27);
  sh_weight = Scalar(10.);
  sh_ub_max_.setZero();
  psh.setZero();

  // Implicit integration
  // V+ = V + dt*B*u   ; P+ = P + dt*V+ != explicit : P+ = P + dt*V
  implicit_integration = true;
  offset_com = offset_CoM;  // x, y, z offset
}

template <typename Scalar>
ActionModelQuadrupedNonLinearTpl<Scalar>::~ActionModelQuadrupedNonLinearTpl() {}

template <typename Scalar>
void ActionModelQuadrupedNonLinearTpl<Scalar>::calc(
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

  ActionDataQuadrupedNonLinearTpl<Scalar>* d = static_cast<ActionDataQuadrupedNonLinearTpl<Scalar>*>(data.get());

  //  Update B :
  for (int i = 0; i < 4; i = i + 1) {
    if (gait(i, 0) != 0) {
      lever_tmp = lever_arms.block(0, i, 3, 1) - x.block(0, 0, 3, 1);
      R_tmp << 0.0, -lever_tmp[2], lever_tmp[1], lever_tmp[2], 0.0, -lever_tmp[0], -lever_tmp[1], lever_tmp[0], 0.0;
      B.block(9, 3 * i, 3, 3) << dt_ * I_inv * R_tmp;

      // Compute pdistance of the shoulder wrt contact point
      psh.block(0, i, 3, 1) << x[0] - offset_com(0, 0) + pshoulder_0(0, i) - pshoulder_0(1, i) * x[5] -
                                   lever_arms(0, i),
          x[1] - offset_com(1, 0) + pshoulder_0(1, i) + pshoulder_0(0, i) * x[5] - lever_arms(1, i),
          x[2] - offset_com(2, 0) + pshoulder_0(1, i) * x[3] - pshoulder_0(0, i) * x[4];
    } else {
      // Compute pdistance of the shoulder wrt contact point
      psh.block(0, i, 3, 1).setZero();
    }
  };

  // Discrete dynamic : A*x + B*u + g
  d->xnext << A.diagonal().cwiseProduct(x) + g;
  d->xnext.template head<6>() =
      d->xnext.template head<6>() + A.topRightCorner(6, 6).diagonal().cwiseProduct(x.tail(6));
  d->xnext.template tail<6>() = d->xnext.template tail<6>() + B.block(6, 0, 6, 12) * u;

  // Residual cost on the state and force norm
  d->r.template head<12>() = state_weights_.cwiseProduct(x - xref_);
  d->r.template tail<12>() = force_weights_.cwiseProduct(u - uref_);

  // Friction cone + shoulder height
  for (int i = 0; i < 4; i = i + 1) {
    Fa_x_u.segment(6 * i, 6) << u(3 * i) - mu * u(3 * i + 2), -u(3 * i) - mu * u(3 * i + 2),
        u(3 * i + 1) - mu * u(3 * i + 2), -u(3 * i + 1) - mu * u(3 * i + 2), -u(3 * i + 2), u(3 * i + 2);
  }
  rub_max_ = (Fa_x_u - ub).cwiseMax(Scalar(0.));

  // Shoulder height weight
  sh_ub_max_ << psh.block(0, 0, 3, 1).squaredNorm() - sh_hlim * sh_hlim,
      psh.block(0, 1, 3, 1).squaredNorm() - sh_hlim * sh_hlim, psh.block(0, 2, 3, 1).squaredNorm() - sh_hlim * sh_hlim,
      psh.block(0, 3, 3, 1).squaredNorm() - sh_hlim * sh_hlim;

  sh_ub_max_ = sh_ub_max_.cwiseMax(Scalar(0.));

  // Cost computation
  // d->cost = 0.5 * d->r.transpose() * d->r     + friction_weight_ * Scalar(0.5) * rub_max_.squaredNorm() + sh_weight
  // * Scalar(0.5) * sh_ub_max_.squaredNorm() ;
  d->cost = 0.5 * d->r.transpose() * d->r + friction_weight_ * Scalar(0.5) * rub_max_.squaredNorm() +
            sh_weight * Scalar(0.5) * sh_ub_max_.sum();
}

template <typename Scalar>
void ActionModelQuadrupedNonLinearTpl<Scalar>::calcDiff(
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

  ActionDataQuadrupedNonLinearTpl<Scalar>* d = static_cast<ActionDataQuadrupedNonLinearTpl<Scalar>*>(data.get());

  // Cost derivatives : Lx

  d->Lx = (state_weights_.array() * d->r.template head<12>().array()).matrix();

  // Hessian : Lxx
  d->Lxx.block(0, 0, 6, 6).setZero();
  d->Lxx.diagonal() = (state_weights_.array() * state_weights_.array()).matrix();

  // Shoulder height derivative cost
  for (int j = 0; j < 4; j = j + 1) {
    if (sh_ub_max_[j] > Scalar(0.)) {
      d->Lx(0, 0) += sh_weight * psh(0, j);
      d->Lx(1, 0) += sh_weight * psh(1, j);
      d->Lx(2, 0) += sh_weight * psh(2, j);
      d->Lx(3, 0) += sh_weight * pshoulder_0(1, j) * psh(2, j);
      d->Lx(4, 0) += -sh_weight * pshoulder_0(0, j) * psh(2, j);
      d->Lx(5, 0) += sh_weight * (-pshoulder_0(1, j) * psh(0, j) + pshoulder_0(0, j) * psh(1, j));

      d->Lxx(0, 0) += sh_weight;
      d->Lxx(1, 1) += sh_weight;
      d->Lxx(2, 2) += sh_weight;
      d->Lxx(3, 3) += sh_weight * pshoulder_0(1, j) * pshoulder_0(1, j);
      d->Lxx(3, 3) += sh_weight * pshoulder_0(0, j) * pshoulder_0(0, j);
      d->Lxx(5, 5) += sh_weight * (pshoulder_0(1, j) * pshoulder_0(1, j) + pshoulder_0(0, j) * pshoulder_0(0, j));

      d->Lxx(0, 5) += -sh_weight * pshoulder_0(1, j);
      d->Lxx(5, 0) += -sh_weight * pshoulder_0(1, j);

      d->Lxx(1, 5) += sh_weight * pshoulder_0(0, j);
      d->Lxx(5, 1) += sh_weight * pshoulder_0(0, j);

      d->Lxx(2, 3) += sh_weight * pshoulder_0(1, j);
      d->Lxx(2, 4) += -sh_weight * pshoulder_0(0, j);
      d->Lxx(3, 2) += sh_weight * pshoulder_0(1, j);
      d->Lxx(4, 2) += -sh_weight * pshoulder_0(0, j);

      d->Lxx(3, 4) += -sh_weight * pshoulder_0(1, j) * pshoulder_0(0, j);
      d->Lxx(4, 3) += -sh_weight * pshoulder_0(1, j) * pshoulder_0(0, j);
    }
  }

  // Cost derivative : Lu
  for (int i = 0; i < 4; i = i + 1) {
    r = friction_weight_ * rub_max_.segment(6 * i, 6);
    d->Lu.block(i * 3, 0, 3, 1) << r(0) - r(1), r(2) - r(3), -mu * (r(0) + r(1) + r(2) + r(3)) - r(4) + r(5);
  }
  d->Lu = d->Lu + (force_weights_.array() * d->r.template tail<12>().array()).matrix();

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
  d->Fx << A;

  for (int i = 0; i < 4; i = i + 1) {
    if (gait(i, 0) != 0) {
      forces_3d = u.block(3 * i, 0, 3, 1);
      d->Fx.block(9, 0, 3, 1) += -dt_ * I_inv * (base_vector_x.cross(forces_3d));
      d->Fx.block(9, 1, 3, 1) += -dt_ * I_inv * (base_vector_y.cross(forces_3d));
      d->Fx.block(9, 2, 3, 1) += -dt_ * I_inv * (base_vector_z.cross(forces_3d));
    }
  }
  d->Fu << B;
}

template <typename Scalar>
boost::shared_ptr<crocoddyl::ActionDataAbstractTpl<Scalar> > ActionModelQuadrupedNonLinearTpl<Scalar>::createData() {
  return boost::make_shared<ActionDataQuadrupedNonLinearTpl<Scalar> >(this);
}

////////////////////////////////
// get & set parameters ////////
////////////////////////////////

template <typename Scalar>
const typename Eigen::Matrix<Scalar, 12, 1>& ActionModelQuadrupedNonLinearTpl<Scalar>::get_force_weights() const {
  return force_weights_;
}
template <typename Scalar>
void ActionModelQuadrupedNonLinearTpl<Scalar>::set_force_weights(const typename MathBase::VectorXs& weights) {
  if (static_cast<std::size_t>(weights.size()) != state_->get_nx()) {
    throw_pretty("Invalid argument: "
                 << "Weights vector has wrong dimension (it should be " + std::to_string(state_->get_nx()) + ")");
  }
  force_weights_ = weights;
}

template <typename Scalar>
const typename Eigen::Matrix<Scalar, 12, 1>& ActionModelQuadrupedNonLinearTpl<Scalar>::get_state_weights() const {
  return state_weights_;
}
template <typename Scalar>
void ActionModelQuadrupedNonLinearTpl<Scalar>::set_state_weights(const typename MathBase::VectorXs& weights) {
  if (static_cast<std::size_t>(weights.size()) != state_->get_nx()) {
    throw_pretty("Invalid argument: "
                 << "Weights vector has wrong dimension (it should be " + std::to_string(state_->get_nx()) + ")");
  }
  state_weights_ = weights;
}

template <typename Scalar>
const Scalar& ActionModelQuadrupedNonLinearTpl<Scalar>::get_friction_weight() const {
  return friction_weight_;
}
template <typename Scalar>
void ActionModelQuadrupedNonLinearTpl<Scalar>::set_friction_weight(const Scalar& weight) {
  friction_weight_ = weight;
}

template <typename Scalar>
const Scalar& ActionModelQuadrupedNonLinearTpl<Scalar>::get_mu() const {
  return mu;
}
template <typename Scalar>
void ActionModelQuadrupedNonLinearTpl<Scalar>::set_mu(const Scalar& mu_coeff) {
  mu = mu_coeff;
}

template <typename Scalar>
const Scalar& ActionModelQuadrupedNonLinearTpl<Scalar>::get_mass() const {
  return mass;
}
template <typename Scalar>
void ActionModelQuadrupedNonLinearTpl<Scalar>::set_mass(const Scalar& m) {
  // The model need to be updated after this changed
  mass = m;
}

template <typename Scalar>
const Scalar& ActionModelQuadrupedNonLinearTpl<Scalar>::get_dt() const {
  return dt_;
}
template <typename Scalar>
void ActionModelQuadrupedNonLinearTpl<Scalar>::set_dt(const Scalar& dt) {
  // The model need to be updated after this changed
  dt_ = dt;
  g[8] = Scalar(-9.81) * dt_;
  A.topRightCorner(6, 6) << Eigen::Matrix<Scalar, 6, 6>::Identity() * dt_;
}

template <typename Scalar>
const typename Eigen::Matrix<Scalar, 3, 3>& ActionModelQuadrupedNonLinearTpl<Scalar>::get_gI() const {
  return gI;
}
template <typename Scalar>
void ActionModelQuadrupedNonLinearTpl<Scalar>::set_gI(const typename MathBase::Matrix3s& inertia_matrix) {
  // The model need to be updated after this changed
  if (static_cast<std::size_t>(inertia_matrix.size()) != 9) {
    throw_pretty("Invalid argument: "
                 << "gI has wrong dimension : 3x3");
  }
  gI = inertia_matrix;
}

template <typename Scalar>
const Scalar& ActionModelQuadrupedNonLinearTpl<Scalar>::get_min_fz_contact() const {
  // The model need to be updated after this changed
  return min_fz_in_contact;
}
template <typename Scalar>
void ActionModelQuadrupedNonLinearTpl<Scalar>::set_min_fz_contact(const Scalar& min_fz) {
  // The model need to be updated after this changed
  min_fz_in_contact = min_fz;
}

template <typename Scalar>
const Scalar& ActionModelQuadrupedNonLinearTpl<Scalar>::get_max_fz_contact() const {
  // The model need to be updated after this changed
  return max_fz;
}
template <typename Scalar>
void ActionModelQuadrupedNonLinearTpl<Scalar>::set_max_fz_contact(const Scalar& max_fz_) {
  // The model need to be updated after this changed
  max_fz = max_fz_;
  for (int i = 0; i < 4; i = i + 1) {
    ub(6 * i + 5) = max_fz;
  }
}

template <typename Scalar>
const Scalar& ActionModelQuadrupedNonLinearTpl<Scalar>::get_shoulder_hlim() const {
  return sh_hlim;
}
template <typename Scalar>
void ActionModelQuadrupedNonLinearTpl<Scalar>::set_shoulder_hlim(const Scalar& hlim) {
  // The model need to be updated after this changed
  sh_hlim = hlim;
}

template <typename Scalar>
const Scalar& ActionModelQuadrupedNonLinearTpl<Scalar>::get_shoulder_weight() const {
  return sh_weight;
}
template <typename Scalar>
void ActionModelQuadrupedNonLinearTpl<Scalar>::set_shoulder_weight(const Scalar& weight) {
  // The model need to be updated after this changed
  sh_weight = weight;
}

///////////////////////////
//// get A & B matrix /////
///////////////////////////
template <typename Scalar>
const typename Eigen::Matrix<Scalar, 12, 12>& ActionModelQuadrupedNonLinearTpl<Scalar>::get_A() const {
  return A;
}
template <typename Scalar>
const typename Eigen::Matrix<Scalar, 12, 12>& ActionModelQuadrupedNonLinearTpl<Scalar>::get_B() const {
  return B;
}

// to modify the cost on the command : || fz - m*g/nb contact ||^2
// --> set to True
template <typename Scalar>
const bool& ActionModelQuadrupedNonLinearTpl<Scalar>::get_relative_forces() const {
  return relative_forces;
}
template <typename Scalar>
void ActionModelQuadrupedNonLinearTpl<Scalar>::set_relative_forces(const bool& rel_forces) {
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
void ActionModelQuadrupedNonLinearTpl<Scalar>::update_model(
    const Eigen::Ref<const typename MathBase::MatrixXs>& l_feet,
    const Eigen::Ref<const typename MathBase::MatrixXs>& xref,
    const Eigen::Ref<const typename MathBase::MatrixXs>& S) {
  if (static_cast<std::size_t>(l_feet.size()) != 12) {
    throw_pretty("Invalid argument: "
                 << "l_feet matrix has wrong dimension (it should be : 3x4)");
  }
  if (static_cast<std::size_t>(xref.size()) != state_->get_nx()) {
    throw_pretty("Invalid argument: "
                 << "Weights vector has wrong dimension (it should be " + std::to_string(state_->get_nx()) + ")");
  }
  if (static_cast<std::size_t>(S.size()) != 4) {
    throw_pretty("Invalid argument: "
                 << "S vector has wrong dimension (it should be 4x1)");
  }

  xref_ = xref;
  gait = S;

  // Set ref u vector according to nb of contact
  uref_.setZero();
  if (relative_forces) {
    for (int i = 0; i < 4; i = i + 1) {
      if (gait[i] == 1) {
        uref_[3 * i + 2] = (Scalar(9.81) * mass) / (gait.sum());
      }
    }
  }

  R_tmp << cos(xref(5, 0)), -sin(xref(5, 0)), 0, sin(xref(5, 0)), cos(xref(5, 0)), 0, 0, 0, 1.0;

  I_inv = (R_tmp.transpose() * gI * R_tmp).inverse();  // I_inv
  lever_arms.block(0, 0, 2, 4) = l_feet.block(0, 0, 2, 4);

  for (int i = 0; i < 4; i = i + 1) {
    if (S(i, 0) != 0) {
      // set limit for normal force, (foot in contact with the ground)
      ub(6 * i + 4) = -min_fz_in_contact;

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
      ub(6 * i + 4) = 0.0;
      B.block(6, 3 * i, 3, 3).setZero();
      B.block(9, 3 * i, 3, 3).setZero();
    };
  };
}
}  // namespace quadruped_walkgen

#endif
