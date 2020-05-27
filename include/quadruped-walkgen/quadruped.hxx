#ifndef __quadruped_walkgen_quadruped_hxx__
#define __quadruped_walkgen_quadruped_hxx__

#include "crocoddyl/core/utils/exception.hpp"


namespace quadruped_walkgen  {
template <typename Scalar>
ActionModelQuadrupedTpl<Scalar>::ActionModelQuadrupedTpl()
    : crocoddyl::ActionModelAbstractTpl<Scalar>(boost::make_shared<crocoddyl::StateVectorTpl<Scalar> >(12), 12, 24) ,
    cone(Eigen::Matrix<Scalar, 3, 1>(0,0,1) , 1, 4 , false, 0., std::numeric_limits<Scalar>::max())
  {
  mu = 0.8 ; 
  dt_ = 0.02 ; 
  mass = 2.97784899 ; 

  // Weights initialization
  force_weights_ << Eigen::Matrix<Scalar, 12, 1>::Constant(12,1,0.1);
  state_weights_ << 1., 1.,150.,35.,30.,8.,20.,20.,15.,4.,4.,8.  ; 
  friction_weight_ = 1 ;

  // Matrix initialization
  g << 0.,0.,0. ,0.,0.,0., 0.,0.,-9.81*dt_ ,0.,0.,0. ; 
  gI << 0.00578574,0.,0.,
        0., 0.01938108, 0.,
        0.,0.,0.02476124;
  A << Eigen::Matrix<Scalar, 12, 12>::Identity() ; 
  A.block(0,6,6,6) << Eigen::Matrix<Scalar, 6, 6>::Identity()*dt_ ; 
  B << Eigen::Matrix<Scalar, 12, 12>::Zero() ; 
  lever_arms << Eigen::Matrix<Scalar, 3, 4>::Zero() ; 
  R << Eigen::Matrix<Scalar, 3, 3>::Zero() ; 
  lever_tmp << 0.,0.,0. ;
  R_tmp << Eigen::Matrix<Scalar, 3, 3>::Zero() ; 
  
  // Cone initialization
  nsurf << 0.,0.,1 ; 
  cone.update(nsurf , mu , false  ) ; 
  lb << cone.get_lb() , cone.get_lb() , cone.get_lb() , cone.get_lb() ;
  ub << cone.get_ub() , cone.get_ub() , cone.get_ub() , cone.get_ub() ;
  
  // Matrix (20x12) to deal with 4 friction cone at same time
  Fa << Eigen::Matrix<Scalar, 20, 12>::Zero() ; 
  for (int i=0; i<4; i=i+1){
    Fa.block(i*5,i*3 , 5,3) = cone.get_A() ; 
    
  }

  Fa_x_u << Eigen::Matrix<Scalar, 20, 1>::Zero() ; 
  rlb_ << Eigen::Matrix<Scalar, 20, 1>::Zero() ; 
  rub_ << Eigen::Matrix<Scalar, 20, 1>::Zero() ; 
  rlb_min_ << Eigen::Matrix<Scalar, 20, 1>::Zero() ; 
  rub_max_ << Eigen::Matrix<Scalar, 20, 1>::Zero() ;
  Arr << Eigen::Matrix<Scalar, 20, 20>::Zero() ; 
}


template <typename Scalar>
ActionModelQuadrupedTpl<Scalar>::~ActionModelQuadrupedTpl() {}

template <typename Scalar>
void ActionModelQuadrupedTpl<Scalar>::calc(const boost::shared_ptr<crocoddyl::ActionDataAbstractTpl<Scalar> >& data,
                                          const Eigen::Ref<const typename MathBase::VectorXs>& x,
                                          const Eigen::Ref<const typename MathBase::VectorXs>& u) {
  if (static_cast<std::size_t>(x.size()) != state_->get_nx()) {
    throw_pretty("Invalid argument: "
                 << "x has wrong dimension (it should be " + std::to_string(state_->get_nx()) + ")");
  }
  if (static_cast<std::size_t>(u.size()) != nu_) {
    throw_pretty("Invalid argument: "
                 << "u has wrong dimension (it should be " + std::to_string(nu_) + ")");
  }

  ActionDataQuadrupedTpl<Scalar>* d = static_cast<ActionDataQuadrupedTpl<Scalar>*>(data.get());
 
  // // Discrete dynamic
  d->xnext << A.diagonal().cwiseProduct(x) + B*u + g;

  // Residual cost on the state and force norm
  d->r.template head<12>() =  state_weights_.asDiagonal() * (x - xref_);
  d->r.template tail<12>() =  force_weights_.asDiagonal() * u;

  for (int i=0; i<4; i=i+1){
     Fa_x_u.segment(5*i,5) << u(3*i) - mu*u(3*i+2) , -u(3*i) - mu*u(3*i+2),
                              u(3*i+1) - mu*u(3*i+2) , -u(3*i+1) - mu*u(3*i+2),
                              u(3*i+2) ;                             
  }
  rlb_min_ = (Fa_x_u - lb).array().min(0.);
  rub_max_ = (Fa_x_u - ub).array().max(0.);   
  
  // Cost computation 
  d->cost = 0.5 * d->r.transpose() * d->r     + friction_weight_ * (Scalar(0.5) * rlb_min_.matrix().squaredNorm() +
                                                Scalar(0.5) * rub_max_.matrix().squaredNorm()) ;
}







template <typename Scalar>
void ActionModelQuadrupedTpl<Scalar>::calcDiff(const boost::shared_ptr<crocoddyl::ActionDataAbstractTpl<Scalar> >& data,
                                              const Eigen::Ref<const typename MathBase::VectorXs>& x,
                                              const Eigen::Ref<const typename MathBase::VectorXs>& u) {
  if (static_cast<std::size_t>(x.size()) != state_->get_nx()) {
    throw_pretty("Invalid argument: "
                 << "x has wrong dimension (it should be " + std::to_string(state_->get_nx()) + ")");
  }
  if (static_cast<std::size_t>(u.size()) != nu_) {
    throw_pretty("Invalid argument: "
                 << "u has wrong dimension (it should be " + std::to_string(nu_) + ")");
  }

  ActionDataQuadrupedTpl<Scalar>* d = static_cast<ActionDataQuadrupedTpl<Scalar>*>(data.get());

  // // Matrix friction cone hessian
  rlb_ = Fa_x_u - lb ; 
  rub_ = Fa_x_u - ub ; 
  Arr.diagonal() =
        ((rlb_.array() <= 0.) + (rub_.array() >= 0.)).matrix().template cast<Scalar>() ; 
  
  // Cost derivatives
  d->Lx = (state_weights_.array()* d->r.template head<12>().array()).matrix() ;

  Scalar r1 = friction_weight_*(rlb_min_(0) + rub_max_(0)) ; 
  Scalar r2 = friction_weight_*(rlb_min_(1) + rub_max_(1)) ; 
  Scalar r3 = friction_weight_*(rlb_min_(2) + rub_max_(2)) ; 
  Scalar r4 = friction_weight_*(rlb_min_(3) + rub_max_(3)) ; 
  Scalar r5 = friction_weight_*(rlb_min_(4) + rub_max_(4)) ; 
  d->Lu.block(0,0,3,1) << r1 - r2 , r3 - r4 , -mu*(r1 + r2 + r3 + r4 ) + r5 ;  
  for (int i=1; i<4; i=i+1){
    r1 = friction_weight_*(rlb_min_(5*i) + rub_max_(5*i)) ; 
    r2 = friction_weight_*(rlb_min_(5*i+1) + rub_max_(5*i+1)) ; 
    r3 = friction_weight_*(rlb_min_(5*i+2) + rub_max_(5*i+2)) ; 
    r4 = friction_weight_*(rlb_min_(5*i+3) + rub_max_(5*i+3)) ; 
    r5 = friction_weight_*(rlb_min_(5*i+4) + rub_max_(5*i+4)) ; 
    d->Lu.block(i*3,0,3,1) << r1 - r2 , r3 - r4 , -mu*(r1 + r2 + r3 + r4 ) + r5 ; 
  } 
  d->Lu += (force_weights_.array()*d->r.template tail<12>().array()).matrix() ; 
  
  d->Lxx.diagonal() = (state_weights_.array() * state_weights_.array()).matrix() ;  
  
  for (int i=0; i<4; i=i+1){
      r1 = friction_weight_*Arr(5*i,5*i) ; 
      r2 = friction_weight_*Arr(5*i+1,5*i+1) ; 
      r3 = friction_weight_*Arr(5*i+2,5*i+2) ; 
      r4 = friction_weight_*Arr(5*i+3,5*i+3) ; 
      r5 = friction_weight_*Arr(5*i+4,5*i+4) ; 
      d->Luu.block(3*i,3*i,3,3) << r1 + r2 , 0.0 , mu*(r2 - r1 ),
                                    0.0,  r3 + r4 , mu*(r4 - r3 ), 
                                  mu*(r2 - r1 ) , mu*(r4 - r3) , mu*mu*(r1 + r2 + r3 + r4) + r5  ; 
  }
  d->Luu.diagonal() = d->Luu.diagonal() + (force_weights_.array() * force_weights_.array()).matrix() ;

  // Dynamic derivatives
  d->Fx << A;
  d->Fu << B;  
}



template <typename Scalar>
boost::shared_ptr<crocoddyl::ActionDataAbstractTpl<Scalar> > ActionModelQuadrupedTpl<Scalar>::createData() {
  return boost::make_shared<ActionDataQuadrupedTpl<Scalar> >(this);
}

///////////////////////////////////////
// get and set weights vectors ////////
///////////////////////////////////////

template <typename Scalar>
const typename Eigen::Matrix<Scalar, 12, 1>& ActionModelQuadrupedTpl<Scalar>::get_force_weights() const {
  return force_weights_;
}
template <typename Scalar>
void ActionModelQuadrupedTpl<Scalar>::set_force_weights(const typename MathBase::VectorXs& weights) {
  if (static_cast<std::size_t>(weights.size()) != state_->get_nx()) {
    throw_pretty("Invalid argument: "
                 << "Weights vector has wrong dimension (it should be " + std::to_string(state_->get_nx()) + ")");
  }  
  force_weights_ = weights;
}

template <typename Scalar>
const typename Eigen::Matrix<Scalar, 12, 1>& ActionModelQuadrupedTpl<Scalar>::get_state_weights() const {
  return state_weights_;
}
template <typename Scalar>
void ActionModelQuadrupedTpl<Scalar>::set_state_weights(const typename MathBase::VectorXs& weights) {
  if (static_cast<std::size_t>(weights.size()) != state_->get_nx()) {
    throw_pretty("Invalid argument: "
                 << "Weights vector has wrong dimension (it should be " + std::to_string(state_->get_nx()) + ")");
  }
  state_weights_ = weights;
}

template <typename Scalar>
const Scalar& ActionModelQuadrupedTpl<Scalar>::get_friction_weight() const {
  return friction_weight_;
}

template <typename Scalar>
void ActionModelQuadrupedTpl<Scalar>::set_friction_weight(const Scalar& weight) {
  friction_weight_ = weight;
}


///////////////////////////
//// get A & B matrix /////
///////////////////////////
template <typename Scalar>
const typename Eigen::Matrix<Scalar, 12, 12>& ActionModelQuadrupedTpl<Scalar>::get_A() const {
  return A;
}
template <typename Scalar>
const typename Eigen::Matrix<Scalar, 12, 12>& ActionModelQuadrupedTpl<Scalar>::get_B() const {
  return B;
}

template<typename Scalar>
const typename Eigen::Matrix<Scalar, 3, 3> get_skew(
    const typename Eigen::Matrix<Scalar, 3, 1 >& vec) {

  return (Eigen::Matrix<Scalar, 3, 3>() << 0.0, -vec[2], vec[1],
      vec[2], 0.0, -vec[0], -vec[1], vec[0], 0.0);
}

////////////////////////
// Update current model 
////////////////////////

template <typename Scalar>
void ActionModelQuadrupedTpl<Scalar>::update_model(const Eigen::Ref<const typename MathBase::MatrixXs>& l_feet  ,
                    const Eigen::Ref<const typename MathBase::MatrixXs>& xref,
                    const Eigen::Ref<const typename MathBase::MatrixXs>& S ) {
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

  xref_ = xref ; 

  R_tmp << cos(xref(5,0)),-sin(xref(5,0)),0,
      sin(xref(5,0)),cos(xref(5,0)),0,
      0,0,1.0 ; 
  
  R = (R_tmp*gI).inverse() ; // I_inv  
  lever_arms.block(0,0,2,4) = l_feet.block(0,0,2,4) ; 

  for (int i=0; i<4; i=i+1){
    if (S(i,0) != 0) {
      B.block(6 , 3*i  , 3,3).diagonal() << dt_ / mass  ,  dt_ / mass  ,  dt_ / mass  ; 
      lever_tmp = lever_arms.block(0,i,3,1) - xref.block(0,0,3,1) ;
      R_tmp << 0.0, -lever_tmp[2], lever_tmp[1],
      lever_tmp[2], 0.0, -lever_tmp[0], -lever_tmp[1], lever_tmp[0], 0.0 ; 
      B.block(9 , 3*i  , 3,3) << dt_ * R* R_tmp; 
    }
    else{
      B.block(6 , 3*i  , 3,3) << Eigen::Matrix<Scalar, 3, 3>::Zero();
      B.block(9 , 3*i  , 3,3) << Eigen::Matrix<Scalar, 3, 3>::Zero();
    };
  } ; 
}
}

#endif
