#ifndef __quadruped_walkgen_quadruped_step_hxx__
#define __quadruped_walkgen_quadruped_step_hxx__

#include "crocoddyl/core/utils/exception.hpp"


namespace quadruped_walkgen  {
template <typename Scalar>
ActionModelQuadrupedStepTpl<Scalar>::ActionModelQuadrupedStepTpl()
    : crocoddyl::ActionModelAbstractTpl<Scalar>(boost::make_shared<crocoddyl::StateVectorTpl<Scalar> >(20), 4, 24)
  { 

  B.setZero() ; 
  state_weights_ << Scalar(1.)  , Scalar(1.) , Scalar(150.) , Scalar(35.),
                    Scalar(30.) , Scalar(8.) , Scalar(20.)  , Scalar(20.) , 
                    Scalar(15.) , Scalar(4.) , Scalar(4.)   , Scalar(8.)  ; 
  shoulder_weights_.setConstant(Scalar(1)) ; 
  pshoulder_ <<  Scalar(0.1946) ,  Scalar(0.15005),  Scalar(0.1946) ,  Scalar(-0.15005) ,
                 Scalar(-0.1946),  Scalar(0.15005) , Scalar(-0.1946),  Scalar(-0.15005) ; 
  pshoulder_0 <<  Scalar(0.1946) ,   Scalar(0.1946) ,   Scalar(-0.1946),  Scalar(-0.1946) , 
                  Scalar(0.15005) ,  Scalar(-0.15005)  , Scalar(0.15005)  ,  Scalar(-0.15005) ; 
  pshoulder_tmp.setZero() ; 
  pcentrifugal_tmp_1.setZero() ; 
  pcentrifugal_tmp_2.setZero() ; 
  pcentrifugal_tmp.setZero() ; 
  centrifugal_term = true ; 
  symmetry_term = true ; 
  T_gait = Scalar(0.32) ; 
  
  step_weights_.setConstant(Scalar(1)) ;
  
 
}


template <typename Scalar>
ActionModelQuadrupedStepTpl<Scalar>::~ActionModelQuadrupedStepTpl() {}


template <typename Scalar>
void ActionModelQuadrupedStepTpl<Scalar>::calc(const boost::shared_ptr<crocoddyl::ActionDataAbstractTpl<Scalar> >& data,
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

  ActionDataQuadrupedStepTpl<Scalar>* d = static_cast<ActionDataQuadrupedStepTpl<Scalar>*>(data.get());
 
  d->xnext.template head<12>() = x.head(12) ; 
  d->xnext.template tail<8>() = x.tail(8) + B*u;
  
  // Residual cost on the state and force norm
  d->r.template head<12>() =  state_weights_.cwiseProduct(x.head(12) - xref_);
  d->r.template segment<8>(12) = shoulder_weights_.cwiseProduct(x.tail(8) - pshoulder_); 
  d->r.template tail<4>() =  step_weights_.cwiseProduct(u);

  d->cost = Scalar(0.5) * d->r.transpose() * d->r   ;
}


template <typename Scalar>
void ActionModelQuadrupedStepTpl<Scalar>::calcDiff(const boost::shared_ptr<crocoddyl::ActionDataAbstractTpl<Scalar> >& data,
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

  ActionDataQuadrupedStepTpl<Scalar>* d = static_cast<ActionDataQuadrupedStepTpl<Scalar>*>(data.get());  
  
  // Cost derivatives : Lx
  d->Lx.template head<12>() = (state_weights_.array()* d->r.template head<12>().array()).matrix() ;
  d->Lx.template tail<8>() = (shoulder_weights_.array()* d->r.template segment<8>(12).array()).matrix() ;

  d->Lu =(step_weights_.array()*d->r.template tail<4>().array()).matrix() ; 
  
  // Hessian : Lxx
  d->Lxx.diagonal().head(12) = (state_weights_.array() * state_weights_.array()).matrix() ;  
  d->Lxx.diagonal().tail(8) = (shoulder_weights_.array() * shoulder_weights_.array()).matrix() ;  
  
  d->Luu.diagonal() = (step_weights_.array() * step_weights_.array()).matrix() ;

  // Dynamic derivatives
  d->Fx.setIdentity();
  d->Fu.block(12,0,8,4) = B;  
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
  if (static_cast<std::size_t>(weights.size()) != 12 ) {
    throw_pretty("Invalid argument: "
                 << "Weights vector has wrong dimension (it should be 12)");
  }
  state_weights_ = weights;
}

template <typename Scalar>
const typename Eigen::Matrix<Scalar, 4, 1>& ActionModelQuadrupedStepTpl<Scalar>::get_step_weights() const {
  return step_weights_;
}
template <typename Scalar>
void ActionModelQuadrupedStepTpl<Scalar>::set_step_weights(const typename MathBase::VectorXs& weights) {
  if (static_cast<std::size_t>(weights.size()) != 4 ) {
    throw_pretty("Invalid argument: "
                 << "Weights vector has wrong dimension (it should be 4)");
  }
  step_weights_ = weights;
}

template <typename Scalar>
const typename Eigen::Matrix<Scalar, 8, 1>& ActionModelQuadrupedStepTpl<Scalar>::get_shoulder_weights() const {
  return shoulder_weights_;
}
template <typename Scalar>
void ActionModelQuadrupedStepTpl<Scalar>::set_shoulder_weights(const typename MathBase::VectorXs& weights) {
  if (static_cast<std::size_t>(weights.size()) != 8 ) {
    throw_pretty("Invalid argument: "
                 << "Weights vector has wrong dimension (it should be 8)");
  }
  shoulder_weights_ = weights;
}

template <typename Scalar>
const typename Eigen::Matrix<Scalar, 8, 1>& ActionModelQuadrupedStepTpl<Scalar>::get_shoulder_position() const {
  return pshoulder_ ;
}
template <typename Scalar>
void ActionModelQuadrupedStepTpl<Scalar>::set_shoulder_position(const typename MathBase::VectorXs& pos) {
  if (static_cast<std::size_t>(pos.size()) != 8 ) {
    throw_pretty("Invalid argument: "
                 << "Weights vector has wrong dimension (it should be 8)");
  }
  pshoulder_ = pos;
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


////////////////////////
// Update current model 
////////////////////////

template <typename Scalar>
void ActionModelQuadrupedStepTpl<Scalar>::update_model(const Eigen::Ref<const typename MathBase::MatrixXs>& l_feet  ,
                    const Eigen::Ref<const typename MathBase::MatrixXs>& xref,
                    const Eigen::Ref<const typename MathBase::MatrixXs>& S ) {
  if (static_cast<std::size_t>(l_feet.size()) != 12) {
    throw_pretty("Invalid argument: "
                 << "l_feet matrix has wrong dimension (it should be : 3x4)");
  }
  if (static_cast<std::size_t>(xref.size()) !=  12 ) {
    throw_pretty("Invalid argument: "
                 << "Weights vector has wrong dimension (it should be " + std::to_string(state_->get_nx()) + ")");
  }
  if (static_cast<std::size_t>(S.size()) != 4) {
    throw_pretty("Invalid argument: "
                 << "S vector has wrong dimension (it should be 4x1)");
  }

  xref_ = xref ; 

  R_tmp << cos(xref(5,0)) ,-sin(xref(5,0)) , Scalar(0),
      sin(xref(5,0)), cos(xref(5,0)), Scalar(0),
      Scalar(0),Scalar(0),Scalar(1.0) ; 

   // Centrifual term 
  pcentrifugal_tmp_1 = xref.block(6,0,3,1) ; 
  pcentrifugal_tmp_2 = xref.block(9,0,3,1) ; 
  pcentrifugal_tmp = 0.5*std::sqrt(xref(2,0)/9.81) * pcentrifugal_tmp_1.cross(pcentrifugal_tmp_2) ; 
  

  for (int i=0; i<4; i=i+1){
    pshoulder_tmp.block(0,i,2,1) =  R_tmp.block(0,0,2,2)*(pshoulder_0.block(0,i,2,1) +   symmetry_term * 0.25*T_gait*xref.block(6,0,2,1) + centrifugal_term * pcentrifugal_tmp.block(0,0,2,1) );
    pshoulder_[2*i] = pshoulder_tmp(0,i) +   xref(0,0); 
    pshoulder_[2*i+1] = pshoulder_tmp(1,i) +  xref(1,0); 
  }

  B.setZero() ; 

  if (S(0,0) == Scalar(1)) {
    B.block(0,0,2,2).setIdentity() ; 
    B.block(6,2,2,2).setIdentity() ;    
  }
  else {
    B.block(2,0,2,2).setIdentity() ; 
    B.block(4,2,2,2).setIdentity() ;  

  }
      
 
}
}

#endif
