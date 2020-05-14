#ifndef __quadruped_walkgen_quadruped_hxx__
#define __quadruped_walkgen_quadruped_hxx__

#include "crocoddyl/core/utils/exception.hpp"


namespace quadruped_walkgen  {
template <typename Scalar>
ActionModelQuadrupedTpl<Scalar>::ActionModelQuadrupedTpl()
    : crocoddyl::ActionModelAbstractTpl<Scalar>(boost::make_shared<crocoddyl::StateVectorTpl<Scalar> >(12), 12, 24) ,
    activation(crocoddyl::ActivationBoundsTpl<Scalar>() ) ,
    dataCost(crocoddyl::ActivationModelAbstractTpl<Scalar>(20)) // Problem 
  {
  mu = 0.8 ; 
  dt_ = 0.02 ; 
  mass = 2.97784899 ; 
  cost_weights_ << 10. , 1.;
  // Weights
  force_weights_ << Eigen::Matrix<Scalar, 12, 1>::Constant(12,1,0.1);
  state_weights_ << 1., 1.,150.,35.,30.,8.,20.,20.,15.,4.,4.,8.  ; 

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
  xref_ <<  Eigen::Matrix<Scalar, 12, 1>::Zero() ; 
  
  nsurf << 0.,0.,1 ; 
  cone.update(nsurf , mu , 4 , true ) ; 
  lb << cone.get_lb() , cone.get_lb() , cone.get_lb() , cone.get_lb() ;
  ub << cone.get_ub() , cone.get_ub() , cone.get_ub() , cone.get_ub() ;
  
  Fa << Eigen::Matrix<Scalar, 20, 12>::Zero() ; 
  for (int i=0; i<4; i=i+1){
    Fa.block(i*5,i*3 , 5,3) = cone.get_A() ; 
  }

  activation.set_bounds(crocoddyl::ActivationBoundsTpl<Scalar>(lb,ub) ) ; 
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
 
  d->xnext << A*x + B*u + g;
  d->r.template head<12>() =  state_weights_.asDiagonal() * (x - xref_);
  d->r.template tail<12>() =  force_weights_.asDiagonal() * u;
  d->cost = 0.5 * d->r.transpose() * d->r;
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

  // Cost derivatives
  d->Lx = d->r.template head<12>() ;
  d->Lu = d->r.template tail<12>();
  d->Lxx.diagonal() << state_weights_.array() * state_weights_.array() ;
  d->Luu.diagonal() << force_weights_.array() * force_weights_.array();

  // Dynamic derivatives
  d->Fx << A;
  d->Fu << B;
}






template <typename Scalar>
boost::shared_ptr<crocoddyl::ActionDataAbstractTpl<Scalar> > ActionModelQuadrupedTpl<Scalar>::createData() {
  return boost::make_shared<ActionDataQuadrupedTpl<Scalar> >(this);
}

///////////////////////////////
// get and set weights ////////
///////////////////////////////
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
///////////////////////////
//// get A & B
/////////////////////////
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

template <typename Scalar>
void ActionModelQuadrupedTpl<Scalar>::update_model(const Eigen::Ref<const typename MathBase::MatrixXs>& l_feet  ,
                    const Eigen::Ref<const typename MathBase::VectorXs>& xref,
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
  
  R << cos(xref(5)),-sin(xref(5)),0,
      sin(xref(5)),cos(xref(5)),0,
      0,0,1.0 ; 
  
  R.noalias() = (R*gI).inverse() ; // I_inv  
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
