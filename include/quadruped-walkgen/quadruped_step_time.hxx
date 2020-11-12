#ifndef __quadruped_walkgen_quadruped_step_time_hxx__
#define __quadruped_walkgen_quadruped_step_time_hxx__

#include "crocoddyl/core/utils/exception.hpp"


namespace quadruped_walkgen  {
template <typename Scalar>
ActionModelQuadrupedStepTimeTpl<Scalar>::ActionModelQuadrupedStepTimeTpl()
    : crocoddyl::ActionModelAbstractTpl<Scalar>(boost::make_shared<crocoddyl::StateVectorTpl<Scalar> >(21), 4, 25)
  { 

  B.setZero() ; 
  rub_max_.setZero() ; 
  rub_max_bool.setZero() ;

  state_weights_ << Scalar(1.)  , Scalar(1.) , Scalar(150.) , Scalar(35.),
                    Scalar(30.) , Scalar(8.) , Scalar(20.)  , Scalar(20.) , 
                    Scalar(15.) , Scalar(4.) , Scalar(4.)   , Scalar(8.)  ; 
  heuristicWeights.setConstant(Scalar(1)) ; 
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
  T_gait = Scalar(0.64) ; 
  
  step_weights_.setConstant(Scalar(1)) ;

  // Optim dt
  nb_nodes = Scalar(15.) ; 
  vlim = Scalar(2.) ;
  beta_lim = Scalar((64*nb_nodes*nb_nodes*vlim*vlim)/225) ; // apparent speed used in the cost function
  speed_weight = Scalar(10.) ;   

  // Log cost
  cost_.setZero() ; 
  log_cost = true ; 

  // indicates whether it t the 1st step, otherwise the cost function is much simpler (acc, speed = 0)
  first_step = false ; 


  alpha.setLinSpaced(3, Scalar(0.0), Scalar(1.0)) ; 
  alpha2.setZero() ;
  alpha2.col(0) << alpha ;
  alpha2.col(1) << alpha.pow(2) ; 
  alpha2.col(2) << alpha.pow(3) ; 
  alpha2.col(3) << alpha.pow(4) ; 

  b_coeff.setZero();
  b_coeff.col(0) = Scalar(1.0) - Scalar(18.)*alpha2.col(1) + Scalar(32.)*alpha2.col(2) - Scalar(15.)*alpha2.col(3) ; 
  b_coeff.col(1) = alpha2.col(0) - Scalar(4.5)*alpha2.col(1) + Scalar(6.)*alpha2.col(2) - Scalar(2.5)*alpha2.col(3) ; 
  b_coeff.col(2) = Scalar(30.)*alpha2.col(1) - Scalar(60.)*alpha2.col(2) + Scalar(30.)*alpha2.col(3) ; 

  b_coeff2.setZero() ;

  lfeet.setZero() ;
  rub_max_bool_first.setZero() ; 
  rub_max_first.setZero() ; 
  rub_max_first_2.setZero() ;
  
 
}


template <typename Scalar>
ActionModelQuadrupedStepTimeTpl<Scalar>::~ActionModelQuadrupedStepTimeTpl() {}


template <typename Scalar>
void ActionModelQuadrupedStepTimeTpl<Scalar>::calc(const boost::shared_ptr<crocoddyl::ActionDataAbstractTpl<Scalar> >& data,
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

  ActionDataQuadrupedStepTimeTpl<Scalar>* d = static_cast<ActionDataQuadrupedStepTimeTpl<Scalar>*>(data.get());
 
  d->xnext.template head<12>() = x.head(12) ; 
  d->xnext.template segment<8>(12) = x.segment(12,8) + B*u;
  d->xnext.template tail<1>() = x.tail(1) ; 
  
  // Residual cost on the state and force norm
  d->r.template head<12>() =  state_weights_.cwiseProduct(x.head(12) - xref_);
  d->r.template segment<8>(12) = heuristicWeights.cwiseProduct(x.segment(12,8) - pshoulder_); 
  d->r.template tail<4>() =  step_weights_.cwiseProduct(u);

  if (first_step){
    for (int i = 0 ; i < 4 ; i++){
      rub_max_first.col(i) << x(20)*b_coeff2.col(3*i) + x(20)*x(20)*b_coeff2.col(3*i+1) + u(i)*b_coeff2.col(3*i+2) + lfeet(0,i)*b_coeff2.col(3*i+2);
    }

    rub_max_first_2.col(0) << rub_max_first.col(0).pow(2) + rub_max_first.col(1).pow(2) - x(20)*x(20)*vlim*vlim*nb_nodes*nb_nodes ; 
    rub_max_first_2.col(1) << rub_max_first.col(2).pow(2) + rub_max_first.col(3).pow(2) - x(20)*x(20)*vlim*vlim*nb_nodes*nb_nodes ; 
    
    rub_max_bool_first = (rub_max_first_2 > Scalar(0.)).template cast<Scalar>() ;  
    rub_max_first_2 = rub_max_first_2.cwiseMax(Scalar(0.)) ; 

    d->cost = Scalar(0.5) * d->r.transpose() * d->r  ;
    for (int i=0 ; i<3 ; i++){
      d->cost +=  speed_weight * Scalar(0.5) * rub_max_first_2.row(i).sum() ;
    }
  }
  else{
    rub_max_ <<  u[0]*u[0] + u[1]*u[1] - beta_lim*x[20]*x[20] , u[2]*u[2] + u[3]*u[3] - beta_lim*x[20]*x[20] ;
  
    rub_max_bool = (rub_max_.array() >= Scalar(0.)).matrix().template cast<Scalar>() ; 
    rub_max_ = rub_max_.cwiseMax(Scalar(0.)) ; 

    d->cost = Scalar(0.5) * d->r.transpose() * d->r  + speed_weight * Scalar(0.5) * rub_max_.sum();
  }
  

  if (log_cost){
    cost_[3] = 0 ; 
    // Length to be consistent with others models
    cost_[0] = Scalar(0.5)*d->r.head(12).transpose()*d->r.head(12) ; // State cost
    cost_[1] = Scalar(0.5)*d->r.segment(12,8).transpose()*d->r.segment(12,8) ; // heuristic cost
    cost_[2] = Scalar(0.5)*d->r.tail(4).transpose()*d->r.tail(4) ;  // Delta feet cost

    if (first_step){
      for (int i=0 ; i<3 ; i++){
        cost_[3] +=  speed_weight * Scalar(0.5) * rub_max_first_2.row(i).sum() ;
      }
    }
    else{
      cost_[3] = speed_weight * Scalar(0.5) * rub_max_.sum() ;
    }
  }
}


template <typename Scalar>
void ActionModelQuadrupedStepTimeTpl<Scalar>::calcDiff(const boost::shared_ptr<crocoddyl::ActionDataAbstractTpl<Scalar> >& data,
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

  ActionDataQuadrupedStepTimeTpl<Scalar>* d = static_cast<ActionDataQuadrupedStepTimeTpl<Scalar>*>(data.get());  
  
  d->Lx.setZero() ; 
  d->Lu.setZero() ; 
  d->Lxu.setZero() ; 
  d->Lxx.setZero() ; 
  d->Luu.setZero() ; 
  // Cost derivatives : Lx
  d->Lx.template head<12>() = (state_weights_.array()* d->r.template head<12>().array()).matrix() ;
  d->Lx.template segment<8>(12) = (heuristicWeights.array()* d->r.template segment<8>(12).array()).matrix() ;

  if (first_step){
    for (int i=0 ; i <3 ; i++){
      // First foot
      if (rub_max_bool_first(i,0)){
        d->Lx(20) += speed_weight*( b_coeff2(i,0) + Scalar(2)*x(20)*b_coeff2(i,1))*rub_max_first(i,0) +
                     speed_weight*( b_coeff2(i,3) + Scalar(2)*x(20)*b_coeff2(i,4))*rub_max_first(i,1) -
                     speed_weight*x(20)*vlim*vlim*nb_nodes*nb_nodes;
        d->Lu(0) += speed_weight*b_coeff2(i,2)*rub_max_first(i,0)  ;  
        d->Lu(1) += speed_weight*b_coeff2(i,5)*rub_max_first(i,1)  ;  

        d->Luu(0,0) += speed_weight*b_coeff2(i,2)*b_coeff2(i,2)  ; 
        d->Luu(1,1) += speed_weight*b_coeff2(i,5)*b_coeff2(i,5)  ; 
        d->Lxu(20,0) += speed_weight*( b_coeff2(i,0) + Scalar(2)*x(20)*b_coeff2(i,1))*b_coeff2(i,2) ;
        d->Lxu(20,1) += speed_weight*( b_coeff2(i,3) + Scalar(2)*x(20)*b_coeff2(i,4))*b_coeff2(i,5) ; 
        d->Lxx(20,20) += speed_weight*std::pow( b_coeff2(i,0) + Scalar(2)*x(20)*b_coeff2(i,1) , 2) + 
                        speed_weight*Scalar(2)*b_coeff2(i,1)*rub_max_first(i,0) +
                        speed_weight*std::pow( b_coeff2(i,3) + Scalar(2)*x(20)*b_coeff2(i,4) , 2) + 
                        speed_weight*Scalar(2)*b_coeff2(i,4)*rub_max_first(i,1) -
                        speed_weight*vlim*vlim*nb_nodes*nb_nodes ; 

      } 
    
      //  Second foot
      if (rub_max_bool_first(i,1)){
        d->Lx(20) += speed_weight*( b_coeff2(i,6) + Scalar(2)*x(20)*b_coeff2(i,7))*rub_max_first(i,2) +
                        speed_weight*( b_coeff2(i,9) + Scalar(2)*x(20)*b_coeff2(i,10))*rub_max_first(i,3) -
                        speed_weight*x(20)*vlim*vlim*nb_nodes*nb_nodes;
        d->Lu(2) += speed_weight*b_coeff2(i,8)*rub_max_first(i,2)  ;  
        d->Lu(3) += speed_weight*b_coeff2(i,11)*rub_max_first(i,3)  ;  

        d->Luu(2,2) += speed_weight*b_coeff2(i,8)*b_coeff2(i,8)  ; 
        d->Luu(3,3) += speed_weight*b_coeff2(i,11)*b_coeff2(i,11)  ; 
        d->Lxu(20,2) += speed_weight*( b_coeff2(i,6) + Scalar(2)*x(20)*b_coeff2(i,7))*b_coeff2(i,8) ; 
        d->Lxu(20,3) += speed_weight*( b_coeff2(i,9) + Scalar(2)*x(20)*b_coeff2(i,10))*b_coeff2(i,11) ;
        d->Lxx(20,20) += speed_weight*std::pow( b_coeff2(i,6) + Scalar(2)*x(20)*b_coeff2(i,7) , 2) + 
                        speed_weight*Scalar(2)*b_coeff2(i,7)*rub_max_first(i,2) +
                        speed_weight*std::pow( b_coeff2(i,9) + Scalar(2)*x(20)*b_coeff2(i,10) , 2) + 
                        speed_weight*Scalar(2)*b_coeff2(i,10)*rub_max_first(i,3) -
                        speed_weight*vlim*vlim*nb_nodes*nb_nodes ; 
      }
    }
  }
  else{
    d->Lx.template tail<1>() << - beta_lim*speed_weight*x(20)*rub_max_bool[0] - beta_lim*speed_weight*x(20)*rub_max_bool[1] ;
 
    d->Lu << speed_weight*u[0]*rub_max_bool[0] , speed_weight*u[1]*rub_max_bool[0] , 
            speed_weight*u[2]*rub_max_bool[1],  speed_weight*u[3]*rub_max_bool[1] ;

    d->Lxx(20,20) = - beta_lim*speed_weight*rub_max_bool[0] - beta_lim*speed_weight*rub_max_bool[1] ; 
  

    d->Luu.diagonal() << speed_weight*rub_max_bool[0] ,
                        speed_weight*rub_max_bool[0] ,
                        speed_weight*rub_max_bool[1] , 
                        speed_weight*rub_max_bool[1]  ; 
    }
  
  
  d->Lu +=(step_weights_.array()*d->r.template tail<4>().array()).matrix() ; 
  
  // Hessian : Lxx
  d->Lxx.diagonal().head(12) = (state_weights_.array() * state_weights_.array()).matrix() ;  
  d->Lxx.diagonal().segment(12,8) = (heuristicWeights.array() * heuristicWeights.array()).matrix() ;  

  
  
  d->Luu.diagonal() += (step_weights_.array() * step_weights_.array()).matrix() ;

  // Dynamic derivatives
  d->Fx.setIdentity();
  d->Fu.block(12,0,8,4) = B;  
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
  if (static_cast<std::size_t>(weights.size()) != 12 ) {
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
  if (static_cast<std::size_t>(weights.size()) != 4 ) {
    throw_pretty("Invalid argument: "
                 << "Weights vector has wrong dimension (it should be 4)");
  }
  step_weights_ = weights;
}

template <typename Scalar>
const typename Eigen::Matrix<Scalar, 8, 1>& ActionModelQuadrupedStepTimeTpl<Scalar>::get_heuristic_weights() const {
  return heuristicWeights;
}
template <typename Scalar>
void ActionModelQuadrupedStepTimeTpl<Scalar>::set_heuristic_weights(const typename MathBase::VectorXs& weights) {
  if (static_cast<std::size_t>(weights.size()) != 8 ) {
    throw_pretty("Invalid argument: "
                 << "Weights vector has wrong dimension (it should be 8)");
  }
  heuristicWeights = weights;
}

template <typename Scalar>
const typename Eigen::Matrix<Scalar, 8, 1>& ActionModelQuadrupedStepTimeTpl<Scalar>::get_shoulder_position() const {
  return pshoulder_ ;
}
template <typename Scalar>
void ActionModelQuadrupedStepTimeTpl<Scalar>::set_shoulder_position(const typename MathBase::VectorXs& pos) {
  if (static_cast<std::size_t>(pos.size()) != 8 ) {
    throw_pretty("Invalid argument: "
                 << "Weights vector has wrong dimension (it should be 8)");
  }
  pshoulder_ = pos;
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
  beta_lim = Scalar((64*nb_nodes*nb_nodes*vlim*vlim)/225) ; ;
}

template <typename Scalar>
const Scalar& ActionModelQuadrupedStepTimeTpl<Scalar>::get_vlim() const {
  return vlim;
}
template <typename Scalar>
void ActionModelQuadrupedStepTimeTpl<Scalar>::set_vlim(const Scalar& vlim_) {
  vlim = vlim_; 
  beta_lim = Scalar((64*nb_nodes*nb_nodes*vlim*vlim)/225) ; ;
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
void ActionModelQuadrupedStepTimeTpl<Scalar>::update_model(const Eigen::Ref<const typename MathBase::MatrixXs>& l_feet  ,
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

  // In this model, l_feet contain : 
  // pheur_x  -  p0_x1   ;  pheur_y  -  p0_y1    ;  p0_x2   ;  p0_y2
  // v0_x1   ;  v0_y1   ;  v0_x2   ;  v0_y2
  // acc0_x1 ;  acc0_y1 ;  acc0_x2 ;  acc0_y2 ;
  lfeet = l_feet ;

  for (int i=0; i<4; i=i+1){
    b_coeff2.col(3*i) = nb_nodes*lfeet(1,i)*b_coeff.col(0) ; 
    b_coeff2.col(3*i+1) = nb_nodes*nb_nodes*lfeet(2,i)*b_coeff.col(1) ; 
    b_coeff2.col(3*i+2) = b_coeff.col(2) ;
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
