#ifndef __quadruped_walkgen_quadruped_step_hpp__
#define __quadruped_walkgen_quadruped_step_hpp__
#include <stdexcept>

#include "crocoddyl/core/fwd.hpp"
#include "crocoddyl/core/action-base.hpp"
#include "crocoddyl/core/states/euclidean.hpp"
#include "crocoddyl/multibody/friction-cone.hpp"

#include "crocoddyl/core/utils/timer.hpp"

namespace quadruped_walkgen {
template <typename _Scalar>
class ActionModelQuadrupedStepTpl : public crocoddyl::ActionModelAbstractTpl<_Scalar> {
 public:
  typedef _Scalar Scalar;
  typedef crocoddyl::ActionDataAbstractTpl<Scalar> ActionDataAbstract;
  typedef crocoddyl::ActionModelAbstractTpl<Scalar> Base;
  typedef crocoddyl::MathBaseTpl<Scalar> MathBase;
  
  ActionModelQuadrupedStepTpl();
  ~ActionModelQuadrupedStepTpl();

  virtual void calc(const boost::shared_ptr<ActionDataAbstract>& data,
                    const Eigen::Ref<const typename MathBase::VectorXs>& x,
                    const Eigen::Ref<const typename MathBase::VectorXs>& u);
  virtual void calcDiff(const boost::shared_ptr<ActionDataAbstract>& data,
                        const Eigen::Ref<const typename MathBase::VectorXs>& x,
                        const Eigen::Ref<const typename MathBase::VectorXs>& u);
  virtual boost::shared_ptr<ActionDataAbstract> createData();

  // Get and Set weights vectors : state , force & friction cone :
  const typename Eigen::Matrix<Scalar, 12, 1>& get_state_weights() const;
  void set_state_weights(const typename MathBase::VectorXs& weights);

  const typename Eigen::Matrix<Scalar, 4, 1>& get_step_weights() const;
  void set_step_weights(const typename MathBase::VectorXs& weights);

  const typename Eigen::Matrix<Scalar, 8, 1>& get_shoulder_weights() const;
  void set_shoulder_weights(const typename MathBase::VectorXs& weights);

  const typename Eigen::Matrix<Scalar, 8, 1>& get_shoulder_position() const;
  void set_shoulder_position(const typename MathBase::VectorXs& weights);


  // Update the model depending if the foot in contact with the ground 
  // or the new lever arms
  void update_model(const Eigen::Ref<const typename MathBase::MatrixXs>& l_feet  ,
                    const Eigen::Ref<const typename MathBase::MatrixXs>& xref,
                    const Eigen::Ref<const typename MathBase::MatrixXs>& S ) ;


 protected:
  using Base::has_control_limits_;  //!< Indicates whether any of the control limits
  using Base::nr_;                  //!< Dimension of the cost residual
  using Base::nu_;                  //!< Control dimension
  using Base::state_;               //!< Model of the state
  using Base::u_lb_;                //!< Lower control limits
  using Base::u_ub_;                //!< Upper control limits
  using Base::unone_;               //!< Neutral state
 
 
 private:

  typename Eigen::Matrix<Scalar, 12, 1> state_weights_;
  typename Eigen::Matrix<Scalar, 4, 1> step_weights_;  
  typename Eigen::Matrix<Scalar, 8, 1> shoulder_weights_;

  
  
  typename Eigen::Matrix<Scalar, 8, 4 > B;

  typename MathBase::MatrixXs xref_;  

  typename Eigen::Matrix<Scalar, 8 , 1 > pshoulder_;  
};

template <typename _Scalar>
struct ActionDataQuadrupedStepTpl : public crocoddyl::ActionDataAbstractTpl<_Scalar> {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  typedef _Scalar Scalar;
  typedef crocoddyl::MathBaseTpl<Scalar> MathBase;
  typedef crocoddyl::ActionDataAbstractTpl<Scalar> Base;
  using Base::cost;
  using Base::Fu;
  using Base::Fx;
  using Base::Lu;
  using Base::Luu;
  using Base::Lx;
  using Base::Lxu;
  using Base::Lxx;
  using Base::r;
  using Base::xnext;

  template <template <typename Scalar> class Model>
  explicit ActionDataQuadrupedStepTpl(Model<Scalar>* const model) : crocoddyl::ActionDataAbstractTpl<Scalar>(model) {}
};

  
 

/* --- Details -------------------------------------------------------------- */
/* --- Details -------------------------------------------------------------- */
/* --- Details -------------------------------------------------------------- */

typedef ActionModelQuadrupedStepTpl<double> ActionModelQuadrupedStep;
typedef ActionDataQuadrupedStepTpl<double> ActionDataQuadrupedStep;

}

#include "quadruped_step.hxx"


#endif
