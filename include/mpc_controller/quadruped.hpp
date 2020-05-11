#ifndef __mpc_controller_handtracking_hpp__
#define __mpc_controller_handtracking_hpp__

#include <stdexcept>

#include "crocoddyl/core/fwd.hpp"
#include "crocoddyl/core/action-base.hpp"
#include "crocoddyl/core/states/euclidean.hpp"

namespace mpc_controller {
template <typename _Scalar>
class ActionModelQuadrupedTpl : public crocoddyl::ActionModelAbstractTpl<_Scalar> {
 public:
  typedef _Scalar Scalar;
  typedef crocoddyl::ActionDataAbstractTpl<Scalar> ActionDataAbstract;
  typedef crocoddyl::ActionModelAbstractTpl<Scalar> Base;
  typedef crocoddyl::MathBaseTpl<Scalar> MathBase;

  ActionModelQuadrupedTpl();
  ~ActionModelQuadrupedTpl();

  virtual void calc(const boost::shared_ptr<ActionDataAbstract>& data,
                    const Eigen::Ref<const typename MathBase::VectorXs>& x,
                    const Eigen::Ref<const typename MathBase::VectorXs>& u);
  virtual void calcDiff(const boost::shared_ptr<ActionDataAbstract>& data,
                        const Eigen::Ref<const typename MathBase::VectorXs>& x,
                        const Eigen::Ref<const typename MathBase::VectorXs>& u);
  virtual boost::shared_ptr<ActionDataAbstract> createData();

  const typename MathBase::Vector2s& get_cost_weights() const;
  void set_cost_weights(const typename MathBase::Vector2s& weights);

 protected:
  using Base::has_control_limits_;  //!< Indicates whether any of the control limits
  using Base::nr_;                  //!< Dimension of the cost residual
  using Base::nu_;                  //!< Control dimension
  using Base::state_;               //!< Model of the state
  using Base::u_lb_;                //!< Lower control limits
  using Base::u_ub_;                //!< Upper control limits
  using Base::unone_;               //!< Neutral state

 private:
  typename MathBase::Vector2s cost_weights_;
  Scalar dt_;
};

template <typename _Scalar>
struct ActionDataQuadrupedTpl : public crocoddyl::ActionDataAbstractTpl<_Scalar> {
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
  explicit ActionDataQuadrupedTpl(Model<Scalar>* const model) : crocoddyl::ActionDataAbstractTpl<Scalar>(model) {}
};


typedef ActionModelQuadrupedTpl<double> ActionModelQuadruped;
typedef ActionDataQuadrupedTpl<double> ActionDataQuadruped;

}


#endif
