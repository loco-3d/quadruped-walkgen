#ifndef __quadruped_walkgen_quadruped_step_period_hpp__
#define __quadruped_walkgen_quadruped_step_period_hpp__
#include <stdexcept>

#include "crocoddyl/core/fwd.hpp"
#include "crocoddyl/core/action-base.hpp"
#include "crocoddyl/core/states/euclidean.hpp"
#include "crocoddyl/multibody/friction-cone.hpp"

#include "crocoddyl/core/utils/timer.hpp"

namespace quadruped_walkgen {
template <typename _Scalar>
class ActionModelQuadrupedStepPeriodTpl : public crocoddyl::ActionModelAbstractTpl<_Scalar> {
 public:
  typedef _Scalar Scalar;
  typedef crocoddyl::ActionDataAbstractTpl<Scalar> ActionDataAbstract;
  typedef crocoddyl::ActionModelAbstractTpl<Scalar> Base;
  typedef crocoddyl::MathBaseTpl<Scalar> MathBase;

  ActionModelQuadrupedStepPeriodTpl();
  ~ActionModelQuadrupedStepPeriodTpl();

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
  void update_model(const Eigen::Ref<const typename MathBase::MatrixXs>& l_feet,
                    const Eigen::Ref<const typename MathBase::MatrixXs>& xref,
                    const Eigen::Ref<const typename MathBase::MatrixXs>& S);

  const bool& get_symmetry_term() const;
  void set_symmetry_term(const bool& sym_term);

  const bool& get_centrifugal_term() const;
  void set_centrifugal_term(const bool& cent_term);

  const Scalar& get_T_gait() const;
  void set_T_gait(const Scalar& T_gait_);

  const Scalar& get_dt_ref() const;
  void set_dt_ref(const Scalar& dt);

  const Scalar& get_dt_min() const;
  void set_dt_min(const Scalar& dt);

  const Scalar& get_dt_max() const;
  void set_dt_max(const Scalar& dt);

  const Scalar& get_dt_weight() const;
  void set_dt_weight(const Scalar& weight_);

  const Scalar& get_dt_bound_weight() const;
  void set_dt_bound_weight(const Scalar& weight_);

  const Scalar& get_nb_nodes() const;
  void set_nb_nodes(const Scalar& nodes_);

  const Scalar& get_vlim() const;
  void set_vlim(const Scalar& vlim_);

  const Scalar& get_speed_weight() const;
  void set_speed_weight(const Scalar& weight_);

 protected:
  using Base::has_control_limits_;  //!< Indicates whether any of the control limits
  using Base::nr_;                  //!< Dimension of the cost residual
  using Base::nu_;                  //!< Control dimension
  using Base::state_;               //!< Model of the state
  using Base::u_lb_;                //!< Lower control limits
  using Base::u_ub_;                //!< Upper control limits
  using Base::unone_;               //!< Neutral state

 private:
  Scalar T_gait;
  Scalar dt_weight_;
  Scalar dt_bound_weight;
  Scalar speed_weight;
  bool centrifugal_term;
  bool symmetry_term;

  Scalar nb_nodes;
  Scalar vlim;
  Scalar beta_lim;

  typename Eigen::Matrix<Scalar, 12, 1> state_weights_;
  typename Eigen::Matrix<Scalar, 4, 1> step_weights_;
  typename Eigen::Matrix<Scalar, 8, 1> shoulder_weights_;
  typename MathBase::Matrix3s R_tmp;

  typename Eigen::Matrix<Scalar, 8, 4> B;

  typename MathBase::MatrixXs xref_;

  typename Eigen::Matrix<Scalar, 8, 1> pshoulder_;
  typename Eigen::Matrix<Scalar, 2, 4> pshoulder_0;
  typename Eigen::Matrix<Scalar, 2, 4> pshoulder_tmp;

  typename Eigen::Matrix<Scalar, 3, 1> pcentrifugal_tmp;
  typename Eigen::Matrix<Scalar, 3, 1> pcentrifugal_tmp_1;
  typename Eigen::Matrix<Scalar, 3, 1> pcentrifugal_tmp_2;

  typename Eigen::Matrix<Scalar, 1, 1> dt_ref_;
  typename Eigen::Matrix<Scalar, 1, 1> dt_min_;
  typename Eigen::Matrix<Scalar, 1, 1> dt_max_;

  typename Eigen::Matrix<Scalar, 4, 1> rub_max_;
  typename Eigen::Matrix<Scalar, 4, 1> rub_max_bool;
};

template <typename _Scalar>
struct ActionDataQuadrupedStepPeriodTpl : public crocoddyl::ActionDataAbstractTpl<_Scalar> {
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
  explicit ActionDataQuadrupedStepPeriodTpl(Model<Scalar>* const model)
      : crocoddyl::ActionDataAbstractTpl<Scalar>(model) {}
};

/* --- Details -------------------------------------------------------------- */
/* --- Details -------------------------------------------------------------- */
/* --- Details -------------------------------------------------------------- */

typedef ActionModelQuadrupedStepPeriodTpl<double> ActionModelQuadrupedStepPeriod;
typedef ActionDataQuadrupedStepPeriodTpl<double> ActionDataQuadrupedStepPeriod;

}  // namespace quadruped_walkgen

#include "quadruped_step_period.hxx"

#endif
