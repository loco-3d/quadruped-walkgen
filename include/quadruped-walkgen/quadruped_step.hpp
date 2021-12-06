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

  const typename Eigen::Matrix<Scalar, 8, 1>& get_step_weights() const;
  void set_step_weights(const typename MathBase::VectorXs& weights);

  const typename Eigen::Matrix<Scalar, 8, 1>& get_heuristic_weights() const;
  void set_heuristic_weights(const typename MathBase::VectorXs& weights);

  // Update the model depending if the foot in contact with the ground
  // or the new lever arms
  void update_model(const Eigen::Ref<const typename MathBase::MatrixXs>& l_feet,
                    const Eigen::Ref<const typename MathBase::MatrixXs>& xref,
                    const Eigen::Ref<const typename MathBase::VectorXs>& S,
                    const Eigen::Ref<const typename MathBase::MatrixXs>& position,
                    const Eigen::Ref<const typename MathBase::MatrixXs>& velocity,
                    const Eigen::Ref<const typename MathBase::MatrixXs>& acceleration,
                    const Eigen::Ref<const typename MathBase::MatrixXs>& jerk,
                    const Eigen::Ref<const typename MathBase::MatrixXs>& oRh,
                    const Eigen::Ref<const typename MathBase::MatrixXs>& oTh, const Scalar& delta_T);

  const bool& get_symmetry_term() const;
  void set_symmetry_term(const bool& sym_term);

  const bool& get_centrifugal_term() const;
  void set_centrifugal_term(const bool& cent_term);

  const Scalar& get_T_gait() const;
  void set_T_gait(const Scalar& T_gait_);

  const bool& get_acc_activated() const;
  void set_acc_activated(const bool& is_activated);

  const typename Eigen::Matrix<Scalar, 2, 1>& get_acc_lim() const;
  void set_acc_lim(const typename MathBase::VectorXs& acceleration_lim_);

  const Scalar& get_acc_weight() const;
  void set_acc_weight(const Scalar& weight_);

  const bool& get_vel_activated() const;
  void set_vel_activated(const bool& is_activated);

  const typename Eigen::Matrix<Scalar, 2, 1>& get_vel_lim() const;
  void set_vel_lim(const typename MathBase::VectorXs& velocity_lim_);

  const Scalar& get_vel_weight() const;
  void set_vel_weight(const Scalar& weight_);

  void set_sample_feet_traj(const int& n_sample);

  const bool& get_jerk_activated() const;
  void set_jerk_activated(const bool& is_activated);

  const Scalar& get_jerk_weight() const;
  void set_jerk_weight(const Scalar& weight_);

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
  bool centrifugal_term;
  bool symmetry_term;

  typename Eigen::Matrix<Scalar, 12, 1> state_weights_;
  typename Eigen::Matrix<Scalar, 8, 1> step_weights_;
  typename Eigen::Matrix<Scalar, 8, 1> heuristic_weights_;
  typename MathBase::Matrix3s R_tmp;

  typename Eigen::Matrix<Scalar, 8, 8> B;

  typename MathBase::MatrixXs xref_;
  typename Eigen::Matrix<Scalar, 8, 1> pheuristic_;

  // typename Eigen::Matrix<Scalar, 2, 4> pshoulder_0;
  // typename Eigen::Matrix<Scalar, 2, 4> pshoulder_tmp;
  // typename Eigen::Matrix<Scalar, 3, 1> pcentrifugal_tmp;
  // typename Eigen::Matrix<Scalar, 3, 1> pcentrifugal_tmp_1;
  // typename Eigen::Matrix<Scalar, 3, 1> pcentrifugal_tmp_2;

  // Cost on the acceleration of the feet :
  int N_sampling;
  bool is_acc_activated_;                         // Boolean to activate the cost on the acceleration of the feet
  Scalar acc_weight_;                             // Weight on the acceleration cost
  typename Eigen::Matrix<Scalar, 2, 1> acc_lim_;  // Maximum acceleration allowed on x and y axis
  typename Eigen::Matrix<Scalar, 4, 1> S_;        // Containing the moving feet
  typename Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> delta_;
  typename Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> gamma_;
  typename Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> alpha_;
  typename Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> beta_x_;
  typename Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> beta_y_;
  typename Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> tmp_ones_;
  typename Eigen::Array<Scalar, 3, 4> position_;

  typename Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> rb_accx_max_;
  typename Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> rb_accy_max_;
  typename Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> rb_accx_max_bool_;
  typename Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> rb_accy_max_bool_;

  // Cost on the velocity of the feet :
  bool is_vel_activated_;                         // Boolean to activate the cost on the velocity of the feet
  Scalar vel_weight_;                             // Weight on the velocity cost
  typename Eigen::Matrix<Scalar, 2, 1> vel_lim_;  // Maximum velocity allowed on x and y axis
  typename Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> gamma_v;
  typename Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> alpha_v;
  typename Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> beta_x_v;
  typename Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> beta_y_v;

  typename Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> rb_velx_max_;
  typename Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> rb_vely_max_;
  typename Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> rb_velx_max_bool_;
  typename Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> rb_vely_max_bool_;

  // Cost on the jerk of the feet
  bool is_jerk_activated_;
  Scalar jerk_weight_;
  Scalar alpha_j;
  typename Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> beta_j;
  typename Eigen::Array<Scalar, 3, 4> jerk_;

  typename Eigen::Matrix<Scalar, 2, 4> rb_jerk_;

  typename Eigen::Matrix<Scalar, 3, 3> oRh_;
  typename Eigen::Matrix<Scalar, 3, 1> oTh_;
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

}  // namespace quadruped_walkgen

#include "quadruped_step.hxx"

#endif
