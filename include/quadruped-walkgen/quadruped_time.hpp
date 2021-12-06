#ifndef __quadruped_walkgen_quadruped_time_hpp__
#define __quadruped_walkgen_quadruped_time_hpp__
#include <stdexcept>

#include "crocoddyl/core/fwd.hpp"
#include "crocoddyl/core/action-base.hpp"
#include "crocoddyl/core/states/euclidean.hpp"
#include "crocoddyl/multibody/friction-cone.hpp"

#include "crocoddyl/core/utils/timer.hpp"

namespace quadruped_walkgen {
template <typename _Scalar>
class ActionModelQuadrupedTimeTpl : public crocoddyl::ActionModelAbstractTpl<_Scalar> {
 public:
  typedef _Scalar Scalar;
  typedef crocoddyl::ActionDataAbstractTpl<Scalar> ActionDataAbstract;
  typedef crocoddyl::ActionModelAbstractTpl<Scalar> Base;
  typedef crocoddyl::MathBaseTpl<Scalar> MathBase;

  ActionModelQuadrupedTimeTpl();
  ~ActionModelQuadrupedTimeTpl();

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

  const typename Eigen::Matrix<Scalar, 8, 1>& get_heuristic_weights() const;
  void set_heuristic_weights(const typename MathBase::VectorXs& weights);

  // Update the model depending if the foot in contact with the ground
  // or the new lever arms
  void update_model(const Eigen::Ref<const typename MathBase::MatrixXs>& l_feet,
                    const Eigen::Ref<const typename MathBase::MatrixXs>& xref,
                    const Eigen::Ref<const typename MathBase::VectorXs>& S);

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

  // Command weights
  const Scalar& get_dt_weight_cmd() const;
  void set_dt_weight_cmd(const Scalar& weight_);

  const Scalar& get_dt_bound_weight_cmd() const;
  void set_dt_bound_weight_cmd(const Scalar& weight_);

  // get cost
  const typename Eigen::Matrix<Scalar, 7, 1>& get_cost() const;

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
  Scalar dt_weight_cmd;
  Scalar dt_bound_weight_cmd;
  bool centrifugal_term;
  bool symmetry_term;

  // Log cost
  bool log_cost;
  typename Eigen::Matrix<Scalar, 7, 1> cost_;

  typename Eigen::Matrix<Scalar, 12, 1> state_weights_;
  typename Eigen::Matrix<Scalar, 8, 1> heuristic_weights_;

  typename MathBase::Matrix3s R_tmp;

  typename MathBase::MatrixXs xref_;
  typename Eigen::Matrix<Scalar, 8, 1> pheuristic_;
  typename Eigen::Matrix<Scalar, 8, 1> gait_double_;
  // typename Eigen::Matrix<Scalar, 2 , 4 > pshoulder_0;
  // typename Eigen::Matrix<Scalar, 2 , 4 > pshoulder_tmp;

  // typename Eigen::Matrix<Scalar, 3 , 1 > pcentrifugal_tmp;
  // typename Eigen::Matrix<Scalar, 3 , 1 > pcentrifugal_tmp_1;
  // typename Eigen::Matrix<Scalar, 3 , 1 > pcentrifugal_tmp_2;

  typename Eigen::Matrix<Scalar, 1, 1> dt_ref_;
  typename Eigen::Matrix<Scalar, 1, 1> dt_min_;
  typename Eigen::Matrix<Scalar, 1, 1> dt_max_;

  typename Eigen::Matrix<Scalar, 2, 1> rub_max_;
  typename Eigen::Matrix<Scalar, 2, 1> rub_max_bool;
};

template <typename _Scalar>
struct ActionDataQuadrupedTimeTpl : public crocoddyl::ActionDataAbstractTpl<_Scalar> {
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
  explicit ActionDataQuadrupedTimeTpl(Model<Scalar>* const model) : crocoddyl::ActionDataAbstractTpl<Scalar>(model) {}
};

/* --- Details -------------------------------------------------------------- */
/* --- Details -------------------------------------------------------------- */
/* --- Details -------------------------------------------------------------- */

typedef ActionModelQuadrupedTimeTpl<double> ActionModelQuadrupedTime;
typedef ActionDataQuadrupedTimeTpl<double> ActionDataQuadrupedTime;

}  // namespace quadruped_walkgen

#include "quadruped_time.hxx"

#endif
