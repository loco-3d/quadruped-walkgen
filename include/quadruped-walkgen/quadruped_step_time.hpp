#ifndef __quadruped_walkgen_quadruped_step_time_hpp__
#define __quadruped_walkgen_quadruped_step_time_hpp__
#include <stdexcept>

#include "crocoddyl/core/fwd.hpp"
#include "crocoddyl/core/action-base.hpp"
#include "crocoddyl/core/states/euclidean.hpp"
#include "crocoddyl/multibody/friction-cone.hpp"

#include "crocoddyl/core/utils/timer.hpp"

namespace quadruped_walkgen {
template <typename _Scalar>
class ActionModelQuadrupedStepTimeTpl : public crocoddyl::ActionModelAbstractTpl<_Scalar> {
 public:
  typedef _Scalar Scalar;
  typedef crocoddyl::ActionDataAbstractTpl<Scalar> ActionDataAbstract;
  typedef crocoddyl::ActionModelAbstractTpl<Scalar> Base;
  typedef crocoddyl::MathBaseTpl<Scalar> MathBase;

  ActionModelQuadrupedStepTimeTpl();
  ~ActionModelQuadrupedStepTimeTpl();

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

  const typename Eigen::Matrix<Scalar, 8, 1>& get_heuristic_weights() const;
  void set_heuristic_weights(const typename MathBase::VectorXs& weights);

  // Update the model depending if the foot in contact with the ground
  // or the new lever arms
  void update_model(const Eigen::Ref<const typename MathBase::MatrixXs>& l_feet,
                    const Eigen::Ref<const typename MathBase::MatrixXs>& velocity,
                    const Eigen::Ref<const typename MathBase::MatrixXs>& acceleration,
                    const Eigen::Ref<const typename MathBase::MatrixXs>& xref,
                    const Eigen::Ref<const typename MathBase::VectorXs>& S);

  const bool& get_symmetry_term() const;
  void set_symmetry_term(const bool& sym_term);

  const bool& get_centrifugal_term() const;
  void set_centrifugal_term(const bool& cent_term);

  const Scalar& get_T_gait() const;
  void set_T_gait(const Scalar& T_gait_);

  const Scalar& get_nb_nodes() const;
  void set_nb_nodes(const Scalar& nodes_);

  const Scalar& get_vlim() const;
  void set_vlim(const Scalar& vlim_);

  const Scalar& get_speed_weight() const;
  void set_speed_weight(const Scalar& weight_);

  const bool& get_first_step() const;
  void set_first_step(const bool& first);

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
  Scalar speed_weight;
  Scalar nb_nodes;
  Scalar vlim;
  Scalar beta_lim;
  int nb_alpha_;
  bool centrifugal_term;
  bool symmetry_term;
  // indicates whether it t the 1st step, otherwise the cost function is much simpler (acc, speed = 0)
  bool first_step;

  typename Eigen::Matrix<Scalar, 12, 1> state_weights_;
  typename Eigen::Matrix<Scalar, 4, 1> step_weights_;
  typename Eigen::Matrix<Scalar, 8, 1> heuristicWeights;
  typename MathBase::Matrix3s R_tmp;

  // typename  Eigen::Array<Scalar, 3, 1 > alpha ;
  // typename  Eigen::Array<Scalar, 3, 4 > alpha2 ;
  // typename  Eigen::Array<Scalar, 3, 3 > b_coeff ;
  typename MathBase::ArrayXs alpha;
  typename Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> alpha2;
  typename Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> b_coeff;
  typename Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> b_coeff_x0;
  typename Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> b_coeff_y0;
  typename Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> b_coeff_x1;
  typename Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> b_coeff_y1;
  typename Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> b_coeff_x2;
  typename Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> b_coeff_y2;

  typename Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> rub_max_first_x;
  typename Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> rub_max_first_y;
  typename Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> rub_max_first_2;
  typename Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> rub_max_first_bool;

  // typename  Eigen::Array<Scalar, 3, 12 > b_coeff2 ;
  typename Eigen::Matrix<Scalar, 3, 4> lfeet;
  // typename Eigen::Array<Scalar, 3, 4 > rub_max_first ;
  // typename Eigen::Array<Scalar, 3, 2 > rub_max_first_2 ;
  // typename Eigen::Array<Scalar, 3, 2 > rub_max_bool_first ;

  typename Eigen::Matrix<Scalar, 8, 8> B;

  typename MathBase::MatrixXs xref_;
  typename MathBase::VectorXs S_;  // Containing the flying feet
  typename Eigen::Matrix<Scalar, 8, 1> pheuristic_;

  // Compute heuristic inside update Model
  // typename Eigen::Matrix<Scalar, 2 , 4 > pshoulder_0;
  // typename Eigen::Matrix<Scalar, 2 , 4 > pshoulder_tmp;
  // typename Eigen::Matrix<Scalar, 3 , 1 > pcentrifugal_tmp;
  // typename Eigen::Matrix<Scalar, 3 , 1 > pcentrifugal_tmp_1;
  // typename Eigen::Matrix<Scalar, 3 , 1 > pcentrifugal_tmp_2;

  typename Eigen::Matrix<Scalar, 4, 1> rub_max_;
  typename Eigen::Matrix<Scalar, 4, 1> rub_max_bool;

  // Log cost
  bool log_cost;
  typename Eigen::Matrix<Scalar, 7, 1> cost_;
};

template <typename _Scalar>
struct ActionDataQuadrupedStepTimeTpl : public crocoddyl::ActionDataAbstractTpl<_Scalar> {
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
  explicit ActionDataQuadrupedStepTimeTpl(Model<Scalar>* const model)
      : crocoddyl::ActionDataAbstractTpl<Scalar>(model) {}
};

/* --- Details -------------------------------------------------------------- */
/* --- Details -------------------------------------------------------------- */
/* --- Details -------------------------------------------------------------- */

typedef ActionModelQuadrupedStepTimeTpl<double> ActionModelQuadrupedStepTime;
typedef ActionDataQuadrupedStepTimeTpl<double> ActionDataQuadrupedStepTime;

}  // namespace quadruped_walkgen

#include "quadruped_step_time.hxx"

#endif
