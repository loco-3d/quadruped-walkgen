#ifndef __quadruped_walkgen_quadruped_augmented_time_hpp__
#define __quadruped_walkgen_quadruped_augmented_time_hpp__
#include <stdexcept>

#include "crocoddyl/core/fwd.hpp"
#include "crocoddyl/core/action-base.hpp"
#include "crocoddyl/core/states/euclidean.hpp"
#include "crocoddyl/multibody/friction-cone.hpp"

#include "crocoddyl/core/utils/timer.hpp"

namespace quadruped_walkgen {
template <typename _Scalar>
class ActionModelQuadrupedAugmentedTimeTpl : public crocoddyl::ActionModelAbstractTpl<_Scalar> {
 public:
  typedef _Scalar Scalar;
  typedef crocoddyl::ActionDataAbstractTpl<Scalar> ActionDataAbstract;
  typedef crocoddyl::ActionModelAbstractTpl<Scalar> Base;
  typedef crocoddyl::MathBaseTpl<Scalar> MathBase;

  ActionModelQuadrupedAugmentedTimeTpl();
  ~ActionModelQuadrupedAugmentedTimeTpl();

  virtual void calc(const boost::shared_ptr<ActionDataAbstract>& data,
                    const Eigen::Ref<const typename MathBase::VectorXs>& x,
                    const Eigen::Ref<const typename MathBase::VectorXs>& u);
  virtual void calcDiff(const boost::shared_ptr<ActionDataAbstract>& data,
                        const Eigen::Ref<const typename MathBase::VectorXs>& x,
                        const Eigen::Ref<const typename MathBase::VectorXs>& u);
  virtual boost::shared_ptr<ActionDataAbstract> createData();

  // Get and Set weights vectors : state , force & friction cone :
  const typename Eigen::Matrix<Scalar, 12, 1>& get_force_weights() const;
  void set_force_weights(const typename MathBase::VectorXs& weights);

  const typename Eigen::Matrix<Scalar, 12, 1>& get_state_weights() const;
  void set_state_weights(const typename MathBase::VectorXs& weights);

  const typename Eigen::Matrix<Scalar, 8, 1>& get_heuristic_weights() const;
  void set_heuristic_weights(const typename MathBase::VectorXs& weights);

  const typename Eigen::Matrix<Scalar, 8, 1>& get_stop_weights() const;
  void set_stop_weights(const typename MathBase::VectorXs& weights);

  const Scalar& get_friction_weight() const;
  void set_friction_weight(const Scalar& weight);

  const Scalar& get_mu() const;
  void set_mu(const Scalar& mu_coeff);

  const Scalar& get_mass() const;
  void set_mass(const Scalar& m);

  const Scalar& get_dt_ref() const;
  void set_dt_ref(const Scalar& dt);

  const Scalar& get_dt_min() const;
  void set_dt_min(const Scalar& dt);

  const Scalar& get_dt_max() const;
  void set_dt_max(const Scalar& dt);

  const typename Eigen::Matrix<Scalar, 3, 3>& get_gI() const;
  void set_gI(const typename MathBase::Matrix3s& inertia_matrix);

  const Scalar& get_min_fz_contact() const;
  void set_min_fz_contact(const Scalar& min_fz);

  const Scalar& get_max_fz_contact() const;
  void set_max_fz_contact(const Scalar& max_fz_);

  const typename Eigen::Matrix<Scalar, 8, 1>& get_shoulder_position() const;
  void set_shoulder_position(const typename MathBase::VectorXs& weights);

  const bool& get_symmetry_term() const;
  void set_symmetry_term(const bool& sym_term);

  const bool& get_centrifugal_term() const;
  void set_centrifugal_term(const bool& cent_term);

  const Scalar& get_T_gait() const;
  void set_T_gait(const Scalar& T_gait_);

  const Scalar& get_dt_weight() const;
  void set_dt_weight(const Scalar& weight_);

  const Scalar& get_dt_bound_weight() const;
  void set_dt_bound_weight(const Scalar& weight_);

  const bool& get_relative_forces() const;
  void set_relative_forces(const bool& rel_forces);

  // Set parameter relative to the shoulder height cost
  const Scalar& get_shoulder_hlim() const;
  void set_shoulder_hlim(const Scalar& hlim);

  const Scalar& get_shoulder_contact_weight() const;
  void set_shoulder_contact_weight(const Scalar& weight);

  // Update the model depending if the foot in contact with the ground
  // or the new lever arms
  void update_model(const Eigen::Ref<const typename MathBase::MatrixXs>& l_feet,
                    const Eigen::Ref<const typename MathBase::MatrixXs>& l_stop,
                    const Eigen::Ref<const typename MathBase::MatrixXs>& xref,
                    const Eigen::Ref<const typename MathBase::MatrixXs>& S);

  // Get A & B matrix
  const typename Eigen::Matrix<Scalar, 12, 12>& get_A() const;
  const typename Eigen::Matrix<Scalar, 12, 12>& get_B() const;

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
  Scalar dt_weight_;
  Scalar mass;
  Scalar mu;
  Scalar friction_weight_;
  Scalar min_fz_in_contact;
  Scalar max_fz;
  Scalar T_gait;
  Scalar dt_bound_weight;
  bool centrifugal_term;
  bool symmetry_term;

  bool relative_forces;
  typename Eigen::Matrix<Scalar, 12, 1> uref_;

  typename Eigen::Matrix<Scalar, 12, 1> force_weights_;
  typename Eigen::Matrix<Scalar, 12, 1> state_weights_;
  typename Eigen::Matrix<Scalar, 8, 1> heuristicWeights;
  typename Eigen::Matrix<Scalar, 8, 1> last_position_weights_;

  typename Eigen::Matrix<Scalar, 12, 12> A;
  typename Eigen::Matrix<Scalar, 12, 12> B;
  typename Eigen::Matrix<Scalar, 12, 1> g;
  typename Eigen::Matrix<Scalar, 3, 3> R;
  typename MathBase::Matrix3s R_tmp;
  typename Eigen::Matrix<Scalar, 3, 3> gI;

  typename Eigen::Matrix<Scalar, 3, 4> lever_arms;
  typename MathBase::Vector3s lever_tmp;
  typename MathBase::MatrixXs xref_;

  typename Eigen::Matrix<Scalar, 8, 1> pshoulder_;
  typename Eigen::Matrix<Scalar, 8, 1> pheuristic_;
  typename Eigen::Matrix<Scalar, 2, 4> pshoulder_0;
  typename Eigen::Matrix<Scalar, 2, 4> pshoulder_tmp;

  typename Eigen::Matrix<Scalar, 3, 1> pcentrifugal_tmp;
  typename Eigen::Matrix<Scalar, 3, 1> pcentrifugal_tmp_1;
  typename Eigen::Matrix<Scalar, 3, 1> pcentrifugal_tmp_2;

  typename Eigen::Matrix<Scalar, 8, 1> pref_;

  typename Eigen::Matrix<Scalar, 24, 1> ub;

  typename Eigen::Matrix<Scalar, 24, 1> Fa_x_u;
  typename Eigen::Matrix<Scalar, 24, 1> rub_max_;
  typename Eigen::Matrix<Scalar, 24, 24> Arr;
  typename Eigen::Matrix<Scalar, 6, 1> r;
  typename Eigen::Matrix<Scalar, 4, 1> gait;
  typename Eigen::Matrix<Scalar, 8, 1> gait_double;

  typename Eigen::Matrix<Scalar, 3, 1> base_vector_x;
  typename Eigen::Matrix<Scalar, 3, 1> base_vector_y;
  typename Eigen::Matrix<Scalar, 3, 1> base_vector_z;
  typename Eigen::Matrix<Scalar, 3, 1> forces_3d;

  typename Eigen::Matrix<Scalar, 1, 1> dt_min_;
  typename Eigen::Matrix<Scalar, 1, 1> dt_max_;

  typename Eigen::Matrix<Scalar, 2, 1> rub_max_dt;
  typename Eigen::Matrix<Scalar, 2, 1> rub_max_dt_bool;

  // Log cost
  bool log_cost;
  typename Eigen::Matrix<Scalar, 7, 1> cost_;

  // Cost relative to the shoulder height
  typename Eigen::Matrix<Scalar, 3, 4> psh;
  typename Eigen::Matrix<Scalar, 4, 1> sh_ub_max_;
  Scalar sh_weight;
  Scalar sh_hlim;
};

template <typename _Scalar>
struct ActionDataQuadrupedAugmentedTimeTpl : public crocoddyl::ActionDataAbstractTpl<_Scalar> {
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
  explicit ActionDataQuadrupedAugmentedTimeTpl(Model<Scalar>* const model)
      : crocoddyl::ActionDataAbstractTpl<Scalar>(model) {}
};

/* --- Details -------------------------------------------------------------- */
/* --- Details -------------------------------------------------------------- */
/* --- Details -------------------------------------------------------------- */

typedef ActionModelQuadrupedAugmentedTimeTpl<double> ActionModelQuadrupedAugmentedTime;
typedef ActionDataQuadrupedAugmentedTimeTpl<double> ActionDataQuadrupedAugmentedTime;

}  // namespace quadruped_walkgen

#include "quadruped_augmented_time.hxx"

#endif
