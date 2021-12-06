#ifndef __quadruped_walkgen_quadruped_hpp__
#define __quadruped_walkgen_quadruped_hpp__
#include <stdexcept>

#include "crocoddyl/core/fwd.hpp"
#include "crocoddyl/core/action-base.hpp"
#include "crocoddyl/core/states/euclidean.hpp"
#include "crocoddyl/multibody/friction-cone.hpp"

#include "crocoddyl/core/utils/timer.hpp"

namespace quadruped_walkgen {
template <typename _Scalar>
class ActionModelQuadrupedTpl : public crocoddyl::ActionModelAbstractTpl<_Scalar> {
 public:
  typedef _Scalar Scalar;
  typedef crocoddyl::ActionDataAbstractTpl<Scalar> ActionDataAbstract;
  typedef crocoddyl::ActionModelAbstractTpl<Scalar> Base;
  typedef crocoddyl::MathBaseTpl<Scalar> MathBase;

  ActionModelQuadrupedTpl(typename Eigen::Matrix<Scalar, 3, 1> offset_CoM = Eigen::Matrix<Scalar, 3, 1>::Zero());
  ~ActionModelQuadrupedTpl();

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

  const Scalar& get_friction_weight() const;
  void set_friction_weight(const Scalar& weight);

  const Scalar& get_mu() const;
  void set_mu(const Scalar& mu_coeff);

  const Scalar& get_mass() const;
  void set_mass(const Scalar& m);

  const Scalar& get_dt() const;
  void set_dt(const Scalar& dt);

  const typename Eigen::Matrix<Scalar, 3, 3>& get_gI() const;
  void set_gI(const typename MathBase::Matrix3s& inertia_matrix);

  const Scalar& get_min_fz_contact() const;
  void set_min_fz_contact(const Scalar& min_fz);

  const Scalar& get_max_fz_contact() const;
  void set_max_fz_contact(const Scalar& max_fz_);

  // Set parameter relative to the shoulder height cost
  const Scalar& get_shoulder_hlim() const;
  void set_shoulder_hlim(const Scalar& hlim);

  const Scalar& get_shoulder_weight() const;
  void set_shoulder_weight(const Scalar& weight);

  const bool& get_relative_forces() const;
  void set_relative_forces(const bool& rel_forces);

  const bool& get_implicit_integration() const;
  void set_implicit_integration(const bool& implicit);

  // Update the model depending if the foot in contact with the ground
  // or the new lever arms
  void update_model(const Eigen::Ref<const typename MathBase::MatrixXs>& l_feet,
                    const Eigen::Ref<const typename MathBase::MatrixXs>& xref,
                    const Eigen::Ref<const typename MathBase::MatrixXs>& S);

  // Get A & B matrix
  const typename Eigen::Matrix<Scalar, 12, 12>& get_A() const;
  const typename Eigen::Matrix<Scalar, 12, 12>& get_B() const;

 protected:
  using Base::has_control_limits_;  //!< Indicates whether any of the control limits
  using Base::nr_;                  //!< Dimension of the cost residual
  using Base::nu_;                  //!< Control dimension
  using Base::state_;               //!< Model of the state
  using Base::u_lb_;                //!< Lower control limits
  using Base::u_ub_;                //!< Upper control limits
  using Base::unone_;               //!< Neutral state

 private:
  Scalar dt_;
  Scalar mass;
  Scalar mu;
  Scalar friction_weight_;
  Scalar min_fz_in_contact;
  Scalar max_fz;
  bool relative_forces;
  bool implicit_integration;

  typename Eigen::Matrix<Scalar, 12, 1> uref_;

  typename Eigen::Matrix<Scalar, 12, 1> force_weights_;
  typename Eigen::Matrix<Scalar, 12, 1> state_weights_;

  typename Eigen::Matrix<Scalar, 12, 12> A;
  typename Eigen::Matrix<Scalar, 12, 12> B;
  typename Eigen::Matrix<Scalar, 12, 1> g;
  typename Eigen::Matrix<Scalar, 3, 3> I_inv;
  typename MathBase::Matrix3s R_tmp;
  typename Eigen::Matrix<Scalar, 3, 3> gI;

  typename Eigen::Matrix<Scalar, 3, 4> lever_arms;
  typename MathBase::Vector3s lever_tmp;
  typename MathBase::MatrixXs xref_;

  typename Eigen::Matrix<Scalar, 24, 1> ub;

  typename Eigen::Matrix<Scalar, 24, 1> Fa_x_u;
  typename Eigen::Matrix<Scalar, 24, 1> rub_max_;
  typename Eigen::Matrix<Scalar, 24, 24> Arr;
  typename Eigen::Matrix<Scalar, 6, 1> r;

  // Cost relative to the shoulder height
  typename Eigen::Matrix<Scalar, 2, 4> pshoulder_0;
  typename Eigen::Matrix<Scalar, 3, 4> psh;
  typename Eigen::Matrix<Scalar, 4, 1> sh_ub_max_;
  typename Eigen::Matrix<Scalar, 4, 1> gait;
  typename Eigen::Matrix<Scalar, 3, 1> offset_com;
  Scalar sh_weight;
  Scalar sh_hlim;
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

/* --- Details -------------------------------------------------------------- */
/* --- Details -------------------------------------------------------------- */
/* --- Details -------------------------------------------------------------- */

typedef ActionModelQuadrupedTpl<double> ActionModelQuadruped;
typedef ActionDataQuadrupedTpl<double> ActionDataQuadruped;
typedef crocoddyl::ActionModelAbstractTpl<double> ActionModelAbstract;
typedef crocoddyl::ActionDataAbstractTpl<double> ActionDataAbstract;
typedef crocoddyl::StateAbstractTpl<double> StateAbstract;

typedef crocoddyl::ActionModelUnicycleTpl<double> ActionModelUnicycle;
typedef crocoddyl::ActionDataUnicycleTpl<double> ActionDataUnicycle;
}  // namespace quadruped_walkgen

#include "quadruped.hxx"

#endif
