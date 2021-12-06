///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (C) 2018-2020, LAAS-CNRS, University of Edinburgh
// Copyright note valid unless otherwise stated in individual files.
// All rights reserved.
///////////////////////////////////////////////////////////////////////////////

#ifndef BINDINGS_PYTHON_CROCODDYL_CORE_ACTION_BASE_HPP_
#define BINDINGS_PYTHON_CROCODDYL_CORE_ACTION_BASE_HPP_

#include <quadruped-walkgen/quadruped.hpp>
#include <quadruped-walkgen/quadruped_nl.hpp>
#include <quadruped-walkgen/quadruped_augmented.hpp>
#include <quadruped-walkgen/quadruped_augmented_time.hpp>
#include <quadruped-walkgen/quadruped_step.hpp>
#include <quadruped-walkgen/quadruped_time.hpp>
#include <quadruped-walkgen/quadruped_step_time.hpp>
#include <quadruped-walkgen/quadruped_step_period.hpp>
#include "crocoddyl/core/utils/exception.hpp"

namespace quadruped_walkgen {
namespace python {

class ActionModelAbstract_wrap : public ActionModelAbstract, public bp::wrapper<ActionModelAbstract> {
 public:
  ActionModelAbstract_wrap(boost::shared_ptr<StateAbstract> state, const std::size_t& nu, const std::size_t& nr = 1)
      : ActionModelAbstract(state, nu, nr), bp::wrapper<ActionModelAbstract>() {}

  void calc(const boost::shared_ptr<ActionDataAbstract>& data, const Eigen::Ref<const Eigen::VectorXd>& x,
            const Eigen::Ref<const Eigen::VectorXd>& u) {
    if (static_cast<std::size_t>(x.size()) != state_->get_nx()) {
      throw_pretty("Invalid argument: "
                   << "x has wrong dimension (it should be " + std::to_string(state_->get_nx()) + ")");
    }
    if (static_cast<std::size_t>(u.size()) != nu_) {
      throw_pretty("Invalid argument: "
                   << "u has wrong dimension (it should be " + std::to_string(nu_) + ")");
    }
    return bp::call<void>(this->get_override("calc").ptr(), data, (Eigen::VectorXd)x, (Eigen::VectorXd)u);
  }

  void calcDiff(const boost::shared_ptr<ActionDataAbstract>& data, const Eigen::Ref<const Eigen::VectorXd>& x,
                const Eigen::Ref<const Eigen::VectorXd>& u) {
    if (static_cast<std::size_t>(x.size()) != state_->get_nx()) {
      throw_pretty("Invalid argument: "
                   << "x has wrong dimension (it should be " + std::to_string(state_->get_nx()) + ")");
    }
    if (static_cast<std::size_t>(u.size()) != nu_) {
      throw_pretty("Invalid argument: "
                   << "u has wrong dimension (it should be " + std::to_string(nu_) + ")");
    }
    return bp::call<void>(this->get_override("calcDiff").ptr(), data, (Eigen::VectorXd)x, (Eigen::VectorXd)u);
  }
};

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(ActionModel_quasiStatic_wraps, ActionModelAbstract::quasiStatic_x, 2, 4)

}  // namespace python
}  // namespace quadruped_walkgen

#endif  // BINDINGS_PYTHON_CROCODDYL_CORE_ACTION_BASE_HPP_