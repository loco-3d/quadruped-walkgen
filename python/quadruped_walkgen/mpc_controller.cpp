#include <eigenpy/eigenpy.hpp>
#include <boost/python.hpp>
#include <boost/python/enum.hpp>

#include <quadruped-walkgen/quadruped.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "crocoddyl/core/action-base.hpp"
#include "crocoddyl/core/utils/exception.hpp"


namespace mpc_controller {
namespace python {

namespace bp = boost::python;

/////////////////////////////  WRAP FOR ActionModelAbstract ///////////////////////////////
////////////////////////////   to use bp:bases              ///////////////////////////////

class ActionModelAbstract_wrap : public crocoddyl::ActionModelAbstract, public bp::wrapper<crocoddyl::ActionModelAbstract> {
 public:
  ActionModelAbstract_wrap(boost::shared_ptr<crocoddyl::StateAbstract> state, const std::size_t& nu, const std::size_t& nr = 1)
      : crocoddyl::ActionModelAbstract(state, nu, nr), bp::wrapper<crocoddyl::ActionModelAbstract>() {}

  void calc(const boost::shared_ptr<crocoddyl::ActionDataAbstract>& data, const Eigen::Ref<const Eigen::VectorXd>& x,
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

  void calcDiff(const boost::shared_ptr<crocoddyl::ActionDataAbstract>& data, const Eigen::Ref<const Eigen::VectorXd>& x,
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

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(ActionModel_quasiStatic_wraps, crocoddyl::ActionModelAbstract::quasiStatic_x, 2, 4)

/////////////////////////////  bp::class_   for ActionModelAbstract ///////////////////////////////


BOOST_PYTHON_MODULE(libmpc_controller_pywrap) {

    bp::register_ptr_to_python<boost::shared_ptr<crocoddyl::ActionModelAbstract> >();

  bp::class_<ActionModelAbstract_wrap, boost::noncopyable>(
      "ActionModelAbstract",
      "Abstract class for action models.\n\n"
      "In crocoddyl, an action model combines dynamics and cost data. Each node, in our optimal\n"
      "control problem, is described through an action model. Every time that we want to describe\n"
      "a problem, we need to provide ways of computing the dynamics, cost functions and their\n"
      "derivatives. These computations are mainly carry on inside calc() and calcDiff(),\n"
      "respectively.",
      bp::init<boost::shared_ptr<crocoddyl::StateAbstract>, int, bp::optional<int> >(
          bp::args("self", "state", "nu", "nr"),
          "Initialize the action model.\n\n"
          "You can also describe autonomous systems by setting nu = 0.\n"
          ":param state: state description,\n"
          ":param nu: dimension of control vector,\n"
          ":param nr: dimension of the cost-residual vector (default 1)"))

      .def("calc", pure_virtual(&ActionModelAbstract_wrap::calc), bp::args("self", "data", "x", "u"),
           "Compute the next state and cost value.\n\n"
           "It describes the time-discrete evolution of our dynamical system\n"
           "in which we obtain the next state. Additionally it computes the\n"
           "cost value associated to this discrete state and control pair.\n"
           ":param data: action data\n"
           ":param x: time-discrete state vector\n"
           ":param u: time-discrete control input")
      
      .def("calcDiff", pure_virtual(&ActionModelAbstract_wrap::calcDiff), bp::args("self", "data", "x", "u"),
           "Compute the derivatives of the dynamics and cost functions.\n\n"
           "It computes the partial derivatives of the dynamical system and the\n"
           "cost function. It assumes that calc has been run first.\n"
           "This function builds a quadratic approximation of the\n"
           "action model (i.e. linear dynamics and quadratic cost).\n"
           ":param data: action data\n"
           ":param x: time-discrete state vector\n"
           ":param u: time-discrete control input\n")
      
      .def("createData", &ActionModelAbstract_wrap::createData, bp::args("self"),
           "Create the action data.\n\n"
           "Each action model (AM) has its own data that needs to be allocated.\n"
           "This function returns the allocated data for a predefined AM.\n"
           ":return AM data.")
      
      .add_property(
          "nu", bp::make_function(&ActionModelAbstract_wrap::get_nu, bp::return_value_policy<bp::return_by_value>()),
          "dimension of control vector")
      .add_property(
          "nr", bp::make_function(&ActionModelAbstract_wrap::get_nr, bp::return_value_policy<bp::return_by_value>()),
          "dimension of cost-residual vector")
      .add_property(
          "state",
          bp::make_function(&ActionModelAbstract_wrap::get_state, bp::return_value_policy<bp::return_by_value>()),
          "state")
      .add_property("has_control_limits",
                    bp::make_function(&ActionModelAbstract_wrap::get_has_control_limits,
                                      bp::return_value_policy<bp::return_by_value>()),
                    "indicates whether problem has finite control limits")
      .add_property("u_lb", bp::make_function(&ActionModelAbstract_wrap::get_u_lb, bp::return_internal_reference<>()),
                    &ActionModelAbstract_wrap::set_u_lb, "lower control limits")
      .add_property("u_ub", bp::make_function(&ActionModelAbstract_wrap::get_u_ub, bp::return_internal_reference<>()),
                    &ActionModelAbstract_wrap::set_u_ub, "upper control limits");

  bp::register_ptr_to_python<boost::shared_ptr<crocoddyl::ActionDataAbstract> >();

  bp::class_<crocoddyl::ActionDataAbstract, boost::noncopyable>(
      "ActionDataAbstract",
      "Abstract class for action data.\n\n"
      "In crocoddyl, an action data contains all the required information for processing an\n"
      "user-defined action model. The action data typically is allocated onces by running\n"
      "model.createData() and contains the first- and second- order derivatives of the dynamics\n"
      "and cost function, respectively.",
      bp::init<crocoddyl::ActionModelAbstract*>(bp::args("self", "model"),
                                     "Create common data shared between AMs.\n\n"
                                     "The action data uses the model in order to first process it.\n"
                                     ":param model: action model"))
      .add_property("cost", bp::make_getter(&crocoddyl::ActionDataAbstract::cost, bp::return_value_policy<bp::return_by_value>()),
                    bp::make_setter(&crocoddyl::ActionDataAbstract::cost), "cost value")
      .add_property("xnext", bp::make_getter(&crocoddyl::ActionDataAbstract::xnext, bp::return_internal_reference<>()),
                    bp::make_setter(&crocoddyl::ActionDataAbstract::xnext), "next state")
      .add_property("r", bp::make_getter(&crocoddyl::ActionDataAbstract::r, bp::return_internal_reference<>()),
                    bp::make_setter(&crocoddyl::ActionDataAbstract::r), "cost residual")
      .add_property("Fx", bp::make_getter(&crocoddyl::ActionDataAbstract::Fx, bp::return_internal_reference<>()),
                    bp::make_setter(&crocoddyl::ActionDataAbstract::Fx), "Jacobian of the dynamics")
      .add_property("Fu", bp::make_getter(&crocoddyl::ActionDataAbstract::Fu, bp::return_internal_reference<>()),
                    bp::make_setter(&crocoddyl::ActionDataAbstract::Fu), "Jacobian of the dynamics")
      .add_property("Lx", bp::make_getter(&crocoddyl::ActionDataAbstract::Lx, bp::return_internal_reference<>()),
                    bp::make_setter(&crocoddyl::ActionDataAbstract::Lx), "Jacobian of the cost")
      .add_property("Lu", bp::make_getter(&crocoddyl::ActionDataAbstract::Lu, bp::return_internal_reference<>()),
                    bp::make_setter(&crocoddyl::ActionDataAbstract::Lu), "Jacobian of the cost")
      .add_property("Lxx", bp::make_getter(&crocoddyl::ActionDataAbstract::Lxx, bp::return_internal_reference<>()),
                    bp::make_setter(&crocoddyl::ActionDataAbstract::Lxx), "Hessian of the cost")
      .add_property("Lxu", bp::make_getter(&crocoddyl::ActionDataAbstract::Lxu, bp::return_internal_reference<>()),
                    bp::make_setter(&crocoddyl::ActionDataAbstract::Lxu), "Hessian of the cost")
      .add_property("Luu", bp::make_getter(&crocoddyl::ActionDataAbstract::Luu, bp::return_internal_reference<>()),
                    bp::make_setter(&crocoddyl::ActionDataAbstract::Luu), "Hessian of the cost");

    
    bp::register_ptr_to_python<boost::shared_ptr<ActionModelQuadruped> >();
    //////////////////////////////////////  bp:class_ ActionModelQuadruped ///////////////////////////////////////

    bp::class_<ActionModelQuadruped, bp::bases<crocoddyl::ActionModelAbstract> >(
      "ActionModelQuadruped",
      "Unicycle action model.\n\n"
      "The transition model of an unicycle system is described as\n"
      "    xnext = [v*cos(theta); v*sin(theta); w],\n"
      "where the position is defined by (x, y, theta) and the control input\n"
      "by (v,w). Note that the state is defined only with the position. On the\n"
      "other hand, we define the quadratic cost functions for the state and\n"
      "control.",
      bp::init<>(bp::args("self"), "Initialize the unicycle action model."))

      .def("calc", &ActionModelQuadruped::calc, bp::args("self", "data", "x", "u"),
          "Compute the next state and cost value.\n\n"
          "It describes the time-discrete evolution of the unicycle system.\n"
          "Additionally it computes the cost value associated to this discrete\n"
          "state and control pair.\n"
          ":param data: action data\n"
          ":param x: time-discrete state vector\n"
          ":param u: time-discrete control input")
     
      .def("calcDiff", &ActionModelQuadruped::calcDiff, bp::args("self", "data", "x", "u"),
          "Compute the derivatives of the unicycle dynamics and cost functions.\n\n"
          "It computes the partial derivatives of the unicycle system and the\n"
          "cost function. It assumes that calc has been run first.\n"
          "This function builds a quadratic approximation of the\n"
          "action model (i.e. dynamical system and cost function).\n"
          ":param data: action data\n"
          ":param x: time-discrete state vector\n"
          ":param u: time-discrete control input\n")

      .def("createData", &ActionModelQuadruped::createData, bp::args("self"), "Create the unicycle action data.")
      .add_property("costWeights",
                    bp::make_function(&ActionModelQuadruped::get_cost_weights, bp::return_internal_reference<>()),
                    bp::make_function(&ActionModelQuadruped::set_cost_weights), "cost weights");

  bp::register_ptr_to_python<boost::shared_ptr<ActionDataQuadruped> >();

  bp::class_<ActionDataQuadruped, bp::bases<crocoddyl::ActionDataAbstract> >(
      "ActionDataUnicycle",
      "Action data for the Unicycle system.\n\n"
      "The unicycle data, apart of common one, contains the cost residuals used\n"
      "for the computation of calc and calcDiff.",
      bp::init<ActionModelQuadruped*>(bp::args("self", "model"),
                                     "Create unicycle data.\n\n"
                                     ":param model: unicycle action model"));





}
//////////////////////////////////////////////////////////////////////////////////////////////

 

}
}
