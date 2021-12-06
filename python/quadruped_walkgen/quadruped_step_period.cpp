
#include "core.hpp"
#include "action-base.hpp"

namespace quadruped_walkgen {
namespace python {

void exposeActionQuadrupedStepPeriod() {
  bp::class_<ActionModelQuadrupedStepPeriod, bp::bases<ActionModelAbstract> >(
      "ActionModelQuadrupedStepPeriod",
      "Quadruped action model, non linear.\n\n"
      "The model is based on simplified dynamic model\n"
      "    xnext = Ax + Bu + g,\n"
      "The lever arms assumption is not take into account : the B matrix and its derivative depends on the state x\n"
      "where x i the state vector defined as : \n"
      "x = [x , y, z, rool, pitch, yaw, Vx, Vy, Vz, Vrool, Vpitch, Vyaw] , 12x \n\n"
      "and u is the groud reaction forces at each 4 foot, defined as : \n"
      "u = [fx1 , fy1, fz1, ... fz4], 12x",
      bp::init<>(bp::args("self"), "Initialize the quadruped action model."))
      .def("calc", &ActionModelQuadrupedStepPeriod::calc, bp::args("self", "data", "x", "u"),
           "Compute the next state and cost value.\n\n"
           "It describes the time-discrete evolution of the quadruped system.\n"
           "Additionally it computes the cost value associated to this discrete\n"
           "state and control pair.\n"
           ":param data: action data\n"
           ":param x: time-discrete state vector\n"
           ":param u: time-discrete control input")
      .def<void (ActionModelQuadrupedStepPeriod::*)(const boost::shared_ptr<ActionDataAbstract>&,
                                                    const Eigen::Ref<const Eigen::VectorXd>&)>(
          "calc", &ActionModelAbstract::calc, bp::args("self", "data", "x"))
      .def("calcDiff", &ActionModelQuadrupedStepPeriod::calcDiff, bp::args("self", "data", "x", "u"),
           "Compute the derivatives of the quadruped dynamics and cost functions.\n\n"
           "It computes the partial derivatives of the quadruped system and the\n"
           "cost function. It assumes that calc has been run first.\n"
           "This function builds a quadratic approximation of the\n"
           "action model (i.e. dynamical system and cost function).\n"
           ":param data: action data\n"
           ":param x: time-discrete state vector\n"
           ":param u: time-discrete control input\n")
      .def<void (ActionModelQuadrupedStepPeriod::*)(const boost::shared_ptr<ActionDataAbstract>&,
                                                    const Eigen::Ref<const Eigen::VectorXd>&)>(
          "calcDiff", &ActionModelAbstract::calcDiff, bp::args("self", "data", "x"))
      .def("createData", &ActionModelQuadrupedStepPeriod::createData, bp::args("self"),
           "Create the quadruped action data.")
      .def("updateModel", &ActionModelQuadrupedStepPeriod::update_model, bp::args("self", "l_feet", "xref", "S"),
           "Update the quadruped model depending on the position of the foot in the local frame\n\n"
           ":param l_feet : 3x4, Matrix representing the position of the foot in the local frame \n "
           "                Each colum represents the position of one foot : x,y,z"
           ":param xref : 12x1, Vector representing the reference state."
           ":param S : 4x1, Vector representing the foot in contact with the ground."
           "                S = [1 0 0 1] --> Foot 1 and 4 in contact.")
      .add_property(
          "stateWeights",
          bp::make_function(&ActionModelQuadrupedStepPeriod::get_state_weights, bp::return_internal_reference<>()),
          bp::make_function(&ActionModelQuadrupedStepPeriod::set_state_weights), "Weights on the state vector")
      .add_property(
          "shoulderWeights",
          bp::make_function(&ActionModelQuadrupedStepPeriod::get_shoulder_weights, bp::return_internal_reference<>()),
          bp::make_function(&ActionModelQuadrupedStepPeriod::set_shoulder_weights), "Weights on the shoulder term")
      .add_property(
          "stepWeights",
          bp::make_function(&ActionModelQuadrupedStepPeriod::get_step_weights, bp::return_internal_reference<>()),
          bp::make_function(&ActionModelQuadrupedStepPeriod::set_step_weights), "Weights on the command norm")
      .add_property("symmetry_term",
                    bp::make_function(&ActionModelQuadrupedStepPeriod::get_symmetry_term,
                                      bp::return_value_policy<bp::return_by_value>()),
                    bp::make_function(&ActionModelQuadrupedStepPeriod::set_symmetry_term),
                    "symmetry term for the foot position heuristic")
      .add_property("centrifugal_term",
                    bp::make_function(&ActionModelQuadrupedStepPeriod::get_centrifugal_term,
                                      bp::return_value_policy<bp::return_by_value>()),
                    bp::make_function(&ActionModelQuadrupedStepPeriod::set_centrifugal_term),
                    "centrifugal term for the foot position heuristic")
      .add_property("T_gait",
                    bp::make_function(&ActionModelQuadrupedStepPeriod::get_T_gait,
                                      bp::return_value_policy<bp::return_by_value>()),
                    bp::make_function(&ActionModelQuadrupedStepPeriod::set_T_gait),
                    "Gait period, used to compute the symmetry term")
      .add_property("dt_ref",
                    bp::make_function(&ActionModelQuadrupedStepPeriod::get_dt_ref,
                                      bp::return_value_policy<bp::return_by_value>()),
                    bp::make_function(&ActionModelQuadrupedStepPeriod::set_dt_ref),
                    "dt_ref used on the dt reference : || dt - dt_ref||^2 ")
      .add_property("dt_min",
                    bp::make_function(&ActionModelQuadrupedStepPeriod::get_dt_min,
                                      bp::return_value_policy<bp::return_by_value>()),
                    bp::make_function(&ActionModelQuadrupedStepPeriod::set_dt_min), "Minimum dt")
      .add_property("dt_max",
                    bp::make_function(&ActionModelQuadrupedStepPeriod::get_dt_max,
                                      bp::return_value_policy<bp::return_by_value>()),
                    bp::make_function(&ActionModelQuadrupedStepPeriod::set_dt_max), "Maximum dt")
      .add_property("dt_weight",
                    bp::make_function(&ActionModelQuadrupedStepPeriod::get_dt_weight,
                                      bp::return_value_policy<bp::return_by_value>()),
                    bp::make_function(&ActionModelQuadrupedStepPeriod::set_dt_weight),
                    "Weight on the dt reference : || dt - dt_ref||^2 ")
      .add_property("dt_weight_bound",
                    bp::make_function(&ActionModelQuadrupedStepPeriod::get_dt_bound_weight,
                                      bp::return_value_policy<bp::return_by_value>()),
                    bp::make_function(&ActionModelQuadrupedStepPeriod::set_dt_bound_weight),
                    "Weight on the dt limit : || dt_min - dt||^2 etc..")
      .add_property("nb_nodes",
                    bp::make_function(&ActionModelQuadrupedStepPeriod::get_nb_nodes,
                                      bp::return_value_policy<bp::return_by_value>()),
                    bp::make_function(&ActionModelQuadrupedStepPeriod::set_nb_nodes),
                    "Nb of nodes during flying phase")
      .add_property(
          "vlim",
          bp::make_function(&ActionModelQuadrupedStepPeriod::get_vlim, bp::return_value_policy<bp::return_by_value>()),
          bp::make_function(&ActionModelQuadrupedStepPeriod::set_vlim), "Speed limit of the x and y foot trajectory")
      .add_property("speed_weight",
                    bp::make_function(&ActionModelQuadrupedStepPeriod::get_speed_weight,
                                      bp::return_value_policy<bp::return_by_value>()),
                    bp::make_function(&ActionModelQuadrupedStepPeriod::set_speed_weight),
                    "Speed limit of the x and y foot trajectory");

  bp::register_ptr_to_python<boost::shared_ptr<ActionDataQuadrupedStepPeriod> >();

  bp::class_<ActionDataQuadrupedStepPeriod, bp::bases<ActionDataAbstract> >(
      "ActionDataQuadrupedStepPeriod",
      "Action data for the non linear quadruped system.\n\n"
      "The quadruped data, apart of common one, contains the cost residuals used\n"
      "for the computation of calc and calcDiff.",
      bp::init<ActionModelQuadrupedStepPeriod*>(bp::args("self", "model"),
                                                "Create quadruped data.\n\n"
                                                ":param model: quadruped action model"));
}

}  // namespace python
}  // namespace quadruped_walkgen
