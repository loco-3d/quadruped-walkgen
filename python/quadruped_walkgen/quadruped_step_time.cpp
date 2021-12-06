
#include "core.hpp"
#include "action-base.hpp"

namespace quadruped_walkgen {
namespace python {

void exposeActionQuadrupedStepTime() {
  bp::class_<ActionModelQuadrupedStepTime, bp::bases<ActionModelAbstract> >(
      "ActionModelQuadrupedStepTime",
      "Quadruped action model, non linear.\n\n"
      "The model is based on simplified dynamic model\n"
      "    xnext = Ax + Bu + g,\n"
      "The lever arms assumption is not take into account : the B matrix and its derivative depends on the state x\n"
      "where x i the state vector defined as : \n"
      "x = [x , y, z, rool, pitch, yaw, Vx, Vy, Vz, Vrool, Vpitch, Vyaw] , 12x \n\n"
      "and u is the groud reaction forces at each 4 foot, defined as : \n"
      "u = [fx1 , fy1, fz1, ... fz4], 12x",
      bp::init<>(bp::args("self"), "Initialize the quadruped action model."))
      .def("calc", &ActionModelQuadrupedStepTime::calc, bp::args("self", "data", "x", "u"),
           "Compute the next state and cost value.\n\n"
           "It describes the time-discrete evolution of the quadruped system.\n"
           "Additionally it computes the cost value associated to this discrete\n"
           "state and control pair.\n"
           ":param data: action data\n"
           ":param x: time-discrete state vector\n"
           ":param u: time-discrete control input")
      .def<void (ActionModelQuadrupedStepTime::*)(const boost::shared_ptr<ActionDataAbstract>&,
                                                  const Eigen::Ref<const Eigen::VectorXd>&)>(
          "calc", &ActionModelAbstract::calc, bp::args("self", "data", "x"))
      .def("calcDiff", &ActionModelQuadrupedStepTime::calcDiff, bp::args("self", "data", "x", "u"),
           "Compute the derivatives of the quadruped dynamics and cost functions.\n\n"
           "It computes the partial derivatives of the quadruped system and the\n"
           "cost function. It assumes that calc has been run first.\n"
           "This function builds a quadratic approximation of the\n"
           "action model (i.e. dynamical system and cost function).\n"
           ":param data: action data\n"
           ":param x: time-discrete state vector\n"
           ":param u: time-discrete control input\n")
      .def<void (ActionModelQuadrupedStepTime::*)(const boost::shared_ptr<ActionDataAbstract>&,
                                                  const Eigen::Ref<const Eigen::VectorXd>&)>(
          "calcDiff", &ActionModelAbstract::calcDiff, bp::args("self", "data", "x"))
      .def("createData", &ActionModelQuadrupedStepTime::createData, bp::args("self"),
           "Create the quadruped action data.")
      .def("updateModel", &ActionModelQuadrupedStepTime::update_model,
           bp::args("self", "l_feet", "velocity", "acceleration", "xref", "S"),
           "Update the quadruped model depending on the position of the foot in the local frame\n\n"
           ":param l_feet : 3x4, Matrix representing the position of the foot in the local frame \n "
           "                Each colum represents the position of one foot : x,y,z"
           ":param xref : 12x1, Vector representing the reference state."
           ":param S : 4x1, Vector representing the foot in contact with the ground."
           "                S = [1 0 0 1] --> Foot 1 and 4 in contact.")
      .add_property(
          "stateWeights",
          bp::make_function(&ActionModelQuadrupedStepTime::get_state_weights, bp::return_internal_reference<>()),
          bp::make_function(&ActionModelQuadrupedStepTime::set_state_weights), "Weights on the state vector")
      .add_property(
          "heuristicWeights",
          bp::make_function(&ActionModelQuadrupedStepTime::get_heuristic_weights, bp::return_internal_reference<>()),
          bp::make_function(&ActionModelQuadrupedStepTime::set_heuristic_weights), "Weights on the shoulder term")
      .add_property(
          "stepWeights",
          bp::make_function(&ActionModelQuadrupedStepTime::get_step_weights, bp::return_internal_reference<>()),
          bp::make_function(&ActionModelQuadrupedStepTime::set_step_weights), "Weights on the command norm")
      .add_property("symmetry_term",
                    bp::make_function(&ActionModelQuadrupedStepTime::get_symmetry_term,
                                      bp::return_value_policy<bp::return_by_value>()),
                    bp::make_function(&ActionModelQuadrupedStepTime::set_symmetry_term),
                    "symmetry term for the foot position heuristic")
      .add_property("centrifugal_term",
                    bp::make_function(&ActionModelQuadrupedStepTime::get_centrifugal_term,
                                      bp::return_value_policy<bp::return_by_value>()),
                    bp::make_function(&ActionModelQuadrupedStepTime::set_centrifugal_term),
                    "centrifugal term for the foot position heuristic")
      .add_property(
          "T_gait",
          bp::make_function(&ActionModelQuadrupedStepTime::get_T_gait, bp::return_value_policy<bp::return_by_value>()),
          bp::make_function(&ActionModelQuadrupedStepTime::set_T_gait),
          "Gait period, used to compute the symmetry term")
      .add_property("speed_weight",
                    bp::make_function(&ActionModelQuadrupedStepTime::get_speed_weight,
                                      bp::return_value_policy<bp::return_by_value>()),
                    bp::make_function(&ActionModelQuadrupedStepTime::set_speed_weight),
                    "Gait period, used to compute the symmetry term")
      .add_property(
          "vlim",
          bp::make_function(&ActionModelQuadrupedStepTime::get_vlim, bp::return_value_policy<bp::return_by_value>()),
          bp::make_function(&ActionModelQuadrupedStepTime::set_vlim), "Gait period, used to compute the symmetry term")
      .add_property("nb_nodes",
                    bp::make_function(&ActionModelQuadrupedStepTime::get_nb_nodes,
                                      bp::return_value_policy<bp::return_by_value>()),
                    bp::make_function(&ActionModelQuadrupedStepTime::set_nb_nodes),
                    "Gait period, used to compute the symmetry term")
      .add_property("first_step",
                    bp::make_function(&ActionModelQuadrupedStepTime::get_first_step,
                                      bp::return_value_policy<bp::return_by_value>()),
                    bp::make_function(&ActionModelQuadrupedStepTime::set_first_step),
                    "bool thats indicates xhether it is the first step foot, otherwise function cost is much simpler")
      .add_property("Cost",
                    bp::make_function(&ActionModelQuadrupedStepTime::get_cost, bp::return_internal_reference<>()),
                    "get log cost");

  bp::register_ptr_to_python<boost::shared_ptr<ActionDataQuadrupedStepTime> >();

  bp::class_<ActionDataQuadrupedStepTime, bp::bases<ActionDataAbstract> >(
      "ActionDataQuadrupedStepTime",
      "Action data for the non linear quadruped system.\n\n"
      "The quadruped data, apart of common one, contains the cost residuals used\n"
      "for the computation of calc and calcDiff.",
      bp::init<ActionModelQuadrupedStepTime*>(bp::args("self", "model"),
                                              "Create quadruped data.\n\n"
                                              ":param model: quadruped action model"));
}

}  // namespace python
}  // namespace quadruped_walkgen
