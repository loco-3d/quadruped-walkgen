
#include "core.hpp"
#include "action-base.hpp"

namespace quadruped_walkgen {
namespace python {

void exposeActionQuadrupedAugmented() {
  bp::class_<ActionModelQuadrupedAugmented, bp::bases<ActionModelAbstract>>(
      "ActionModelQuadrupedAugmented",
      "Quadruped action model, non linear.\n\n"
      "The model is based on simplified dynamic model\n"
      "    xnext = Ax + Bu + g,\n"
      "The lever arms assumption is not take into account : the B matrix and its derivative depends on the state x\n"
      "where x i the state vector defined as : \n"
      "x = [x , y, z, rool, pitch, yaw, Vx, Vy, Vz, Vrool, Vpitch, Vyaw] , 12x \n\n"
      "and u is the groud reaction forces at each 4 foot, defined as : \n"
      "u = [fx1 , fy1, fz1, ... fz4], 12x",
      bp::init<Eigen::Matrix<double, 3, 1>>(bp::args("self", "offset_CoM"), "Initialize the quadruped action model."))
      .def("calc", &ActionModelQuadrupedAugmented::calc, bp::args("self", "data", "x", "u"),
           "Compute the next state and cost value.\n\n"
           "It describes the time-discrete evolution of the quadruped system.\n"
           "Additionally it computes the cost value associated to this discrete\n"
           "state and control pair.\n"
           ":param data: action data\n"
           ":param x: time-discrete state vector\n"
           ":param u: time-discrete control input")
      .def<void (ActionModelQuadrupedAugmented::*)(const boost::shared_ptr<ActionDataAbstract>&,
                                                   const Eigen::Ref<const Eigen::VectorXd>&)>(
          "calc", &ActionModelAbstract::calc, bp::args("self", "data", "x"))
      .def("calcDiff", &ActionModelQuadrupedAugmented::calcDiff, bp::args("self", "data", "x", "u"),
           "Compute the derivatives of the quadruped dynamics and cost functions.\n\n"
           "It computes the partial derivatives of the quadruped system and the\n"
           "cost function. It assumes that calc has been run first.\n"
           "This function builds a quadratic approximation of the\n"
           "action model (i.e. dynamical system and cost function).\n"
           ":param data: action data\n"
           ":param x: time-discrete state vector\n"
           ":param u: time-discrete control input\n")
      .def<void (ActionModelQuadrupedAugmented::*)(const boost::shared_ptr<ActionDataAbstract>&,
                                                   const Eigen::Ref<const Eigen::VectorXd>&)>(
          "calcDiff", &ActionModelAbstract::calcDiff, bp::args("self", "data", "x"))
      .def("createData", &ActionModelQuadrupedAugmented::createData, bp::args("self"),
           "Create the quadruped action data.")
      .def("updateModel", &ActionModelQuadrupedAugmented::update_model,
           bp::args("self", "l_feet", "l_stop", "xref", "S"),
           "Update the quadruped model depending on the position of the foot in the local frame\n\n"
           ":param l_feet : 3x4, Matrix representing the position of the foot in the local frame \n "
           "                Each colum represents the position of one foot : x,y,z"
           ":param xref : 12x1, Vector representing the reference state."
           ":param S : 4x1, Vector representing the foot in contact with the ground."
           "                S = [1 0 0 1] --> Foot 1 and 4 in contact.")
      .add_property(
          "forceWeights",
          bp::make_function(&ActionModelQuadrupedAugmented::get_force_weights, bp::return_internal_reference<>()),
          bp::make_function(&ActionModelQuadrupedAugmented::set_force_weights),
          "Weights on the control input : ground reaction forces")
      .add_property(
          "stateWeights",
          bp::make_function(&ActionModelQuadrupedAugmented::get_state_weights, bp::return_internal_reference<>()),
          bp::make_function(&ActionModelQuadrupedAugmented::set_state_weights), "Weights on the state vector")
      .add_property(
          "heuristicWeights",
          bp::make_function(&ActionModelQuadrupedAugmented::get_heuristic_weights, bp::return_internal_reference<>()),
          bp::make_function(&ActionModelQuadrupedAugmented::set_heuristic_weights), "Weights on the heuristic term")
      .add_property(
          "stopWeights",
          bp::make_function(&ActionModelQuadrupedAugmented::get_stop_weights, bp::return_internal_reference<>()),
          bp::make_function(&ActionModelQuadrupedAugmented::set_stop_weights), "Weights on the last position term")
      .add_property("frictionWeights",
                    bp::make_function(&ActionModelQuadrupedAugmented::get_friction_weight,
                                      bp::return_value_policy<bp::return_by_value>()),
                    bp::make_function(&ActionModelQuadrupedAugmented::set_friction_weight),
                    "Weight on friction cone term")
      .add_property(
          "mu",
          bp::make_function(&ActionModelQuadrupedAugmented::get_mu, bp::return_value_policy<bp::return_by_value>()),
          bp::make_function(&ActionModelQuadrupedAugmented::set_mu), "Friction coefficient")
      .add_property(
          "mass",
          bp::make_function(&ActionModelQuadrupedAugmented::get_mass, bp::return_value_policy<bp::return_by_value>()),
          bp::make_function(&ActionModelQuadrupedAugmented::set_mass),
          "Mass \n Warning : The model needs to be updated")
      .add_property(
          "dt",
          bp::make_function(&ActionModelQuadrupedAugmented::get_dt, bp::return_value_policy<bp::return_by_value>()),
          bp::make_function(&ActionModelQuadrupedAugmented::set_dt),
          "Minimum normal force allowed for feet in contact with the ground \n Warning : The model needs to be "
          "updated")
      .add_property("min_fz",
                    bp::make_function(&ActionModelQuadrupedAugmented::get_min_fz_contact,
                                      bp::return_value_policy<bp::return_by_value>()),
                    bp::make_function(&ActionModelQuadrupedAugmented::set_min_fz_contact),
                    "dt \n Warning : The model needs to be updated")
      .add_property("max_fz",
                    bp::make_function(&ActionModelQuadrupedAugmented::get_max_fz_contact,
                                      bp::return_value_policy<bp::return_by_value>()),
                    bp::make_function(&ActionModelQuadrupedAugmented::set_max_fz_contact),
                    "dt \n Warning : The model needs to be updated")
      .add_property(
          "gI",
          bp::make_function(&ActionModelQuadrupedAugmented::get_gI, bp::return_value_policy<bp::return_by_value>()),
          bp::make_function(&ActionModelQuadrupedAugmented::set_gI),
          "Inertia matrix of the robot in body frame (found in urdf) \n Warning : The model needs to be updated")
      .add_property("A", bp::make_function(&ActionModelQuadrupedAugmented::get_A, bp::return_internal_reference<>()),
                    "get A matrix")
      .add_property("B", bp::make_function(&ActionModelQuadrupedAugmented::get_B, bp::return_internal_reference<>()),
                    "get B matrix")
      .add_property("shoulder_hlim",
                    bp::make_function(&ActionModelQuadrupedAugmented::get_shoulder_hlim,
                                      bp::return_value_policy<bp::return_by_value>()),
                    bp::make_function(&ActionModelQuadrupedAugmented::set_shoulder_hlim), "Shoulder height limit ")
      .add_property("shoulderContactWeight",
                    bp::make_function(&ActionModelQuadrupedAugmented::get_shoulder_contact_weight,
                                      bp::return_value_policy<bp::return_by_value>()),
                    bp::make_function(&ActionModelQuadrupedAugmented::set_shoulder_contact_weight),
                    "shoulder Weights terms (array of size 8) ")
      .add_property("symmetry_term",
                    bp::make_function(&ActionModelQuadrupedAugmented::get_symmetry_term,
                                      bp::return_value_policy<bp::return_by_value>()),
                    bp::make_function(&ActionModelQuadrupedAugmented::set_symmetry_term),
                    "symmetry term for the foot position heuristic")
      .add_property("relative_forces",
                    bp::make_function(&ActionModelQuadrupedAugmented::get_relative_forces,
                                      bp::return_value_policy<bp::return_by_value>()),
                    bp::make_function(&ActionModelQuadrupedAugmented::set_relative_forces),
                    "activate relative forces ||fz-mg/nb_contact||^2")
      .add_property("centrifugal_term",
                    bp::make_function(&ActionModelQuadrupedAugmented::get_centrifugal_term,
                                      bp::return_value_policy<bp::return_by_value>()),
                    bp::make_function(&ActionModelQuadrupedAugmented::set_centrifugal_term),
                    "centrifugal term for the foot position heuristic")
      .add_property("T_gait",
                    bp::make_function(&ActionModelQuadrupedAugmented::get_T_gait,
                                      bp::return_value_policy<bp::return_by_value>()),
                    bp::make_function(&ActionModelQuadrupedAugmented::set_T_gait),
                    "Gait period, used to compute the symmetry term")
      .add_property("shoulderReferencePosition",
                    bp::make_function(&ActionModelQuadrupedAugmented::get_shoulder_reference_position,
                                      bp::return_value_policy<bp::return_by_value>()),
                    bp::make_function(&ActionModelQuadrupedAugmented::set_shoulder_reference_position),
                    "Bool, reference trajectory of the shoulder to compute the shoulder-to-contact cost");

  bp::register_ptr_to_python<boost::shared_ptr<ActionDataQuadrupedAugmented>>();

  bp::class_<ActionDataQuadrupedAugmented, bp::bases<ActionDataAbstract>>(
      "ActionDataQuadrupedAugmented",
      "Action data for the non linear quadruped system.\n\n"
      "The quadruped data, apart of common one, contains the cost residuals used\n"
      "for the computation of calc and calcDiff.",
      bp::init<ActionModelQuadrupedAugmented*>(bp::args("self", "model"),
                                               "Create quadruped data.\n\n"
                                               ":param model: quadruped action model"));
}

}  // namespace python
}  // namespace quadruped_walkgen
