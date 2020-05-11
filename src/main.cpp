#include <pinocchio/parsers/urdf.hpp>
#include <example-robot-data/path.hpp>

int main(){

  // Load Pinocchio model
  pinocchio::Model model;
  pinocchio::urdf::buildReducedModel(EXAMPLE_ROBOT_DATA_MODEL_DIR "/talos_data/robots/talos_reduced.urdf", model);
