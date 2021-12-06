///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (C) 2018-2019, LAAS-CNRS
// Copyright note valid unless otherwise stated in individual files.
// All rights reserved.
///////////////////////////////////////////////////////////////////////////////

#include "crocoddyl/core/actions/unicycle.hpp"
#include "crocoddyl/core/utils/callbacks.hpp"
#include "crocoddyl/core/solvers/ddp.hpp"
#include "crocoddyl/core/utils/timer.hpp"
#include <quadruped-walkgen/quadruped.hpp>

void updateModel(std::vector<boost::shared_ptr<quadruped_walkgen::ActionModelQuadruped> > running_models_2,
                 boost::shared_ptr<quadruped_walkgen::ActionModelQuadruped> terminal_model_2,
                 Eigen::Matrix<double, 6, 5> gait, Eigen::Matrix<double, 12, 17> xref,
                 Eigen::Matrix<double, 6, 13> fsteps, int N) {
  // Iterate over all the phases of the gait matrix
  // The first column of xref correspond to the current state = x0
  // Tmp is needed to use .data(), transformation of a column into a vector
  Eigen::Array<double, 1, 12> tmp = Eigen::Array<double, 1, 12>::Zero();
  int max_index = int(gait.block(0, 0, 6, 1).array().min(1.).matrix().sum());
  int k_cum = 0;

  for (int j = 0; j < max_index; j++) {
    for (int k = k_cum; k < k_cum + int(gait(j, 0)); ++k) {
      if (k < int(N)) {
        // Update model :
        tmp = fsteps.block(j, 1, 1, 12).array();
        running_models_2[k]->update_model(
            Eigen::Map<Eigen::Matrix<double, 3, 4> >(tmp.data(), 3, 4),
            Eigen::Map<Eigen::Matrix<double, 12, 1> >(xref.block(0, k + 1, 12, 1).data(), 12, 1),
            Eigen::Map<Eigen::Matrix<double, 4, 1> >(gait.block(j, 1, 1, 4).data(), 4, 1));
      }
    }
    k_cum += int(gait(j, 0));
  }

  tmp = fsteps.block(max_index - 1, 1, 1, 12).array();
  Eigen::Array<double, 1, 4> gait_tmp = Eigen::Array<double, 1, 4>::Zero();
  gait_tmp = gait.block(max_index - 1, 1, 1, 4).array();

  terminal_model_2->update_model(Eigen::Map<Eigen::Matrix<double, 3, 4> >(tmp.data(), 3, 4),
                                 Eigen::Map<Eigen::Matrix<double, 12, 1> >(xref.block(0, 16, 12, 1).data(), 12, 1),
                                 Eigen::Map<Eigen::Matrix<double, 4, 1> >(gait_tmp.data(), 4, 1));
  terminal_model_2->set_force_weights(Eigen::Matrix<double, 12, 1>::Zero());
  terminal_model_2->set_friction_weight(0);
}

int main(int argc, char* argv[]) {
  // The time of the cycle contol is 0.02s, and last 0.32s --> 16nodes
  // Control cycle during one gait period
  unsigned int N = 16;    // number of nodes
  unsigned int T = 1000;  // number of trials
  unsigned int MAXITER = 1;
  if (argc > 1) {
    T = atoi(argv[1]);
    MAXITER = atoi(argv[2]);
    ;
  }

  // Creating the initial state vector (size x12) [x,y,z,Roll,Pitch,Yaw,Vx,Vy,Vz,Wroll,Wpitch,Wyaw]
  // Perturbation of Vx = 0.2m.s-1
  Eigen::Matrix<double, 12, 1> x0;
  x0 << 0, 0, 0.2, 0, 0, 0, 0.2, 0, 0, 0, 0, 0;
  Eigen::Matrix<double, 4, 1> S;
  S << 1, 0, 0, 1;

  // Creating the reference state vector (size 12x16) to follow during the control cycle
  // Nullifying the Vx speed.
  Eigen::Matrix<double, 12, 1> xref_vector;
  xref_vector << 0, 0, 0.2, 0, 0, 0, 0, 0, 0, 0, 0, 0;
  Eigen::Matrix<double, 12, 17> xref;
  xref.block(0, 0, 12, 1) = x0;  // first vector is the initial state
  xref.block(0, 1, 12, 16) = xref_vector.replicate<1, 16>();

  // Creating the gait matrix : The number at the beginning represents the number of node spent in that position
  // 1 -> foot in contact with the ground :  0-> foot in the air
  Eigen::Matrix<double, 6, 5> gait;
  gait << 1, 1, 1, 1, 1, 7, 1, 0, 0, 1, 1, 1, 1, 1, 1, 7, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

  // Creating the fsteps matrix that represents the position of the feet during the whole control cycle (0.32s).
  // [nb , x1,y1,z1,  x2,y2,z2 ...] in local frame
  // The number at the beginning represents the number of node spent in that position
  // Here, the robot starts with 4 feet on the ground at the first node, then during 7 nodes (0.02s * 7)
  // The leg right front leg and left back leg are in the air ...etc

  Eigen::Matrix<double, 6, 13> fsteps;
  fsteps << 1, 0.19, 0.15, 0.0, 0.19, -0.15, 0.0, -0.19, 0.15, 0.0, -0.19, -0.15, 0.0, 7, 0.19, 0.15, 0.0, 0, 0, 0, 0,
      0, 0, -0.19, -0.15, 0.0, 1, 0.19, 0.15, 0.0, 0.19, -0.15, 0.0, -0.19, 0.15, 0.0, -0.19, -0.15, 0.0, 7, 0, 0, 0,
      0.19, -0.15, 0.0, -0.19, 0.15, 0.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0;

  // Creating the Shoting problem that needs boost::shared_ptr<crocoddyl::ActionModelAbstract>

  // Cannot use 1 model for the whole control cycle, because each model depends on the position of the feet
  // And the inertia matrix depends on the reference state (approximation )
  std::vector<boost::shared_ptr<crocoddyl::ActionModelAbstract> > running_models;

  for (int i = 0; i < int(N); ++i) {  // 16 nodes
    boost::shared_ptr<crocoddyl::ActionModelAbstract> model =
        boost::make_shared<quadruped_walkgen::ActionModelQuadruped>();
    running_models.push_back(model);
  }

  boost::shared_ptr<crocoddyl::ActionModelAbstract> terminal_model;
  terminal_model = boost::make_shared<quadruped_walkgen::ActionModelQuadruped>();

  // Update each model and set to 0 the weight ont the command for the terminal node
  // For that, the internal method of quadruped_walkgen::ActionModelQuadruped needs to be accessed
  // -> Creation of a 2nd list using dynamic_cast

  int k_cum = 0;
  std::vector<boost::shared_ptr<quadruped_walkgen::ActionModelQuadruped> > running_models_2;

  // Iterate over all the phases of the gait matrix
  // The first column of xref correspond to the current state = x0
  // Tmp is needed to use .data(), transformation of a column into a vector
  Eigen::Array<double, 1, 12> tmp = Eigen::Array<double, 1, 12>::Zero();
  int max_index = int(gait.block(0, 0, 6, 1).array().min(1.).matrix().sum());

  for (int j = 0; j < max_index; j++) {
    for (int k = k_cum; k < k_cum + int(gait(j, 0)); ++k) {
      if (k < int(N)) {
        boost::shared_ptr<quadruped_walkgen::ActionModelQuadruped> model2 =
            boost::dynamic_pointer_cast<quadruped_walkgen::ActionModelQuadruped>(running_models[k]);
        running_models_2.push_back(model2);

        // Update model :
        tmp = fsteps.block(j, 1, 1, 12).array();
        model2->update_model(Eigen::Map<Eigen::Matrix<double, 3, 4> >(tmp.data(), 3, 4),
                             Eigen::Map<Eigen::Matrix<double, 12, 1> >(xref.block(0, k + 1, 12, 1).data(), 12, 1),
                             Eigen::Map<Eigen::Matrix<double, 4, 1> >(gait.block(j, 1, 1, 4).data(), 4, 1));
      }
    }
    k_cum += int(gait(j, 0));
  }

  boost::shared_ptr<quadruped_walkgen::ActionModelQuadruped> terminal_model_2 =
      boost::dynamic_pointer_cast<quadruped_walkgen::ActionModelQuadruped>(terminal_model);

  tmp = fsteps.block(max_index - 1, 1, 1, 12).array();
  Eigen::Array<double, 1, 4> gait_tmp = Eigen::Array<double, 1, 4>::Zero();
  gait_tmp = gait.block(max_index - 1, 1, 1, 4).array();

  terminal_model_2->update_model(Eigen::Map<Eigen::Matrix<double, 3, 4> >(tmp.data(), 3, 4),
                                 Eigen::Map<Eigen::Matrix<double, 12, 1> >(xref.block(0, 16, 12, 1).data(), 12, 1),
                                 Eigen::Map<Eigen::Matrix<double, 4, 1> >(gait_tmp.data(), 4, 1));
  terminal_model_2->set_force_weights(Eigen::Matrix<double, 12, 1>::Zero());
  terminal_model_2->set_friction_weight(0);

  boost::shared_ptr<crocoddyl::ShootingProblem> problem =
      boost::make_shared<crocoddyl::ShootingProblem>(x0, running_models, terminal_model);
  crocoddyl::SolverDDP ddp(problem);

  std::vector<Eigen::VectorXd> xs(int(N) + 1, x0);
  Eigen::Matrix<double, 12, 1> u0;
  u0 << 1, 0.2, 0.5, 1, 1, -0.2, -1, 1, 0.5, -1, -1, -0.5;
  std::vector<Eigen::VectorXd> us(int(N), u0);

  Eigen::ArrayXd duration(T);

  // Updating the problem
  for (unsigned int i = 0; i < T; ++i) {
    crocoddyl::Timer timer;
    updateModel(running_models_2, terminal_model_2, gait, xref, fsteps, N);
    duration[i] = timer.get_duration();
  }
  double avrg_duration = duration.sum() / T;
  double min_duration = duration.minCoeff();
  double max_duration = duration.maxCoeff();
  std::cout << "  UpdateModel [ms]: " << avrg_duration << " (" << min_duration << "-" << max_duration << ")"
            << std::endl;

  // Solving the optimal control problem
  for (unsigned int i = 0; i < T; ++i) {
    crocoddyl::Timer timer;
    ddp.solve(xs, us, MAXITER);
    duration[i] = timer.get_duration();
  }

  avrg_duration = duration.sum() / T;
  min_duration = duration.minCoeff();
  max_duration = duration.maxCoeff();
  std::cout << "  DDP.solve [ms]: " << avrg_duration << " (" << min_duration << "-" << max_duration << ")"
            << std::endl;

  // Solving the optimal control problem
  for (unsigned int i = 0; i < T; ++i) {
    crocoddyl::Timer timer;
    ddp.solve(xs, us, MAXITER);
    duration[i] = timer.get_duration();
  }

  avrg_duration = duration.sum() / T;
  min_duration = duration.minCoeff();
  max_duration = duration.maxCoeff();
  std::cout << "  DDP.solve [ms]: " << avrg_duration << " (" << min_duration << "-" << max_duration << ")"
            << std::endl;

  // Running calc
  for (unsigned int i = 0; i < T; ++i) {
    crocoddyl::Timer timer;
    problem->calc(xs, us);
    duration[i] = timer.get_duration();
  }

  avrg_duration = duration.sum() / T;
  min_duration = duration.minCoeff();
  max_duration = duration.maxCoeff();
  std::cout << "  ShootingProblem.calc [ms]: " << avrg_duration << " (" << min_duration << "-" << max_duration << ")"
            << std::endl;

  // Running calcDiff
  for (unsigned int i = 0; i < T; ++i) {
    crocoddyl::Timer timer;
    problem->calcDiff(xs, us);
    duration[i] = timer.get_duration();
  }

  avrg_duration = duration.sum() / T;
  min_duration = duration.minCoeff();
  max_duration = duration.maxCoeff();
  std::cout << "  ShootingProblem.calcDiff [ms]: " << avrg_duration << " (" << min_duration << "-" << max_duration
            << ")" << std::endl;
}