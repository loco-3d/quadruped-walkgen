///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (C) 2018-2019, LAAS-CNRS
// Copyright note valid unless otherwise stated in individual files.
// All rights reserved.
///////////////////////////////////////////////////////////////////////////////

// #include "crocoddyl/core/codegen/action-base.hpp"

#include "crocoddyl/core/actions/unicycle.hpp"
#include "crocoddyl/core/utils/callbacks.hpp"
#include "crocoddyl/core/solvers/ddp.hpp"
#include "crocoddyl/core/utils/timer.hpp"
#include <quadruped-walkgen/quadruped_augmented.hpp>
#include <quadruped-walkgen/quadruped_augmented_time.hpp>
#include <quadruped-walkgen/quadruped_step.hpp>
#include <quadruped-walkgen/quadruped_time.hpp>
#include <quadruped-walkgen/quadruped_step_time.hpp>

#define STDDEV(vec) std::sqrt(((vec - vec.mean())).square().sum() / (double(vec.size()) - 1.))
#define AVG(vec) (vec.mean())

int main(int argc, char* argv[]) {
  // The time of the cycle contol is 0.02s, and last 0.32s --> 16nodes
  // Control cycle during one gait period
  unsigned int N = 16;     // number of nodes
  unsigned int T = 20000;  // number of trials
  unsigned int MAXITER = 1;
  if (argc > 1) {
    T = atoi(argv[1]);
    MAXITER = atoi(argv[2]);
    ;
  }

  boost::shared_ptr<crocoddyl::ActionModelAbstract> model_test;
  model_test = boost::make_shared<quadruped_walkgen::ActionModelQuadrupedStepTime>();

  // Creating the initial state vector (size x12) [x,y,z,Roll,Pitch,Yaw,Vx,Vy,Vz,Wroll,Wpitch,Wyaw]
  // Perturbation of Vx = 0.2m.s-1
  Eigen::Matrix<double, 21, 1> x0;
  x0 << 0., 0., 0.2, 0., 0., 0., 0.2, 0., 0., 0., 0., 0., 0.1946, 0.15005, 0.204, -0.137, -0.184, 0.14, -0.1946,
      -0.1505, 0.02;
  Eigen::Matrix<double, 4, 1> S;
  S << 0, 1, 1, 0;

  Eigen::Matrix<double, 3, 4> l_feet;  // computed by previous gait cycle
  l_feet << 0.1946, 0.21, -0.18, -0.19, 0.15, -0.16, 0.145, -0.135, 0.0, 0.0, 0.0, 0.0;

  // Creating the reference state vector (size 12x16) to follow during the control cycle
  // Nullifying the Vx speed.
  Eigen::Matrix<double, 12, 1> xref_vector;
  xref_vector << 0., 0., 0.2, 0., 0., 0., 0., 0., 0., 0., 0., 0.;
  Eigen::Matrix<double, 12, 17> xref;
  xref.block(0, 0, 12, 1) << 0., 0., 0.2, 0., 0., 0., 0.1, 0., 0., 0., 0., 0.;  // first vector is the initial state
  xref.block(0, 1, 12, 16) = xref_vector.replicate<1, 16>();

  // Creating the gait matrix : The number at the beginning represents the number of node spent in that position
  // 1 -> foot in contact with the ground :  0-> foot in the air
  Eigen::Matrix<double, 6, 5> gait;
  gait << 7, 0, 1, 1, 0, 1, 1, 1, 1, 1, 7, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

  // Creating the Shoting problem that needs boost::shared_ptr<crocoddyl::ActionModelAbstract>

  // Cannot use 1 model for the whole control cycle, because each model depends on the position of the feet
  // And the inertia matrix depends on the reference state (approximation )
  std::vector<boost::shared_ptr<crocoddyl::ActionModelAbstract> > running_models;

  int max_index = int(gait.block(0, 0, 6, 1).array().min(1.).matrix().sum());
  int k_cum = 0;
  for (int j = 0; j < max_index; j++) {
    for (int k = k_cum; k < k_cum + int(gait(j, 0)); ++k) {
      if (k < int(N)) {
        if (int(gait.block(j, 1, 1, 4).sum()) == 4) {
          boost::shared_ptr<crocoddyl::ActionModelAbstract> model =
              boost::make_shared<quadruped_walkgen::ActionModelQuadrupedStepTime>();
          running_models.push_back(model);
        }
        if (j == 0 and k == 1) {
          boost::shared_ptr<crocoddyl::ActionModelAbstract> model =
              boost::make_shared<quadruped_walkgen::ActionModelQuadrupedTime>();
          running_models.push_back(model);
          std::cout << "okok1" << std::endl;
        }
        // if ( j == 2 and k == 9){
        //   boost::shared_ptr<crocoddyl::ActionModelAbstract> model
        //              = boost::make_shared<quadruped_walkgen::ActionModelQuadrupedTime>() ;
        //   running_models.push_back(model) ;
        //   std::cout<<"okok2"<<std::endl ;
        // }
        boost::shared_ptr<crocoddyl::ActionModelAbstract> model =
            boost::make_shared<quadruped_walkgen::ActionModelQuadrupedAugmentedTime>();
        running_models.push_back(model);
      }
    }
    k_cum += int(gait(j, 0));
  }

  boost::shared_ptr<crocoddyl::ActionModelAbstract> terminal_model;
  terminal_model = boost::make_shared<quadruped_walkgen::ActionModelQuadrupedAugmentedTime>();

  // Update each model and set to 0 the weight ont the command for the terminal node
  // For that, the internal method of quadruped_walkgen::ActionModelQuadruped needs to be accessed
  // -> Creation of a 2nd list using dynamic_cast

  k_cum = 0;
  Eigen::Matrix<double, 12, 1> u0;
  u0 << 1, 0.2, 8, 1, 1, 8, -1, 1, 8, -1, -1, 8;
  std::vector<Eigen::VectorXd> us;
  Eigen::Matrix<double, 4, 1> u0_step;
  u0_step << 0.05, 0.01, 0.02, 0.06;
  Eigen::Matrix<double, 1, 1> u0_time;
  u0_time << 0.02;

  // Iterate over all the phases of the gait matrix
  // The first column of xref correspond to the current state = x0
  // Tmp is needed to use .data(), transformation of a column into a vector
  Eigen::Array<double, 3, 4> tmp = Eigen::Array<double, 3, 4>::Zero();
  Eigen::Matrix<double, 1, 4> S_tmp;
  S_tmp.setZero();

  int gap = 0;
  for (int j = 0; j < max_index; j++) {
    for (int k = k_cum; k < k_cum + int(gait(j, 0)); ++k) {
      std::cout << k << std::endl;
      if (k < int(N)) {
        if (int(gait.block(j, 1, 1, 4).sum()) == 4) {
          boost::shared_ptr<quadruped_walkgen::ActionModelQuadrupedStepTime> model3 =
              boost::dynamic_pointer_cast<quadruped_walkgen::ActionModelQuadrupedStepTime>(running_models[k + gap]);

          tmp = l_feet.array();
          S_tmp = gait.block(j, 1, 1, 4) - gait.block(j - 1, 1, 1, 4);
          model3->update_model(Eigen::Map<Eigen::Matrix<double, 3, 4> >(tmp.data(), 3, 4),
                               Eigen::Map<Eigen::Matrix<double, 3, 4> >(tmp.data(), 3, 4),
                               Eigen::Map<Eigen::Matrix<double, 3, 4> >(tmp.data(), 3, 4),
                               Eigen::Map<Eigen::Matrix<double, 12, 1> >(xref.block(0, k, 12, 1).data(), 12, 1),
                               Eigen::Map<Eigen::Matrix<double, 4, 1> >(S_tmp.data(), 4, 1));

          gap = gap + 1;
          us.push_back(u0_step);
        }

        if (j == 0 and k == 1) {
          boost::shared_ptr<quadruped_walkgen::ActionModelQuadrupedTime> model1 =
              boost::dynamic_pointer_cast<quadruped_walkgen::ActionModelQuadrupedTime>(running_models[k + gap]);

          tmp = l_feet.array();
          S_tmp = gait.block(j, 1, 1, 4) - gait.block(j - 1, 1, 1, 4);
          model1->update_model(Eigen::Map<Eigen::Matrix<double, 3, 4> >(tmp.data(), 3, 4),
                               Eigen::Map<Eigen::Matrix<double, 12, 1> >(xref.block(0, k, 12, 1).data(), 12, 1),
                               Eigen::Map<Eigen::Matrix<double, 4, 1> >(S_tmp.data(), 4, 1));

          std::cout << "ok1" << std::endl;
          gap = gap + 1;
          us.push_back(u0_time);
        }
        // if ( j == 2 and k == 9){
        //   std::cout << "ok2" << std::endl ;
        //   gap = gap + 1 ;
        //   us.push_back(u0_time)  ;
        // }
        // std::cout << "indice :" <<  gap << std::endl ;

        boost::shared_ptr<quadruped_walkgen::ActionModelQuadrupedAugmentedTime> model2 =
            boost::dynamic_pointer_cast<quadruped_walkgen::ActionModelQuadrupedAugmentedTime>(running_models[k + gap]);

        // //Update model :
        // tmp = l_feet.array() ;
        // if (int(gait.block(j,1,1,4).sum()) == 4 and gap == 1) {
        //   model2->set_last_position_weights(Eigen::Matrix<double,8,1>::Constant(1)) ;
        // }
        // else{
        //   model2->set_last_position_weights(Eigen::Matrix<double,8,1>::Zero()) ;

        // }
        model2->update_model(Eigen::Map<Eigen::Matrix<double, 3, 4> >(tmp.data(), 3, 4),
                             Eigen::Map<Eigen::Matrix<double, 3, 4> >(tmp.data(), 3, 4),
                             Eigen::Map<Eigen::Matrix<double, 12, 1> >(xref.block(0, k + 1, 12, 1).data(), 12, 1),
                             Eigen::Map<Eigen::Matrix<double, 4, 1> >(gait.block(j, 1, 1, 4).data(), 4, 1));
        us.push_back(u0);
      }
    }
    k_cum += int(gait(j, 0));
  }

  std::cout << "term before" << std::endl;

  boost::shared_ptr<quadruped_walkgen::ActionModelQuadrupedAugmentedTime> terminal_model_2 =
      boost::dynamic_pointer_cast<quadruped_walkgen::ActionModelQuadrupedAugmentedTime>(terminal_model);

  std::cout << "term after" << std::endl;

  tmp = l_feet.array();
  Eigen::Array<double, 1, 4> gait_tmp = Eigen::Array<double, 1, 4>::Zero();
  gait_tmp = gait.block(max_index - 1, 1, 1, 4).array();

  terminal_model_2->update_model(Eigen::Map<Eigen::Matrix<double, 3, 4> >(tmp.data(), 3, 4),
                                 Eigen::Map<Eigen::Matrix<double, 3, 4> >(tmp.data(), 3, 4),
                                 Eigen::Map<Eigen::Matrix<double, 12, 1> >(xref.block(0, 16, 12, 1).data(), 12, 1),
                                 Eigen::Map<Eigen::Matrix<double, 4, 1> >(gait_tmp.data(), 4, 1));
  terminal_model_2->set_force_weights(Eigen::Matrix<double, 12, 1>::Zero());
  terminal_model_2->set_friction_weight(0);
  // terminal_model_2->set_last_position_weights(Eigen::Matrix<double,8,1>::Zero()) ;

  ////////////////////////////////////
  // Code gen
  // ////////////////////////////////////

  // typedef CppAD::AD<CppAD::cg::CG<double> > ADScalar;

  // // Code generation of the running an terminal models
  // boost::shared_ptr<crocoddyl::ActionModelAbstractTpl<ADScalar> > ad_runningModel, ad_terminalModel;
  // crocoddyl::benchmark::build_arm_action_models(ad_runningModel, ad_terminalModel);
  // boost::shared_ptr<crocoddyl::ActionModelAbstract> cg_runningModel =
  //     boost::make_shared<crocoddyl::ActionModelCodeGen>(ad_runningModel, runningModel,
  //     "arm_manipulation_running_cg");
  // boost::shared_ptr<crocoddyl::ActionModelAbstract> cg_terminalModel =
  //     boost::make_shared<crocoddyl::ActionModelCodeGen>(ad_terminalModel, terminalModel,
  //                                                       "arm_manipulation_terminal_cg");

  // for (int k = 0 ; k < running_models.size() ; j++ ) {

  // }
  // std::cout << "----------------------------------------" << std::endl ;

  // std::cout << running_models.size() << std::endl ;
  // for (int j = 0 ; j < running_models.size() ; j++ ) {
  //   std::cout<<j<<std::endl ;

  //   data = running_models[j]->createModel() ;
  //   if (running_models[j]->nu == 4){
  //     boost::shared_ptr<crocoddyl::ActionDataAbstract> data
  //                    = boost::make_shared<quadruped_walkgen::ActionDataQuadrupedStepTime>() ;
  //     data = running_models[j]->createModel() ;
  //     running_models[j]->calc(data,x0,u0_step)
  //   }
  //   if (running_models[j]->nu == 1){
  //     boost::shared_ptr<crocoddyl::ActionDataAbstract> data
  //                    = boost::make_shared<quadruped_walkgen::ActionDataQuadrupedTime>() ;
  //     data = running_models[j]->createModel() ;
  //     running_models[j]->calc(data,x0,u0_time)
  //   }
  //    if (running_models[j]->nu == 12){
  //     boost::shared_ptr<crocoddyl::ActionDataAbstract> data
  //                    = boost::make_shared<quadruped_walkgen::ActionDataQuadrupedAugmentedTime>() ;
  //     data = running_models[j]->createModel() ;
  //     running_models[j]->calc(data,x0,u0)
  //   }
  // }

  boost::shared_ptr<crocoddyl::ShootingProblem> problem =
      boost::make_shared<crocoddyl::ShootingProblem>(x0, running_models, terminal_model);
  crocoddyl::SolverDDP ddp(problem);

  std::cout << "probel ok" << std::endl;

  std::vector<Eigen::VectorXd> xs(running_models.size() + 1, x0);
  // Eigen::Matrix<double,12,1> u0 ;
  // u0 << 1,0.2,8, 1,1,8, -1,1,8, -1,-1,8;
  // std::vector<Eigen::VectorXd> us(int(N), u0);

  Eigen::ArrayXd duration(T);

  // Solving the optimal control problem
  for (unsigned int i = 0; i < T; ++i) {
    crocoddyl::Timer timer;
    ddp.solve(xs, us, MAXITER);
    duration[i] = timer.get_duration();
  }

  double avrg_duration = duration.sum() / T;
  double min_duration = duration.minCoeff();
  double max_duration = duration.maxCoeff();
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