#ifndef __quadruped_walkgen_gait_problem_hpp__
#define __quadruped_walkgen_gait_problem_hpp__
#include <stdexcept>

#include "crocoddyl/core/fwd.hpp"
#include "crocoddyl/core/action-base.hpp"

#include <crocoddyl/core/optctrl/shooting.hpp>
#include <crocoddyl/core/solver-base.hpp>

namespace quadruped_walkgen {
class gaitProblem {
  
 public:

  //typedef crocoddyl::ActionModelAbstract ActionModelAbstract;
  //boost::shared_ptr<crocoddyl::SolverDDP> ddp;
  //boost::shared_ptr<crocoddyl::ShootingProblem> shooting_problem;
  
  gaitProblem();
  ~gaitProblem();

  //void createProblem() ; 
  //void updateProblem(const Eigen::Matrix<double,6,5>& f_steps ,
    //                  const Eigen::Matrix<double,12,17>& x_ref ,
      //                const Eigen::Matrix<double,12,1>& x_0 ) ; 

  

  
 
 //protected:
  //std::vector<boost::shared_ptr<ActionModelAbstract> > running_models_;
  //boost::shared_ptr<ActionModelAbstract> terminal_model_;

 private:
  double dt_ ; 
  double T_mpc_ ;
  int max_iter_ ; 
  double mu_ ; 


 
};

}



#endif
