#include <quadruped-walkgen/gait_problem.hpp>
#include <quadruped-walkgen/quadruped.hpp>

namespace quadruped_walkgen {

gaitProblem::gaitProblem(){
    dt_ = 0.02 ; 
    mu_ = 0.8 ;
    T_mpc_ = 0.32 ; 
    max_iter_ = 15 ; 
}

}