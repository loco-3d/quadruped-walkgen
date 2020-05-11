# quadruped-walkgen

The objective of this project is the writing of a crocoddyl class in C ++ to use it in a simulation of the quadruped. It is based on a simplified dynamical model of the quadruped and used to solve a MPC problem. 

For now, just trying to create a python binding with an already existing class : unicycle. 
scr : quadruped.cpp  -> same class than unicycle
include : quadruped.hpp -> same class than unicycle

binding : To use bp:bases <> , I need to first create the binding with the ActionModelAbstract because quadruped (unicycle for now) inherits of it. 
