# quadruped-walkgen

The objective of this project is the writing of a crocoddyl class in C ++ to use it in a simulation of the quadruped. It is based on a simplified dynamical model of the quadruped and used to solve a MPC problem. 

To run the exemple in pyhton, from the directory exemple : 
python3 exemple_simple.py

To run the benchmark in cpp, from the build directory in release mode: 
make -s benchmarks-cpp-quadruped INPUT="1000 5"
INPUT="nb of trials , maximum iteration for ddp solver"

To run the benchmark in python, from benchmark folder : 
python3 quadruped.py

The results obtained are summarized in the results.pdf

