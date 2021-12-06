# quadruped-walkgen

The objective of this project is the writing of a crocoddyl class in C ++ to use it in a simulation of the quadruped. It is based on a simplified dynamical model of the quadruped and used to solve a MPC problem.
The quadruped action model class is mainly based on the unicycle example class which derives from on the ActionModelAbstract class.

To install :
```bash
git clone --recursive https://github.com/loco-3d/quadruped-walkgen.git
mkdir quadruped-walkgen/build
cd quadruped-walkgen/build
cmake .. -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=your_path
make install
```

Python binding :
- add `your_path` to python path
- `import quadruped_walkgen`

To run the benchmark in cpp, from the build directory in release mode:
```bash
make -s benchmarks-cpp-quadruped INPUT="1000 5"
INPUT="nb of trials , maximum iteration for ddp solver"
```

To run the benchmark in python, from benchmark folder :
```bash
python3 quadruped.py
```
