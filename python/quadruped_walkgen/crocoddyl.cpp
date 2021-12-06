///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (C) 2018-2020, LAAS-CNRS, University of Edinburgh, INRIA
// Copyright note valid unless otherwise stated in individual files.
// All rights reserved.
///////////////////////////////////////////////////////////////////////////////

#include "fwd.hpp"
#include "crocoddyl/core/utils/version.hpp"
#include "vector-converter.hpp"

namespace quadruped_walkgen {
namespace python {

namespace bp = boost::python;

BOOST_PYTHON_MODULE(quadruped_walkgen_pywrap) {
  bp::scope().attr("__version__") = crocoddyl::printVersion();

  eigenpy::enableEigenPy();

  typedef double Scalar;
  typedef Eigen::Matrix<Scalar, 6, 1> Vector6;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 3> MatrixX3;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VectorX;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixX;

  eigenpy::enableEigenPySpecific<Vector6>();
  eigenpy::enableEigenPySpecific<MatrixX3>();

  // Register converters between std::vector and Python list
  // TODO(cmastalli): figure out how to convert std::vector<double> to Python list
  // bp::to_python_converter<std::vector<double, std::allocator<double> >, vector_to_list<double, false>, true>();
  // bp::to_python_converter<std::vector<VectorX, std::allocator<VectorX> >, vector_to_list<VectorX, false>, true>();
  // bp::to_python_converter<std::vector<MatrixX, std::allocator<MatrixX> >, vector_to_list<MatrixX, false>, true>();
  list_to_vector()
      .from_python<std::vector<double, std::allocator<double> > >()
      .from_python<std::vector<VectorX, std::allocator<VectorX> > >()
      .from_python<std::vector<MatrixX, std::allocator<MatrixX> > >();

  exposeCore();
}

}  // namespace python
}  // namespace quadruped_walkgen
