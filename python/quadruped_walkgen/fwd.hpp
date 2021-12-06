///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (C) 2018-2020, University of Edinburgh
// Copyright note valid unless otherwise stated in individual files.
// All rights reserved.
///////////////////////////////////////////////////////////////////////////////

#ifndef BINDINGS_PYTHON_QUADRUPED_WALKGEN_FWD_HPP_
#define BINDINGS_PYTHON_QUADRUPED_WALKGEN_FWD_HPP_

#include <eigenpy/eigenpy.hpp>
#include <boost/python.hpp>
#include <boost/python/enum.hpp>

namespace quadruped_walkgen {
namespace python {

namespace bp = boost::python;

void exposeCore();

}  // namespace python
}  // namespace quadruped_walkgen

#endif  // BINDINGS_PYTHON_QUADRUPED_WALKGEN_FWD_HPP_