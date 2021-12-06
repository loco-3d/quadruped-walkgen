///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (C) 2018-2020, LAAS-CNRS, University of Edinburgh
// Copyright note valid unless otherwise stated in individual files.
// All rights reserved.
///////////////////////////////////////////////////////////////////////////////

#ifndef BINDINGS_PYTHON_QUADRUPED_WALKGEN_CORE_HPP_
#define BINDINGS_PYTHON_QUADRUPED_WALKGEN_CORE_HPP_

#include "fwd.hpp"

namespace quadruped_walkgen {
namespace python {

namespace bp = boost::python;

void exposeActionAbstract();
void exposeActionQuadruped();
void exposeActionQuadrupedNonLinear();
void exposeActionQuadrupedAugmented();
void exposeActionQuadrupedStep();
void exposeActionQuadrupedAugmentedTime();
void exposeActionQuadrupedStepTime();
void exposeActionQuadrupedTime();
void exposeActionQuadrupedStepPeriod();

void exposeCore();

}  // namespace python
}  // namespace quadruped_walkgen

#endif  // BINDINGS_PYTHON_QUADRUPED_WALKGEN_CORE_HPP_
