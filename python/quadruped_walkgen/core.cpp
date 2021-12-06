///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (C) 2018-2020, University of Edinburgh
// Copyright note valid unless otherwise stated in individual files.
// All rights reserved.
///////////////////////////////////////////////////////////////////////////////

#include "core.hpp"

namespace quadruped_walkgen {
namespace python {

void exposeCore() {
  exposeActionAbstract();
  exposeActionQuadruped();
  exposeActionQuadrupedNonLinear();
  exposeActionQuadrupedAugmented();
  exposeActionQuadrupedStep();
  exposeActionQuadrupedAugmentedTime();
  exposeActionQuadrupedStepTime();
  exposeActionQuadrupedTime();
  exposeActionQuadrupedStepPeriod();
}

}  // namespace python
}  // namespace quadruped_walkgen