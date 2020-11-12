SET(${PROJECT_NAME}_HEADERS
  quadruped.hpp
  quadruped.hxx
  quadruped_nl.hpp
  quadruped_nl.hxx
  quadruped_augmented.hpp
  quadruped_augmented.hxx
  quadruped_step.hpp
  quadruped_step.hxx
  quadruped_augmented_time.hpp
  quadruped_augmented_time.hxx
  quadruped_step_time.hpp
  quadruped_step_time.hxx
  quadruped_step_period.hpp
  quadruped_step_period.hxx
  quadruped_time.hpp
  quadruped_time.hxx
  )

INSTALL(FILES ${${PROJECT_NAME}_HEADERS} DESTINATION include/${PROJECT_NAME})
