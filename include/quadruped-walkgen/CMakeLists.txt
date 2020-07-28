SET(${PROJECT_NAME}_HEADERS
  quadruped.hpp
  quadruped.hxx
  quadruped_nl.hpp
  quadruped_nl.hxx
  quadruped_augmented.hpp
  quadruped_augmented.hxx
  quadruped_step.hpp
  quadruped_step.hxx
  )

INSTALL(FILES ${${PROJECT_NAME}_HEADERS} DESTINATION include/${PROJECT_NAME})
