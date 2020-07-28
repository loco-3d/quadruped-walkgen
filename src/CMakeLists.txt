SET(${PROJECT_NAME}_SOURCES
  quadruped.cpp
  quadruped_nl.cpp
  quadruped_augmented.cpp
  quadruped_step.cpp
  )

ADD_LIBRARY(${PROJECT_NAME} SHARED ${${PROJECT_NAME}_SOURCES})
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${Boost_LIBRARIES} crocoddyl::crocoddyl)

INSTALL(TARGETS ${PROJECT_NAME} DESTINATION lib)
