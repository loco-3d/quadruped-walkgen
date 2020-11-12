SET(${PROJECT_NAME}_BENCHMARK
  quadruped
  quadruped-non-linear
  quadruped-planner
  quadruped-planner-period
  )

FOREACH(BENCHMARK_NAME ${${PROJECT_NAME}_BENCHMARK})
  ADD_EXECUTABLE(${BENCHMARK_NAME} ${BENCHMARK_NAME}.cpp)
  TARGET_LINK_LIBRARIES(${BENCHMARK_NAME} ${PROJECT_NAME} example-robot-data::example-robot-data)
  ADD_CUSTOM_TARGET("benchmarks-cpp-${BENCHMARK_NAME}"
    ${BENCHMARK_NAME} \${INPUT})
ENDFOREACH(BENCHMARK_NAME ${${PROJECT_NAME}_BENCHMARK})
