set(${PROJECT_NAME}_BENCHMARK quadruped quadruped-non-linear quadruped-planner
                              quadruped-planner-period)

foreach(BENCHMARK_NAME ${${PROJECT_NAME}_BENCHMARK})
  add_executable(${BENCHMARK_NAME} ${BENCHMARK_NAME}.cpp)
  target_link_libraries(${BENCHMARK_NAME} ${PROJECT_NAME}
                        example-robot-data::example-robot-data)
  add_custom_target("benchmarks-cpp-${BENCHMARK_NAME}" ${BENCHMARK_NAME}
                                                       \${INPUT})
endforeach(BENCHMARK_NAME ${${PROJECT_NAME}_BENCHMARK})
