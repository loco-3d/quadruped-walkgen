set(${PROJECT_NAME}_PYTHON_BINDINGS_HEADERS
    ${PYTHON_DIR}/core.hpp ${PYTHON_DIR}/action-base.hpp ${PYTHON_DIR}/fwd.hpp
    ${PYTHON_DIR}/vector-converter.hpp)

set(${PROJECT_NAME}_PYTHON_BINDINGS_SOURCES
    ${PYTHON_DIR}/crocoddyl.cpp
    ${PYTHON_DIR}/core.cpp
    ${PYTHON_DIR}/action-base.cpp
    ${PYTHON_DIR}/quadruped.cpp
    ${PYTHON_DIR}/quadruped_nl.cpp
    ${PYTHON_DIR}/quadruped_augmented.cpp
    ${PYTHON_DIR}/quadruped_step.cpp
    ${PYTHON_DIR}/quadruped_augmented_time.cpp
    ${PYTHON_DIR}/quadruped_step_time.cpp
    ${PYTHON_DIR}/quadruped_step_period.cpp
    ${PYTHON_DIR}/quadruped_time.cpp)
add_library(
  ${PYTHON_DIR}_pywrap SHARED ${${PROJECT_NAME}_PYTHON_BINDINGS_SOURCES}
                              ${${PROJECT_NAME}_PYTHON_BINDINGS_HEADERS})
set_target_properties(
  ${PYTHON_DIR}_pywrap
  PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PYTHON_DIR}
             PREFIX ""
             SUFFIX ${PYTHON_EXT_SUFFIX})
target_link_libraries(${PYTHON_DIR}_pywrap PUBLIC ${PROJECT_NAME}
                                                  eigenpy::eigenpy)
target_link_boost_python(${PYTHON_DIR}_pywrap PRIVATE)

target_compile_options(${PYTHON_DIR}_pywrap PRIVATE "-Wno-conversion")

set(${PROJECT_NAME}_PYTHON_BINDINGS_FILES __init__.py)

foreach(python ${${PROJECT_NAME}_PYTHON_BINDINGS_FILES})
  python_install_on_site(${PYTHON_DIR} ${python})
endforeach(python ${${PROJECT_NAME}_PYTHON_BINDINGS_FILES})

install(TARGETS ${PYTHON_DIR}_pywrap
        DESTINATION ${PYTHON_SITELIB}/${PYTHON_DIR})
