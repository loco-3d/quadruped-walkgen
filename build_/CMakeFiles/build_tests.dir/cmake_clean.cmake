file(REMOVE_RECURSE
  "doc/quadruped-walkgen.doxytag"
  "doc/doxygen.log"
  "doc/doxygen-html"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/build_tests.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
