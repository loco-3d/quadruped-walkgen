file(REMOVE_RECURSE
  "quadruped_walkgen/__init__.pyc"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/compile_pyc.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
