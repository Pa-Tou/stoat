file(REMOVE_RECURSE
  "libhandlegraph.a"
  "libhandlegraph.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang CXX)
  include(CMakeFiles/handlegraph_static.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
