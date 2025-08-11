file(REMOVE_RECURSE
  "libhandlegraph.pdb"
  "libhandlegraph.so"
)

# Per-language clean rules from dependency scanning.
foreach(lang CXX)
  include(CMakeFiles/handlegraph_shared.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
