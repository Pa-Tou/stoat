#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "libhandlegraph::handlegraph_shared" for configuration ""
set_property(TARGET libhandlegraph::handlegraph_shared APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(libhandlegraph::handlegraph_shared PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libhandlegraph.so"
  IMPORTED_SONAME_NOCONFIG "libhandlegraph.so"
  )

list(APPEND _cmake_import_check_targets libhandlegraph::handlegraph_shared )
list(APPEND _cmake_import_check_files_for_libhandlegraph::handlegraph_shared "${_IMPORT_PREFIX}/lib/libhandlegraph.so" )

# Import target "libhandlegraph::handlegraph_static" for configuration ""
set_property(TARGET libhandlegraph::handlegraph_static APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(libhandlegraph::handlegraph_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "CXX"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libhandlegraph.a"
  )

list(APPEND _cmake_import_check_targets libhandlegraph::handlegraph_static )
list(APPEND _cmake_import_check_files_for_libhandlegraph::handlegraph_static "${_IMPORT_PREFIX}/lib/libhandlegraph.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
