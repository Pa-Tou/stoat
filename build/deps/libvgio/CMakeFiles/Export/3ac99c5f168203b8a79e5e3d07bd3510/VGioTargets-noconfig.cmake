#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "VGio::VGio" for configuration ""
set_property(TARGET VGio::VGio APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(VGio::VGio PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libvgio.so"
  IMPORTED_SONAME_NOCONFIG "libvgio.so"
  )

list(APPEND _cmake_import_check_targets VGio::VGio )
list(APPEND _cmake_import_check_files_for_VGio::VGio "${_IMPORT_PREFIX}/lib/libvgio.so" )

# Import target "VGio::VGio_static" for configuration ""
set_property(TARGET VGio::VGio_static APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(VGio::VGio_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "CXX"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libvgio.a"
  )

list(APPEND _cmake_import_check_targets VGio::VGio_static )
list(APPEND _cmake_import_check_files_for_VGio::VGio_static "${_IMPORT_PREFIX}/lib/libvgio.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
