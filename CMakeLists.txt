cmake_minimum_required(VERSION 3.10)
project(stoat VERSION 0.1.0 LANGUAGES CXX)

# C++ Standard Settings
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_OSX_DEPLOYMENT_TARGET "12.7") # IOS fix

set(deps_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps")

# Put the executables in a bin directory
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# To change OPTIMIZE use : cmake -DOPTIMIZE=OFF ..
option(OPTIMIZE "Build with optimization [default: on]" ON)
message(STATUS "OPTIMIZE is set to: ${OPTIMIZE}")

# This is the same flag used by libbdsg and libvgio to specify which version of libhandlegraph should be used
OPTION(USE_INSTALLED_LIBHANDLEGRAPH "Use the version of libhandlegraph installed on the system, if it exists, instead of the bundled version [default: off]" OFF)

# DEBUG STEP : 
# cmake -DCMAKE_BUILD_TYPE=Debug .
# make -j 8
# cd ../bin
# gdb --args ./stoat ... (command that you want to test)
# (gdb) run

# Optimize option
if (OPTIMIZE)
  message(STATUS "Optimization compilation mode activated")
  # Use all standard-compliant optimizations
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
else()
  # Debug build: no optimization, include debug symbols
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O0")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0")
  set(CMAKE_BUILD_TYPE Debug)
  message(STATUS "Debug compilation mode activated: -g -O0 flags enabled")
endif()

find_package(PkgConfig REQUIRED)
pkg_check_modules(HTSLIB REQUIRED htslib)

include_directories(${HTSLIB_INCLUDE_DIRS})
link_directories(${HTSLIB_LIBRARY_DIRS})

# Find Packages
find_package(Boost REQUIRED)

#Get eigen3
include_directories(${deps_DIR}/eigen)

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    # Detect Homebrew prefix
    set(HOMEBREW_PREFIX "$ENV{HOMEBREW_PREFIX}")
    if ("${HOMEBREW_PREFIX}" STREQUAL "")
        execute_process(COMMAND brew --prefix OUTPUT_VARIABLE HOMEBREW_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE)
    endif()
    if ("${HOMEBREW_PREFIX}" STREQUAL "")
        set(HOMEBREW_PREFIX "/opt/homebrew")  # fallback
    endif()

    # Set OpenMP manually (Apple Clang needs this)
    set(OpenMP_INCLUDES "${HOMEBREW_PREFIX}/opt/libomp/include")
    set(OpenMP_LIBRARIES "${HOMEBREW_PREFIX}/opt/libomp/lib/libomp.dylib")
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${OpenMP_INCLUDES}")
    set(OpenMP_LINK_FLAGS "-L${HOMEBREW_PREFIX}/opt/libomp/lib -lomp")

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_LINK_FLAGS}")
endif()

find_package(OpenMP REQUIRED)

# Produces handlegraph_shared and handlegraph_static targets (as well as handlegraph_objs)
if (USE_INSTALLED_LIBHANDLEGRAPH)
    find_package(libhandlegraph)
endif()
if (${libhandlegraph_FOUND})
    message("STOAT using installed libhandlegraph")
else ()
    message("STOAT using bundled libhandlegraph")
    add_subdirectory(${deps_DIR}/libhandlegraph)
endif ()

# Add libvgio using its cmake build
add_subdirectory(${deps_DIR}/libvgio)

# libbdsg (full build using its cmake config)
# Produces libbdsg target, dynamic only.
# Comes along with its own libhandlegraph, which it depends on.
# Make sure to never build its Python bindings.
set(BUILD_PYTHON_BINDINGS OFF CACHE BOOL "Build Python bindings for libbdsg" FORCE)
set(RUN_DOXYGEN OFF CACHE BOOL "Run libbdsg Doxygen" FORCE)
add_subdirectory(${deps_DIR}/libbdsg)

# Copied from https://github.com/vgteam/vg/wiki/Programming-with-the-vg-API
# SDSL library ("sdsl" target) already exists from libbdsg
# GBWT and GBWTgraph don't have CMake yet so we need to explain how to build it.

# GBWT
foreach (gbwt_source IN ITEMS algorithms.cpp bwtmerge.cpp cached_gbwt.cpp dynamic_gbwt.cpp fast_locate.cpp files.cpp gbwt.cpp internal.cpp metadata.cpp support.cpp test.cpp utils.cpp variants.cpp)
list(APPEND gbwt_libFiles "deps/gbwt/src/${gbwt_source}")
endforeach()
if(MSVC)
file(GLOB gbwt_headerFiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}  "${deps_DIR}/gbwt/include/gbwt/*.h")
endif()
set( gbwt_SRCS ${gbwt_libFiles} ${gbwt_headerFiles})
add_library( gbwt ${gbwt_SRCS} )
target_include_directories( gbwt PUBLIC "${deps_DIR}/gbwt/include")
target_link_libraries( gbwt PUBLIC sdsl OpenMP::OpenMP_CXX )
install(TARGETS gbwt
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    PUBLIC_HEADER DESTINATION include)

# GBWTGraph
foreach (gbwtgraph_source IN ITEMS algorithms.cpp cached_gbwtgraph.cpp gbwtgraph.cpp gbz.cpp gfa.cpp internal.cpp minimizer.cpp path_cover.cpp utils.cpp)
list(APPEND gbwtgraph_libFiles "deps/gbwtgraph/src/${gbwtgraph_source}")
endforeach()
if(MSVC)
file(GLOB gbwtgraph_headerFiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}  "${deps_DIR}/gbwtgraph/include/gbwt/*.h")
endif()
set( gbwtgraph_SRCS ${gbwtgraph_libFiles} ${gbwtgraph_headerFiles})
add_library( gbwtgraph ${gbwtgraph_SRCS} )
target_include_directories( gbwtgraph PUBLIC "${deps_DIR}/gbwtgraph/include")
target_link_libraries( gbwtgraph PUBLIC gbwt handlegraph_shared )
install(TARGETS gbwtgraph
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    PUBLIC_HEADER DESTINATION include)

add_library(${PROJECT_NAME}_objs OBJECT
    src/utils.cpp
    src/arg_parser.cpp
    src/post_processing.cpp
    src/gaf_creator.cpp
    src/snarl_data_t.cpp
    src/matrix.cpp
    src/log.cpp
    src/binary_table.cpp
    src/quantitative_table.cpp
    src/stats_test.cpp
    src/writer.cpp
    src/snarl_analyzer.cpp
    src/graph_path_association_finder.cpp
    src/partitioner.cpp
    
    src/io/register_loader_saver_gbz.cpp
    src/io/register_loader_saver_gbzgraph.cpp
    src/io/register_loader_saver_hash_graph.cpp
    src/io/register_loader_saver_packed_graph.cpp
    src/io/register_io.cpp

    src/subcommand/vcf.cpp
    src/subcommand/graph.cpp
    src/subcommand/bh_correct.cpp
)

# Include Directories
target_include_directories(${PROJECT_NAME}_objs PRIVATE
${Boost_INCLUDE_DIRS}
${HTSLIB_INCLUDE_DIRS}
)

# Link Libraries
target_link_libraries(${PROJECT_NAME}_objs
    Boost::boost
    ${HTSLIB_LIBRARIES}
    libbdsg
    gbwtgraph
    VGio::vgio
)

# Add Executable
add_executable(${PROJECT_NAME} src/main.cpp)

# Link the code
target_link_libraries(${PROJECT_NAME} ${PROJECT_NAME}_objs)

target_compile_options(${PROJECT_NAME} PRIVATE -w) # fix warning

# Enable testing
enable_testing()

# Add test executable
add_executable(unit_tests
    tests/unittest/main_test.cpp
    tests/unittest/binary_test_unit.cpp
    tests/unittest/gaf_creator_unit.cpp
    tests/unittest/snarl_data_t_unit.cpp
    tests/unittest/matrix_unit.cpp
    tests/unittest/quantitative_test_unit.cpp
    tests/unittest/path_partitioner_unit.cpp
    tests/unittest/arg_parsing_unit.cpp
    tests/unittest/post_processing_unit.cpp
    tests/unittest/snarl_analyzer_unit.cpp
    tests/unittest/utils_unit.cpp
)

target_compile_options(unit_tests PRIVATE -w) # fix warning

target_link_libraries(unit_tests PRIVATE 
    ${PROJECT_NAME}_objs)

include_directories(${deps_DIR}/Catch2/single_include/catch2 )

add_test(NAME UnitTests COMMAND unit_tests)

# Add test executable system tests
add_executable(vcf_simu_test tests/system/compare_files_utils.cpp tests/system/vcf_simu_test.cpp)

# add src target lib + Catch2 lib 
target_link_libraries(vcf_simu_test PRIVATE 
    Boost::boost
    ${HTSLIB_LIBRARIES})

add_test(NAME VcfSimuTest COMMAND vcf_simu_test)

add_executable(graph_simu_test tests/system/compare_files_utils.cpp tests/system/graph_simu_test.cpp)

# add src target lib + Catch2 lib 
target_link_libraries(graph_simu_test PRIVATE 
    Boost::boost
    ${HTSLIB_LIBRARIES})

# Register with CTest
add_test(NAME GraphSimuTest COMMAND graph_simu_test)
