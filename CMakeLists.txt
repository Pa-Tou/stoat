cmake_minimum_required(VERSION 3.10)
project(stoat VERSION 0.0.4 LANGUAGES CXX)

# Needed for the STATIC_BUILD libbdsg workaround further down.
include(ExternalProject)

add_compile_definitions(STOAT_VERSION="v${PROJECT_VERSION}")
message(STATUS "Stoat version: ${PROJECT_VERSION}")

# C++ Standard Settings
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_OSX_DEPLOYMENT_TARGET "12.7") # MACOS fix

include_directories(BEFORE
    ${PROJECT_SOURCE_DIR}/deps/libbdsg/bdsg/deps/sdsl-lite/external/googletest/googletest/include
)

set(deps_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps")

# Put the executables in a bin directory
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# To change to OPTIMIZE=OFF use : cmake -DOPTIMIZE=OFF ..
option(OPTIMIZE "Build with optimization [default: on]" ON)
message(STATUS "OPTIMIZE is set to: ${OPTIMIZE}")

# NOTE: renamed from STOAT_STATIC -> STATIC_BUILD. Dockerfile.static already
# invokes `cmake .. -DSTATIC_BUILD=ON`; with the old name that flag was
# silently ignored (unknown -D args don't error), so the "static" build was
# actually a normal dynamic build that happened to succeed.
option(STATIC_BUILD "Build stoat as a fully static executable [Linux only]" OFF)
message(STATUS "STATIC_BUILD is set to: ${STATIC_BUILD}")

# This is the same flag used by libbdsg and libvgio to specify which version of libhandlegraph should be used
OPTION(USE_INSTALLED_LIBHANDLEGRAPH "Use the version of libhandlegraph installed on the system, if it exists, instead of the bundled version [default: off]" OFF)

# DEBUG STEP : 
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

# ---------------------------------------------------------------------------
# Static-build scaffolding. Everything here is set BEFORE any add_subdirectory()
# call, so every dependency configured below sees the same, consistent
# settings from the moment it is first processed.
# ---------------------------------------------------------------------------
if (STATIC_BUILD)
    if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        message(FATAL_ERROR "STATIC_BUILD is not supported on Darwin")
    endif()

    # Make any nested add_library() call that doesn't hardcode a type
    # (STATIC/SHARED) default to static. CMake already defaults to static
    # when this is unset, but some vendored projects set their own default
    # to ON, so force it explicitly rather than relying on that.
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build static libraries by default" FORCE)

    # Prefer .a over .so when find_package()/find_library() resolve system
    # libraries (Boost, Threads, OpenMP, bz2, lzma, ...).
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")

    # Keep every compiled object PIC-consistent across the static archives.
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)

    # Fully static Linux executable.
    set(CMAKE_EXE_LINKER_FLAGS
        "${CMAKE_EXE_LINKER_FLAGS} -static -static-libgcc -static-libstdc++ -pthread"
    )

    # A "-static" link gets none of a dependency's *transitive* shared-lib
    # requirements for free (there's no dynamic loader resolving them at
    # runtime), so every system library htslib's static archive needs has to
    # be listed explicitly. pkg-config's own "--static" output (queried
    # automatically below via HTSLIB_STATIC_LIBRARIES) should cover most of
    # this; these are kept as an explicit safety net.
    find_library(BZ2_LIB  NAMES bz2  REQUIRED)
    find_library(LZMA_LIB NAMES lzma REQUIRED)
    set(STOAT_STATIC_SYSTEM_LIBS ${BZ2_LIB} ${LZMA_LIB} z dl m)

    message(STATUS "Configured for FULL STATIC Linux build")
endif()

find_package(PkgConfig REQUIRED)
pkg_check_modules(HTSLIB REQUIRED htslib)
# pkg_check_modules() always additionally computes HTSLIB_STATIC_LIBRARIES /
# HTSLIB_STATIC_LDFLAGS / HTSLIB_STATIC_LIBRARY_DIRS from `pkg-config --static`
# (this is a standard FindPkgConfig feature, no extra call needed) -- use
# those instead of the plain ones when linking a static executable.
if (STATIC_BUILD)
    set(STOAT_HTSLIB_LIBRARIES ${HTSLIB_STATIC_LIBRARIES})
    link_directories(${HTSLIB_STATIC_LIBRARY_DIRS})
else()
    set(STOAT_HTSLIB_LIBRARIES ${HTSLIB_LIBRARIES})
    link_directories(${HTSLIB_LIBRARY_DIRS})
endif()

include_directories(${HTSLIB_INCLUDE_DIRS})

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

# handlegraph_shared / vgio (shared) / libbdsg are libhandlegraph's, libvgio's
# and libbdsg's default (dynamic) targets. libhandlegraph and libvgio always
# additionally produce handlegraph_static / vgio_static unconditionally (no
# flag required), which is why only STOAT_BDSG_TARGET needs real work below.
set(STOAT_HANDLEGRAPH_TARGET handlegraph_shared)
set(STOAT_VGIO_TARGET VGio::vgio)
set(STOAT_BDSG_TARGET libbdsg)

if (STATIC_BUILD)
    set(STOAT_HANDLEGRAPH_TARGET handlegraph_static)
    set(STOAT_VGIO_TARGET vgio_static)
    # STOAT_BDSG_TARGET is set to bdsg::static further down, once that
    # target actually exists.
endif()

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

# Add libvgio using its cmake build. It always builds both its shared
# (VGio::vgio) and static (vgio_static) targets, so no extra flag is needed
# here regardless of STATIC_BUILD.
add_subdirectory(${deps_DIR}/libvgio EXCLUDE_FROM_ALL)

# libbdsg
# Make sure to never build its Python bindings.
set(BUILD_PYTHON_BINDINGS OFF CACHE BOOL "Build Python bindings for libbdsg" FORCE)
set(RUN_DOXYGEN OFF CACHE BOOL "Run libbdsg Doxygen" FORCE)

if (STATIC_BUILD)
    # ---------------------------------------------------------------------
    # libbdsg's own CMakeLists.txt exposes a BUILD_STATIC option, but it is
    # not enough on its own: libbdsg vendors sdsl-lite, which is built by a
    # foreign, non-CMake `install.sh` script wrapped by libbdsg's own build
    # rules. That wrapping does not forward this project's CXXFLAGS / PIC /
    # portability settings into the script, and that forwarding can only be
    # added by editing libbdsg's (or sdsl-lite's) own CMakeLists.txt/build
    # files -- which is out of scope here. So, for the static case only,
    # libbdsg and its vendored dependencies are built directly as isolated
    # ExternalProject_Add steps, driven by their own native build systems
    # (install.sh / CMake / make) but with every flag they need supplied
    # explicitly from this file, and the result is exposed as a normal
    # IMPORTED CMake target (bdsg::static). This mirrors exactly how
    # libbdsg is built upstream for a static `vg`, and touches nothing
    # under deps/.
    # ---------------------------------------------------------------------

    set(LIBBDSG_SRC    "${deps_DIR}/libbdsg")
    set(BDSG           "${LIBBDSG_SRC}/bdsg")

    set(DEPS_INCLUDE    "${BDSG}/include")
    set(DEPS_LIBBDSG    "${BDSG}/deps")
    set(DEPS_LIB        "${deps_DIR}/lib")
    set(HANDLEGRAPH_LIB "${DEPS_LIBBDSG}/lib")
    set(SDSL_LIB        "${DEPS_LIBBDSG}/lib")

    file(MAKE_DIRECTORY ${DEPS_LIB})
    file(MAKE_DIRECTORY ${DEPS_INCLUDE})
    include_directories(${DEPS_INCLUDE})
    link_directories(${DEPS_LIB})

    set(SDSL_SRC           "${DEPS_LIBBDSG}/sdsl-lite")
    set(LIBHANDLEGRAPH_SRC "${DEPS_LIBBDSG}/libhandlegraph")
    set(SPARSEPP_SRC       "${DEPS_LIBBDSG}/sparsepp")
    set(DYNAMIC_SRC        "${DEPS_LIBBDSG}/DYNAMIC")
    set(MIO_SRC            "${DEPS_LIBBDSG}/mio")
    set(BBHASH_SRC         "${DEPS_LIBBDSG}/BBHash")

    # Include directories needed to compile against libbdsg + everything it
    # vendors (used by gbwt/gbwtgraph and stoat's own sources below).
    set(STATIC_DEPS_INCLUDE_DIRS
        "${DEPS_INCLUDE}"
        "${DEPS_LIBBDSG}/include"
        "${LIBHANDLEGRAPH_SRC}/src/include"
        "${deps_DIR}/libvgio/include"
        "${DEPS_LIBBDSG}/sparsepp"
        "${DEPS_LIBBDSG}/DYNAMIC/include"
        "${DEPS_LIBBDSG}/mio/include"
        "${DEPS_LIBBDSG}/BBHash"
        "${deps_DIR}/gbwt/include"
        "${deps_DIR}/gbwtgraph/include"
    )

    set(COMMON_CXX_FLAGS "-fPIC -std=c++17 ${CMAKE_CXX_FLAGS}")

    # 1. sdsl-lite (foreign, install.sh-driven build)
    ExternalProject_Add(sdsl_lite
        SOURCE_DIR ${SDSL_SRC}
        CONFIGURE_COMMAND ""
        BUILD_IN_SOURCE 1
        BUILD_COMMAND ${CMAKE_COMMAND} -E env
            BUILD_PORTABLE=1
            CXXFLAGS=${COMMON_CXX_FLAGS}
            ./install.sh ${DEPS_LIBBDSG}
        INSTALL_COMMAND ""
        BUILD_BYPRODUCTS
            ${SDSL_LIB}/libsdsl.a
            ${SDSL_LIB}/libdivsufsort.a
            ${SDSL_LIB}/libdivsufsort64.a
    )

    # 2. libbdsg's own vendored libhandlegraph copy
    ExternalProject_Add(libhandlegraph_bdsg
        SOURCE_DIR ${LIBHANDLEGRAPH_SRC}
        CMAKE_ARGS
            -DCMAKE_INSTALL_PREFIX=${DEPS_LIBBDSG}
            -DCMAKE_INSTALL_LIBDIR=lib
            -DCMAKE_CXX_FLAGS=${COMMON_CXX_FLAGS}
            -DCMAKE_BUILD_TYPE=Release
        BUILD_BYPRODUCTS ${HANDLEGRAPH_LIB}/libhandlegraph.a
    )

    # 3. sparsepp (header-only)
    ExternalProject_Add(sparsepp
        SOURCE_DIR ${SPARSEPP_SRC}
        DOWNLOAD_COMMAND ""
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${SPARSEPP_SRC}/sparsepp ${DEPS_INCLUDE}/sparsepp
    )

    # 4. DYNAMIC
    ExternalProject_Add(dynamic
        SOURCE_DIR ${DYNAMIC_SRC}
        CMAKE_ARGS
            -DCMAKE_CXX_FLAGS=${COMMON_CXX_FLAGS}
        INSTALL_COMMAND
            ${CMAKE_COMMAND} -E copy_directory
                ${DYNAMIC_SRC}/deps/hopscotch_map/include
                ${DEPS_INCLUDE}
            COMMAND ${CMAKE_COMMAND} -E make_directory
                ${DEPS_INCLUDE}/dynamic
            COMMAND ${CMAKE_COMMAND} -E copy_directory
                ${DYNAMIC_SRC}/include/dynamic
                ${DEPS_INCLUDE}/dynamic
    )

    # 5. mio (header-only)
    ExternalProject_Add(mio
        SOURCE_DIR ${MIO_SRC}
        DOWNLOAD_COMMAND ""
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${MIO_SRC}/include/mio ${DEPS_INCLUDE}/mio
    )

    # 6. BBHash (header-only)
    ExternalProject_Add(bbhash
        SOURCE_DIR ${BBHASH_SRC}
        DOWNLOAD_COMMAND ""
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ${CMAKE_COMMAND} -E copy
            ${BBHASH_SRC}/BooPHF.h
            ${DEPS_INCLUDE}/BooPHF.h
    )

    # 7. libbdsg itself, driven by its own Makefile (not its CMakeLists.txt),
    # now that everything it vendors has been built with the right flags.
    ExternalProject_Add(libbdsg
        SOURCE_DIR ${LIBBDSG_SRC}
        DEPENDS
            sdsl_lite
            libhandlegraph_bdsg
            sparsepp
            dynamic
            mio
            bbhash
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND ""
        BUILD_COMMAND
            ${CMAKE_COMMAND} -E env
                CPLUS_INCLUDE_PATH=${DEPS_INCLUDE}:${LIBHANDLEGRAPH_SRC}/src/include:${SDSL_SRC}/include:${DEPS_INCLUDE}/dynamic
                CXXFLAGS=-I${DEPS_INCLUDE}\ -I${LIBHANDLEGRAPH_SRC}/src/include\ -I${SDSL_SRC}/include\ -I${DEPS_INCLUDE}/dynamic\ -fPIC\ ${COMMON_CXX_FLAGS}
                make clean
            COMMAND ${CMAKE_COMMAND} -E env
                CPLUS_INCLUDE_PATH=${DEPS_INCLUDE}:${LIBHANDLEGRAPH_SRC}/src/include:${SDSL_SRC}/include:${DEPS_INCLUDE}/dynamic
                CXXFLAGS=-I${DEPS_INCLUDE}\ -I${LIBHANDLEGRAPH_SRC}/src/include\ -I${SDSL_SRC}/include\ -I${DEPS_INCLUDE}/dynamic\ -fPIC\ ${COMMON_CXX_FLAGS}
                make
        INSTALL_COMMAND
            ${CMAKE_COMMAND} -E copy
                ${LIBBDSG_SRC}/lib/libbdsg.a
                ${DEPS_LIB}/libbdsg.a
            COMMAND ${CMAKE_COMMAND} -E copy_directory
                ${LIBBDSG_SRC}/bdsg/include
                ${DEPS_INCLUDE}
        BUILD_BYPRODUCTS ${DEPS_LIB}/libbdsg.a
    )

    add_library(bdsg_static STATIC IMPORTED GLOBAL)
    set_target_properties(bdsg_static PROPERTIES
        IMPORTED_LOCATION ${DEPS_LIB}/libbdsg.a
        INTERFACE_INCLUDE_DIRECTORIES ${DEPS_INCLUDE}
    )
    add_dependencies(bdsg_static libbdsg)
    add_library(bdsg::static ALIAS bdsg_static)

    set(STOAT_BDSG_TARGET bdsg::static)

else()
    # Unchanged: full dynamic build using libbdsg's own cmake config.
    # Produces libbdsg target (dynamic) and, as a side effect, the "sdsl"
    # target used by gbwt/gbwtgraph below.
    add_subdirectory(${deps_DIR}/libbdsg EXCLUDE_FROM_ALL)
endif()

# Copied from https://github.com/vgteam/vg/wiki/Programming-with-the-vg-API
# GBWT and GBWTgraph don't have CMake yet so we need to explain how to build it.

# GBWT
foreach (gbwt_source IN ITEMS algorithms.cpp bwtmerge.cpp cached_gbwt.cpp dynamic_gbwt.cpp fast_locate.cpp files.cpp gbwt.cpp internal.cpp metadata.cpp support.cpp test.cpp utils.cpp variants.cpp)
list(APPEND gbwt_libFiles "deps/gbwt/src/${gbwt_source}")
endforeach()
if(MSVC)
file(GLOB gbwt_headerFiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}  "${deps_DIR}/gbwt/include/gbwt/*.h")
endif()
set( gbwt_SRCS ${gbwt_libFiles} ${gbwt_headerFiles})

if (STATIC_BUILD)
    add_library( gbwt STATIC ${gbwt_SRCS} )
    target_include_directories( gbwt PUBLIC "${deps_DIR}/gbwt/include" ${STATIC_DEPS_INCLUDE_DIRS})
    target_link_libraries( gbwt PUBLIC "${SDSL_LIB}/libsdsl.a" OpenMP::OpenMP_CXX )
    add_dependencies( gbwt sdsl_lite )
else()
    add_library( gbwt ${gbwt_SRCS} )
    target_include_directories( gbwt PUBLIC "${deps_DIR}/gbwt/include")
    target_link_libraries( gbwt PUBLIC sdsl OpenMP::OpenMP_CXX )
endif()

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

if (STATIC_BUILD)
    add_library( gbwtgraph STATIC ${gbwtgraph_SRCS} )
    target_include_directories( gbwtgraph PUBLIC "${deps_DIR}/gbwtgraph/include" "${deps_DIR}/gbwt/include" ${STATIC_DEPS_INCLUDE_DIRS})
    # Deliberately link libbdsg's own vendored libhandlegraph.a here (not
    # the top-level handlegraph_static target) so that only one compiled
    # copy of libhandlegraph ends up in the final static archive set.
    target_link_libraries( gbwtgraph PUBLIC gbwt "${HANDLEGRAPH_LIB}/libhandlegraph.a")
    add_dependencies( gbwtgraph gbwt libhandlegraph_bdsg )
else()
    add_library( gbwtgraph ${gbwtgraph_SRCS} )
    target_include_directories( gbwtgraph PUBLIC "${deps_DIR}/gbwtgraph/include")
    target_link_libraries( gbwtgraph PUBLIC gbwt ${STOAT_HANDLEGRAPH_TARGET} )
endif()

install(TARGETS gbwtgraph
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    PUBLIC_HEADER DESTINATION include)

add_library(${PROJECT_NAME}_objs OBJECT
    src/utils.cpp
    src/arg_parser.cpp
    src/post_processing.cpp
    src/gaf_creator.cpp
    src/types_and_structs.cpp
    src/matrix.cpp
    src/log.cpp
    src/stats_test.cpp
    src/fast_fishers_exact_test.cpp
    src/writer.cpp
    src/snarl_analyzer.cpp
    src/path_partitioner.cpp
    src/snarl_data_collection.cpp
    src/feature_tables.cpp
    src/vcf_parser.cpp
    src/banner.cpp

    src/io/register_loader_saver_gbz.cpp
    src/io/register_loader_saver_gbzgraph.cpp
    src/io/register_loader_saver_hash_graph.cpp
    src/io/register_loader_saver_packed_graph.cpp
    src/io/register_io.cpp

    src/subcommand/vcf.cpp
    src/subcommand/graph.cpp
    src/subcommand/test.cpp
    src/subcommand/bh_correct.cpp
    src/subcommand/change_reference.cpp
)

# Include Directories
target_include_directories(${PROJECT_NAME}_objs PRIVATE
${Boost_INCLUDE_DIRS}
${HTSLIB_INCLUDE_DIRS}
)
if (STATIC_BUILD)
    target_include_directories(${PROJECT_NAME}_objs PRIVATE ${STATIC_DEPS_INCLUDE_DIRS})
endif()

# Links a target against the Stoat object library plus the right set of
# bdsg/gbwt/gbwtgraph/vgio/htslib/OpenMP libraries for the current build mode.
# Centralizing this avoids repeating the (long) static link-library list on
# every one of the six executables below.
function(stoat_link_target target)
    if (STATIC_BUILD)
        target_link_options(${target} PRIVATE "-Wl,--start-group")
        target_link_libraries(${target} PRIVATE
            ${PROJECT_NAME}_objs
            ${STOAT_BDSG_TARGET}
            "${SDSL_LIB}/libsdsl.a"
            "${SDSL_LIB}/libdivsufsort.a"
            "${SDSL_LIB}/libdivsufsort64.a"
            "${HANDLEGRAPH_LIB}/libhandlegraph.a"
            gbwtgraph
            gbwt
            ${STOAT_VGIO_TARGET}
            Boost::boost
            ${STOAT_HTSLIB_LIBRARIES}
            OpenMP::OpenMP_CXX
            ${STOAT_STATIC_SYSTEM_LIBS}
            ${ARGN}
        )
        target_link_options(${target} PRIVATE "-Wl,--end-group")
        add_dependencies(${target} sdsl_lite libhandlegraph_bdsg libbdsg)
    else()
        target_link_libraries(${target} PRIVATE
            ${PROJECT_NAME}_objs
            Boost::boost
            ${STOAT_HTSLIB_LIBRARIES}
            ${STOAT_BDSG_TARGET}
            gbwtgraph
            ${STOAT_VGIO_TARGET}
            ${ARGN}
        )
    endif()
endfunction()

# Add Executable
add_executable(${PROJECT_NAME} src/main.cpp)
stoat_link_target(${PROJECT_NAME})
target_compile_options(${PROJECT_NAME} PRIVATE -w) # fix warning

# Enable testing
enable_testing()

# Add test executable
add_executable(unit_tests
    tests/unittest/main_test.cpp
    tests/unittest/binary_test_unit.cpp
    tests/unittest/gaf_creator_unit.cpp
    tests/unittest/types_and_structs_unit.cpp
    tests/unittest/matrix_unit.cpp
    tests/unittest/quantitative_test_unit.cpp
    tests/unittest/path_partitioner_unit.cpp
    tests/unittest/arg_parsing_unit.cpp
    tests/unittest/post_processing_unit.cpp
    tests/unittest/snarl_analyzer_unit.cpp
    tests/unittest/utils_unit.cpp
    tests/unittest/snarl_data_collection_unit.cpp
    tests/unittest/feature_tables_unit.cpp
    tests/unittest/vcf_parser_unit.cpp
)

target_compile_options(unit_tests PRIVATE -w) # fix warning
stoat_link_target(unit_tests)

include_directories(${deps_DIR}/Catch2/single_include/catch2 )

add_test(NAME UnitTests COMMAND unit_tests)

# Add test executable system tests
add_executable(vcf_simu_test tests/system/compare_files_utils.cpp tests/system/vcf_simu_test.cpp)
stoat_link_target(vcf_simu_test)
add_test(NAME VcfSimuTest COMMAND vcf_simu_test)

add_executable(graph_simu_test tests/system/compare_files_utils.cpp tests/system/load_tables.cpp tests/system/graph_simu_test.cpp)
stoat_link_target(graph_simu_test)
add_test(NAME GraphSimuTest COMMAND graph_simu_test)

add_executable(test_simu_test tests/system/compare_files_utils.cpp tests/system/load_tables.cpp tests/system/test_simu_test.cpp)
stoat_link_target(test_simu_test)
add_test(NAME TestSimuTest COMMAND test_simu_test)


##### Add changing reference system test
add_executable(change_ref_test tests/system/compare_files_utils.cpp tests/system/change_reference_test.cpp)
stoat_link_target(change_ref_test)
add_test(NAME ChangeRefTest COMMAND change_ref_test)