cmake_minimum_required(VERSION 3.13)
project(stoat VERSION 0.0.4 LANGUAGES CXX)

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

# To change to OPTIMIZE=OFF use: cmake -DOPTIMIZE=OFF ..
option(OPTIMIZE "Build with optimization [default: on]" ON)
message(STATUS "OPTIMIZE is set to: ${OPTIMIZE}")

# Option to produce a fully statically linked executable (Linux only)
option(STATIC_BUILD "Build a statically linked executable and force static libs for subprojects (Linux only)" OFF)
message(STATUS "STATIC_BUILD is set to: ${STATIC_BUILD}")

# Static build configuration
if(STATIC_BUILD)

    # Force CMake subprojects to build static libraries by default
    set(BUILD_SHARED_LIBS OFF CACHE BOOL
        "Build static libraries by default" FORCE)

    # Prefer static archives when resolving dependencies
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")

    # Static libraries should be built with PIC where possible
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)

    if(UNIX AND NOT APPLE)

        # Fully static Linux executable
        set(CMAKE_EXE_LINKER_FLAGS
            "${CMAKE_EXE_LINKER_FLAGS} -static -s -pthread -static-libgcc -static-libstdc++"
        )

        message(STATUS "Configured for FULL STATIC Linux build")

    elseif(APPLE)

        # Build third-party dependencies as .a where possible.
        message(STATUS "Configured for STATIC third-party libraries on macOS")

    endif()

endif()

find_package(PkgConfig REQUIRED)
pkg_check_modules(HTSLIB REQUIRED htslib)
pkg_check_modules(Jansson REQUIRED IMPORTED_TARGET jansson)

include_directories(${HTSLIB_INCLUDE_DIRS})
link_directories(${HTSLIB_LIBRARY_DIRS})

# Find Packages
find_package(Boost REQUIRED)
find_package(OpenMP REQUIRED)

if(STATIC_BUILD)
    find_library(BZ2_LIB NAMES bz2 REQUIRED)

    set(STATIC_SYSTEM_LIBS
        ${BZ2_LIB}
        z
        pthread
        dl
        m
    )
endif()

# Get eigen3
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
    # Try multiple possible locations for libomp
    set(LIBOMP_PATHS "${HOMEBREW_PREFIX}/opt/libomp" "/usr/local/opt/libomp")

    set(LIBOMP_FOUND FALSE)
    foreach(LIBOMP_PATH ${LIBOMP_PATHS})
        if(EXISTS "${LIBOMP_PATH}/include/omp.h")
            set(OpenMP_INCLUDES "${LIBOMP_PATH}/include")
            set(OpenMP_LIBRARIES "${LIBOMP_PATH}/lib/libomp.dylib")
            set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${OpenMP_INCLUDES}")
            set(OpenMP_LINK_FLAGS "-L${LIBOMP_PATH}/lib -lomp")
            set(LIBOMP_FOUND TRUE)
            message(STATUS "Found libomp at: ${LIBOMP_PATH}")
            break()
        endif()
    endforeach()

    if(LIBOMP_FOUND)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_LINK_FLAGS}")
    else()
        message(STATUS "Warning: libomp not found in standard Homebrew locations. CMake will try to find OpenMP using standard search paths.")
    endif()
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

add_subdirectory(${deps_DIR}/libvgio EXCLUDE_FROM_ALL)

# Comes along with its own libhandlegraph, which it depends on.
# Make sure to never build its Python bindings.
set(BUILD_PYTHON_BINDINGS OFF CACHE BOOL "Build Python bindings for libbdsg" FORCE)

set(RUN_DOXYGEN OFF CACHE BOOL "Run libbdsg Doxygen" FORCE)

if (STATIC_BUILD)

    set(DEPS           "${CMAKE_SOURCE_DIR}/deps")
    set(LIBBDSG_SRC    "${DEPS}/libbdsg")
    set(BDSG           "${LIBBDSG_SRC}/bdsg")

    set(DEPS_INCLUDE    "${BDSG}/include")
    set(DEPS_LIBBDSG    "${BDSG}/deps")
    set(DEPS_LIB        "${DEPS}/lib")
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

    # Include directories for static dependencies
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

    set(COMMON_CXX_FLAGS "-fPIC -O3 -std=c++17 ${CMAKE_CXX_FLAGS}")

    # 1. sdsl-lite
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

    # 2. libhandlegraph
    ExternalProject_Add(libhandlegraph
        SOURCE_DIR ${LIBHANDLEGRAPH_SRC}
        CMAKE_ARGS
            -DCMAKE_INSTALL_PREFIX=${DEPS_LIBBDSG}
            -DCMAKE_INSTALL_LIBDIR=lib
            -DCMAKE_CXX_FLAGS=${COMMON_CXX_FLAGS}
            -DCMAKE_BUILD_TYPE=Release
        BUILD_BYPRODUCTS ${HANDLEGRAPH_LIB}/libhandlegraph.a
    )

    # 3. sparsepp
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

    # 5. mio
    ExternalProject_Add(mio
        SOURCE_DIR ${MIO_SRC}
        DOWNLOAD_COMMAND ""
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${MIO_SRC}/include/mio ${DEPS_INCLUDE}/mio
    )

    # 6. BBHash
    ExternalProject_Add(bbhash
        SOURCE_DIR ${BBHASH_SRC}
        DOWNLOAD_COMMAND ""
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ${CMAKE_COMMAND} -E copy
            ${BBHASH_SRC}/BooPHF.h
            ${DEPS_INCLUDE}/BooPHF.h
    )

    # 7. libbdsg
    ExternalProject_Add(libbdsg
        SOURCE_DIR ${LIBBDSG_SRC}
        DEPENDS
            sdsl_lite
            libhandlegraph
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

else()

    # Dynamic libbdsg and its dependencies
    add_subdirectory(
        ${deps_DIR}/libbdsg
        EXCLUDE_FROM_ALL
    )

    set(DYNAMIC_DEPS_INCLUDE_DIRS
        "${deps_DIR}/libbdsg/bdsg/include"
        "${deps_DIR}/libbdsg/bdsg/deps/libhandlegraph/src/include"
        "${deps_DIR}/libbdsg/bdsg/deps/sdsl-lite/include"
        "${deps_DIR}/libbdsg/bdsg/deps/sparsepp"
        "${deps_DIR}/libbdsg/bdsg/deps/mio/include"
        "${deps_DIR}/libbdsg/bdsg/deps/BBHash"
        "${deps_DIR}/libbdsg/bdsg/deps/DYNAMIC/include"
        "${deps_DIR}/libbdsg/bdsg/deps/DYNAMIC/deps/hopscotch_map/include"
        "${deps_DIR}/libhandlegraph/src/include"
        "${deps_DIR}/libvgio/include"
        "${deps_DIR}/gbwt/include"
        "${deps_DIR}/gbwtgraph/include"
    )

    set(DYNAMIC_DEPS_LIBRARIES
        Boost::boost
        libbdsg
        gbwtgraph
        gbwt
        handlegraph_shared
        vgio
        sdsl
        OpenMP::OpenMP_CXX
        ${HTSLIB_LIBRARIES}
    )

endif()

if(STATIC_BUILD)
    set(STOAT_STATIC_LINK_LIBS
        bdsg::static
        ${SDSL_LIB}/libsdsl.a
        ${SDSL_LIB}/libdivsufsort.a
        ${SDSL_LIB}/libdivsufsort64.a
        ${HANDLEGRAPH_LIB}/libhandlegraph.a
        gbwtgraph
        gbwt
        vgio_static
        ${HTSLIB_LIBRARIES}
        OpenMP::OpenMP_CXX
        ${STATIC_SYSTEM_LIBS}
    )
endif()

# Links target against the Stoat object library
function(stoat_link_target target)
    if(STATIC_BUILD)
        target_link_options(${target} PRIVATE "-Wl,--start-group")

        target_link_libraries(${target} PRIVATE
            ${PROJECT_NAME}_objs
            ${STOAT_STATIC_LINK_LIBS}
            ${ARGN}
        )

        target_link_options(${target} PRIVATE "-Wl,--end-group")
        add_dependencies(${target} sdsl_lite libhandlegraph libbdsg)
        target_include_directories(${target} PRIVATE ${STATIC_DEPS_INCLUDE_DIRS})
    
    else()
        target_link_libraries(${target} PRIVATE
            ${PROJECT_NAME}_objs
            ${DYNAMIC_DEPS_LIBRARIES}
        )

        target_include_directories(${target} PRIVATE ${DYNAMIC_DEPS_INCLUDE_DIRS})

    endif()
endfunction()

# Declares a test executable from the given sources, links it via
# stoat_link_target(), and registers it with CTest as `test_name`.
function(stoat_add_test_executable target test_name)
    add_executable(${target} ${ARGN})
    stoat_link_target(${target})
    add_test(NAME ${test_name} COMMAND ${target})
endfunction()

# Same as stoat_add_test_executable(), but also links Boost (needed by the
# shared tests/system/*.cpp helpers) when STATIC_BUILD is on.
function(stoat_add_system_test target test_name)
    add_executable(${target} ${ARGN})
    stoat_link_target(${target} Boost::boost)
    add_test(NAME ${test_name} COMMAND ${target})
endfunction()

# GBWT
foreach (gbwt_source IN ITEMS algorithms.cpp bwtmerge.cpp cached_gbwt.cpp dynamic_gbwt.cpp fast_locate.cpp files.cpp gbwt.cpp internal.cpp metadata.cpp support.cpp test.cpp utils.cpp variants.cpp)
    list(APPEND gbwt_libFiles "deps/gbwt/src/${gbwt_source}")
endforeach()

if(MSVC)
    file(GLOB gbwt_headerFiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "${deps_DIR}/gbwt/include/gbwt/*.h")
endif()

set(gbwt_SRCS ${gbwt_libFiles} ${gbwt_headerFiles})
add_library(gbwt ${gbwt_SRCS})

if(STATIC_BUILD)
    target_include_directories(gbwt PUBLIC
        "${deps_DIR}/gbwt/include"
        ${STATIC_DEPS_INCLUDE_DIRS}
    )

    target_link_libraries(gbwt PUBLIC
        "${SDSL_LIB}/libsdsl.a"
        OpenMP::OpenMP_CXX
    )

    add_dependencies(gbwt sdsl_lite)

else()

    target_include_directories(gbwt PUBLIC
        "${deps_DIR}/gbwt/include"
        "${deps_DIR}/libbdsg/bdsg/deps/sdsl-lite/include"
    )

    target_link_libraries(gbwt PUBLIC
        sdsl
        OpenMP::OpenMP_CXX
    )

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
    file(GLOB gbwtgraph_headerFiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "${deps_DIR}/gbwtgraph/include/gbwt/*.h")
endif()

set(gbwtgraph_SRCS ${gbwtgraph_libFiles} ${gbwtgraph_headerFiles})
add_library(gbwtgraph ${gbwtgraph_SRCS})
target_include_directories(gbwtgraph PUBLIC "${deps_DIR}/gbwtgraph/include")

if(STATIC_BUILD)
    target_include_directories(gbwtgraph PUBLIC
        "${deps_DIR}/gbwtgraph/include"
        "${deps_DIR}/gbwt/include"
        ${STATIC_DEPS_INCLUDE_DIRS}
    )

    target_link_libraries(gbwtgraph PUBLIC
        gbwt
        "${HANDLEGRAPH_LIB}/libhandlegraph.a"
    )

    add_dependencies(gbwtgraph gbwt)

else()
    target_link_libraries(gbwtgraph PUBLIC
        gbwt
        handlegraph_shared
    )
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

if(STATIC_BUILD)

    target_include_directories(${PROJECT_NAME}_objs PRIVATE
        ${STATIC_DEPS_INCLUDE_DIRS}
    )

else()

    target_include_directories(${PROJECT_NAME}_objs PRIVATE
        # libbdsg
        "${deps_DIR}/libbdsg/bdsg/include"

        # libhandlegraph used by libbdsg
        "${deps_DIR}/libbdsg/bdsg/deps/libhandlegraph/src/include"

        # SDSL
        "${deps_DIR}/libbdsg/bdsg/deps/sdsl-lite/include"

        # libbdsg header-only dependencies
        "${deps_DIR}/libbdsg/bdsg/deps/sparsepp"
        "${deps_DIR}/libbdsg/bdsg/deps/mio/include"
        "${deps_DIR}/libbdsg/bdsg/deps/BBHash"
        "${deps_DIR}/libbdsg/bdsg/deps/DYNAMIC/include"
        "${deps_DIR}/libbdsg/bdsg/deps/DYNAMIC/deps/hopscotch_map/include"

        # Other Stoat dependencies
        "${deps_DIR}/libhandlegraph/src/include"
        "${deps_DIR}/libvgio/include"
        "${deps_DIR}/libvgio/deps/libhandlegraph/src/include"
        "${deps_DIR}/gbwt/include"
        "${deps_DIR}/gbwtgraph/include"
    )

endif()

target_link_libraries(${PROJECT_NAME}_objs PRIVATE
    OpenMP::OpenMP_CXX
)

# Add Executable
add_executable(${PROJECT_NAME} src/main.cpp)
stoat_link_target(${PROJECT_NAME})
target_compile_options(${PROJECT_NAME} PRIVATE -w) # fix warning

# Enable testing
enable_testing()

# Add test executable
stoat_add_test_executable(unit_tests UnitTests
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

target_compile_options(unit_tests PRIVATE -w)

include_directories(${deps_DIR}/Catch2/single_include/catch2)

# Add test executable system tests
stoat_add_system_test(vcf_simu_test VcfSimuTest
    tests/system/compare_files_utils.cpp
    tests/system/vcf_simu_test.cpp
)

stoat_add_system_test(graph_simu_test GraphSimuTest
    tests/system/compare_files_utils.cpp
    tests/system/load_tables.cpp
    tests/system/graph_simu_test.cpp
)

stoat_add_system_test(test_simu_test TestSimuTest
    tests/system/compare_files_utils.cpp
    tests/system/load_tables.cpp
    tests/system/test_simu_test.cpp
)

stoat_add_system_test(bh_correct_test BHCorrectTest
    tests/system/compare_files_utils.cpp
    tests/system/load_tables.cpp
    tests/system/bh_correct_test.cpp
)

stoat_add_system_test(change_ref_test ChangeRefTest
    tests/system/compare_files_utils.cpp
    tests/system/change_reference_test.cpp
)