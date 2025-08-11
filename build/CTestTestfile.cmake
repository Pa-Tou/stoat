# CMake generated Testfile for 
# Source directory: /home/mbagarre/Bureau/STOAT
# Build directory: /home/mbagarre/Bureau/STOAT/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(UnitTests "/home/mbagarre/Bureau/STOAT/bin/unit_tests")
set_tests_properties(UnitTests PROPERTIES  _BACKTRACE_TRIPLES "/home/mbagarre/Bureau/STOAT/CMakeLists.txt;203;add_test;/home/mbagarre/Bureau/STOAT/CMakeLists.txt;0;")
add_test(VcfSimuTest "/home/mbagarre/Bureau/STOAT/bin/vcf_simu_test")
set_tests_properties(VcfSimuTest PROPERTIES  _BACKTRACE_TRIPLES "/home/mbagarre/Bureau/STOAT/CMakeLists.txt;213;add_test;/home/mbagarre/Bureau/STOAT/CMakeLists.txt;0;")
add_test(GraphSimuTest "/home/mbagarre/Bureau/STOAT/bin/graph_simu_test")
set_tests_properties(GraphSimuTest PROPERTIES  _BACKTRACE_TRIPLES "/home/mbagarre/Bureau/STOAT/CMakeLists.txt;223;add_test;/home/mbagarre/Bureau/STOAT/CMakeLists.txt;0;")
subdirs("deps/libhandlegraph")
subdirs("deps/libvgio")
subdirs("deps/libbdsg")
