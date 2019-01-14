##### WARNING!! ##################################
# This is not yet functional
# it is designed for me (MES) to be able to use the CLion IDE
# and its refactorization tools.  If anyone would like to
# help get this functional, please contact us!
#
# In particular: one currently needs the built include files
# from our build tree (see CMakeLists.txt in Macaulay2/e)
##################################################

cmake_minimum_required(VERSION 3.10)
project(Macaulay2 VERSION 1.12.0.1)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

####################################################
# header directory includes (CLion takes a while to index)
####################################################
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/include)
#
include_directories(${CMAKE_SOURCE_DIR}/usr-host/include)
include_directories(${CMAKE_SOURCE_DIR}/usr-host/include/gtest/include)
#include_directories(/usr/local/include)  # not sure why not included by default
#
include_directories(${CMAKE_SOURCE_DIR}/Macaulay2/c)
include_directories(${CMAKE_SOURCE_DIR}/Macaulay2/d)
include_directories(${CMAKE_SOURCE_DIR}/Macaulay2/e)
include_directories(${CMAKE_SOURCE_DIR}/Macaulay2/e/f4)
include_directories(${CMAKE_SOURCE_DIR}/Macaulay2/kernel/bibasis)
include_directories(${CMAKE_SOURCE_DIR}/Macaulay2/system)

####################################################
# locate and include external libraries
####################################################
#find_package(PkgConfig REQUIRED)
#pkg_search_module(GTest REQUIRED gtest_main)
##
#if(GTest_FOUND)
#    message("GTEST IS FOUND!")
#else()
#    message("GTEST IS NOT FOUND")
#endif()

# target_include_directories(m2 PUBLIC
#     "BUILD/mike/builds.tmp/darwin64-gcc4.9/include"
#     "BUILD/mike/builds.tmp/darwin64-gcc4.9/Macaulay2/d"
#     "BUILD/mike/builds.tmp/darwin64-gcc4.9/Macaulay2/e-includes"
#     include
#     Macaulay2/e
#     Macaulay2/c
#     "${PROJECT_BINARY_DIR}"
#     )                            

# set(E-FILES
#     Macaulay2/e/comp.cpp
# )
# set(SOURCE_FILES
#     Macaulay2/e/comp.cpp
#     Macaulay2/e/unit-tests/testMain.cpp
#     )

# add_library(m2 STATIC ${E-FILES})
# add_library(m2 STATIC class.cpp)
# target_link_libraries(m2 STATIC )



# enable_testing()
add_executable(macaulay2-test Macaulay2/e/unit-tests/testMain.cpp)
target_link_libraries(macaulay2-test ${GTest_BOTH_LIBRARIES})
add_test(AllTestsInFoo macaulay2-test)

# find_package(GTest REQUIRED)
# include_directories(${GTEST_INCLUDE_DIRS})
# target_link_libraries(macaulay2-test ${GTEST_BOTH_LIBRARIES})
# add_test(AllTestsInFoo macaulay2-test)


# add_executable(proj-test
#   ${source_files}
#   ${test_files}
#   )

# set(GTEST_INCLUDE_DIRS
#   gtest-1.7.0/include
#   gtest-1.7.0
#   )

# target_include_directories(proj-test PUBLIC
#   ${GTEST_INCLUDE_DIRS}
#   include
#   )


# include_directories("${PROJECT_BINARY_DIR}")
# include_directories(
#     "BUILD/mike/builds.tmp/darwin64-gcc4.9/include"
#     "BUILD/mike/builds.tmp/darwin64-gcc4.9/Macaulay2/d"
#     "BUILD/mike/builds.tmp/darwin64-gcc4.9/Macaulay2/e-includes"
#     include
#     Macaulay2/e
#     Macaulay2/c
# )

