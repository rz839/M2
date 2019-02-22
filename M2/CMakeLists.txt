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
project(Macaulay2 VERSION 1.14.0)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

####################################################
# INCLUDE DIRECTORIES
####################################################

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/include)

include_directories(${CMAKE_SOURCE_DIR}/usr-host/include)
#include_directories(/usr/local/include)  # not sure why not included by default

include_directories(${CMAKE_SOURCE_DIR}/Macaulay2/c)
include_directories(${CMAKE_SOURCE_DIR}/Macaulay2/d)
include_directories(${CMAKE_SOURCE_DIR}/Macaulay2/e)
include_directories(${CMAKE_SOURCE_DIR}/Macaulay2/e/f4)
include_directories(${CMAKE_SOURCE_DIR}/Macaulay2/kernel/bibasis)
include_directories(${CMAKE_SOURCE_DIR}/Macaulay2/system)

####################################################
# SET UP GTEST
####################################################

enable_testing()

set(GTEST_ROOT "${CMAKE_SOURCE_DIR}/libraries/gtest/build/googletest-release-1.8.0/googletest" CACHE PATH "")
add_subdirectory(${GTEST_ROOT})
find_package(GTest REQUIRED gtest gtest_main)
include_directories(${GTEST_INCLUDE_DIRS})

####################################################
# UNIT TESTS - put in target lists instead
####################################################

#add_subdirectory(Macaulay2/e/scratch)
add_subdirectory(Macaulay2/e/unit-tests)
add_subdirectory(Macaulay2/e/more-tests)

####################################################
# PHONY TARGET FOR CLION
####################################################

file(GLOB SOURCES
        Macaulay2/e/*.cpp
        Macaulay2/e/*.hpp
        Macaulay2/e/unit-tests/*.cpp
        Macaulay2/e/unit-tests/*.hpp)
add_executable(CLionPhony EXCLUDE_FROM_ALL ${SOURCES})

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


