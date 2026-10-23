# calx: numerical methods and the computation stack.
#
# Every module is a target calx::<module>. Built on its own, calx also builds its checks
# (CALX_CHECKS); added to a larger project (the infrastructure's root, physis, a domain),
# it builds them only when asked.

cmake_minimum_required(VERSION 3.20)
project(calx VERSION 0.1.0 LANGUAGES CXX)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(CALX_TOP_LEVEL ON)
else()
  set(CALX_TOP_LEVEL OFF)
endif()
option(CALX_CHECKS "Build calx's checks" ${CALX_TOP_LEVEL})

if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 20)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
endif()
if(CALX_CHECKS)
  enable_testing()
endif()

# numerics/verification first: it defines calx_add_check() for the others.
add_subdirectory(numerics/verification)
add_subdirectory(numerics/random)
add_subdirectory(numerics/linalg)
add_subdirectory(numerics/integrators)
add_subdirectory(numerics/special_functions)
add_subdirectory(numerics/quadrature)
add_subdirectory(numerics/solvers)
add_subdirectory(compute/scheduler)
