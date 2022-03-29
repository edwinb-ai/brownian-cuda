cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
message(STATUS "CMake version: ${CMAKE_VERSION}")

project(cufort Fortran)

# Assume the compiler is NVHPC
set(dialect "-cuda -fast -cudalib=curand")

# Compilation steps
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")
message(STATUS "Fortran compiler flags=${CMAKE_Fortran_FLAGS}")

file(GLOB_RECURSE sources include/*.f90 src/*.f90)
add_executable(cufort ${sources})