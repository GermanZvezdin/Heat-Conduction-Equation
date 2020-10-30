cmake_minimum_required(VERSION 3.15)
project(Heat-Conduction-Equation)

set(CMAKE_CXX_STANDARD 20)

add_executable(Heat-Conduction-Equation z-arrays.hpp heat.cpp)