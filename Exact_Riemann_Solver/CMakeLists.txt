cmake_minimum_required(VERSION 3.29)
project(WorkSpace)

set(CMAKE_CXX_STANDARD 20)

add_executable(WorkSpace main.cpp
        RiemannSolver.cpp
        RiemannSolver.h
)
