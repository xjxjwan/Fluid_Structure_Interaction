cmake_minimum_required(VERSION 3.29)
project(WorkSpace)

set(CMAKE_CXX_STANDARD 20)

add_executable(WorkSpace main.cpp
        RiemannSolver.cpp
        RiemannSolver.h
        DataReconstruct.cpp
        DataReconstruct.h
        CalTimeStep.cpp
        CalTimeStep.h
        HalfTimeStepUpdate.cpp
        HalfTimeStepUpdate.h
        CalFlux.cpp
        CalFlux.h
        SetDomainBoundary.cpp
        SetDomainBoundary.h
        LevelSetFunctions.cpp
        LevelSetFunctions.h
        FastSweeping.cpp
        FastSweeping.h
        ConstantExtrapolation.cpp
        ConstantExtrapolation.h
        RiemannGFM.cpp
        RiemannGFM.h
        AuxiliaryFunctions.cpp
        AuxiliaryFunctions.h
)
