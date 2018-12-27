# Created by Christian Bach, christiangoebelbach@gmail.com

cmake_minimum_required (VERSION 3.10)
project (Solver)

# Find Armadillo on system
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

# Compiler flags to optimise performance
set(CMAKE_CXX_FLAGS "-Ofast -march=native")

# Recursively builds all the needed components in subdirectories
add_subdirectory(Interpolation)
add_subdirectory(Domain)

# Create Solver target
set(Solver_Source
    main.cpp
    Domain/Quadratures.h)

add_executable(Solver ${Solver_Source})

# Link all libs and target
target_link_libraries(Solver ${ARMADILLO_LIBRARIES} Interpolation Element)
