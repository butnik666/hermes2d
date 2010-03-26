if(NOT REAL)
    return()
endif(NOT REAL)
project(singular-perturbation)

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)
