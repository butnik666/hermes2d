if(NOT REAL)
    return()
endif(NOT REAL)
project(benchmark-lshape)

add_executable(${PROJECT_NAME} main.cpp)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(benchmark-lshape ${BIN})

