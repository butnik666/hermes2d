if(NOT REAL)
    return()
endif(NOT REAL)

project(sea-breeze)

find_package(PythonLibs REQUIRED)
find_package(NumPy REQUIRED)
include(UseCython)

include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_PATH})
include_directories(${hermes2d_SOURCE_DIR}/python/hermes2d)

add_executable(${PROJECT_NAME} main.cpp forms.cpp numerical_flux.cpp)
set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS_DEBUG})
include (../CMake.common)

target_link_libraries(${PROJECT_NAME} ${PYTHON_LIBRARIES})
add_dependencies(${PROJECT_NAME} _hermes2d)

cython_add_module(_forms forms.cpp)
target_link_libraries(_forms ${HERMES_REAL_BIN})
add_dependencies(_forms _hermes2d)

cython_add_module(_numerical_flux numerical_flux.cpp forms.cpp)
target_link_libraries(_numerical_flux ${HERMES_REAL_BIN})
add_dependencies(_numerical_flux _hermes2d)

