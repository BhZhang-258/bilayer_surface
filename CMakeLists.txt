cmake_minimum_required(VERSION 3.10)
project(simDER)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -fopenmp")

find_package(Eigen3 REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

set(MKLROOT /opt/intel/oneapi/mkl/latest)

include_directories(
    ${EIGEN3_INCLUDE_DIR}
    ${MKLROOT}/include
)

link_directories(
    ${MKLROOT}/lib/intel64
)

file(GLOB SOURCES "*.cpp")

add_executable(simDER ${SOURCES})

target_link_libraries(simDER
    ${OPENGL_LIBRARIES}
    ${GLUT_LIBRARIES}
    mkl_intel_lp64
    mkl_intel_thread
    mkl_core
    iomp5
    pthread
    m
)