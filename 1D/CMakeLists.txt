cmake_minimum_required(VERSION 2.8)
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
add_executable(fdtd_1D FDHeatEquation1D.cc)
target_link_libraries(fdtd_1D ${ARMADILLO_LIBRARIES})
