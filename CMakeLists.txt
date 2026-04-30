cmake_minimum_required(VERSION 3.13.4)
project(WaveEquation)

# Define the build type (Release for performance, Debug for troubleshooting)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

# Find the deal.II library
find_package(deal.II 9.0.0 REQUIRED HINTS ${DEAL_II_DIR})

# Define the executable and its source files
add_executable(wave_equation main.cpp WaveEquation.cpp)
deal_ii_setup_target(wave_equation)