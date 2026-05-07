cmake_minimum_required(VERSION 3.10)
project(wave_equation3D)

find_package(deal.II 9.0 REQUIRED)

add_executable(wave_equation3D main.cpp WaveEquation.cpp)

deal_ii_setup_target(wave_equation3D)

target_link_libraries(wave_equation3D ${DEAL_II_LIBRARIES})