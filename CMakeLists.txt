cmake_minimum_required(VERSION 3.10)
project(wave_equation)

find_package(deal.II 9.0 REQUIRED)

add_executable(wave_equation main.cpp WaveEquation.cpp)

deal_ii_setup_target(wave_equation)

target_link_libraries(wave_equation ${DEAL_II_LIBRARIES})