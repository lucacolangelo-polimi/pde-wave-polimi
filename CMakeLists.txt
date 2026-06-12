cmake_minimum_required(VERSION 3.13)
project(WaveEquation)
 
find_package(deal.II 9.3 REQUIRED
  HINTS ${DEAL_II_DIR} $ENV{DEAL_II_DIR}
)
deal_ii_initialize_cached_variables()

add_executable(wave_equation
  mainParallel.cpp
  WaveEquationParallel.cpp
)
deal_ii_setup_target(wave_equation)

add_executable(wave_equation_standard
  main.cpp
  WaveEquation.cpp
)
deal_ii_setup_target(wave_equation_standard)

add_executable(wave_equation_optimized
  mainOptimized.cpp
  WaveEquationOptimized.cpp
)
deal_ii_setup_target(wave_equation_optimized)

message(STATUS "deal.II trovato in: ${DEAL_II_DIR}")
message(STATUS "MPI abilitato: ${DEAL_II_WITH_MPI}")
message(STATUS "Trilinos abilitato: ${DEAL_II_WITH_TRILINOS}")
message(STATUS "p4est abilitato: ${DEAL_II_WITH_P4EST}")

if(NOT DEAL_II_WITH_MPI OR NOT DEAL_II_WITH_TRILINOS OR NOT DEAL_II_WITH_P4EST)
  message(FATAL_ERROR
    "deal.II deve essere compilato con MPI + Trilinos + p4est per questo solver.\n"
    "Controlla la tua installazione della libreria."
  )
endif()
add_executable(test_st test_st2.cpp)
deal_ii_setup_target(test_st)
