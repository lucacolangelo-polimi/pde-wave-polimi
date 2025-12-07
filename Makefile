cmake_minimum_required(VERSION 3.13)

# ----------------------------------------------------------------------------
# TODO: You can change the project name here if desired
# ----------------------------------------------------------------------------
project(wave_solver_dealii)

# ----------------------------------------------------------------------------
# BLOCK: deal.II Configuration
# Finds the library installed on the system. Do not edit.
# ----------------------------------------------------------------------------
find_package(deal.II 9.0.0 REQUIRED)
deal_ii_initialize_cached_variables()

# ----------------------------------------------------------------------------
# BLOCK: Sources and Executable Definition
# Add any other .cc files here if the project expands
# ----------------------------------------------------------------------------
add_executable(wave_solver
    main.cc
    source/WaveEquation.cc
)

# Link deal.II libraries to the executable
target_link_libraries(wave_solver ${DEAL_II_LIBRARIES})

# Copy .msh files to the build directory (optional, useful if using GMSH)
# file(COPY ${CMAKE_SOURCE_DIR}/mesh.msh DESTINATION ${CMAKE_BINARY_DIR})