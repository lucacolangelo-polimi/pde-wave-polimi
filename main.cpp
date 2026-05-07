#include "WaveEquation.hpp"

int main(int argc, char *argv[])
{
  using namespace dealii;

  try
  {

    Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
    WaveEquation<3> wave_problem3D; // Create an instance of the WaveEquation class for 3D
    wave_problem3D.run();
  }
  
  catch (std::exception &exc)
  {
    std::cerr << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  }
  catch (...)
  {
    std::cerr << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  }

  return 0;
}