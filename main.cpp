#include "WaveEquation.hpp"

int main(int argc, char *argv[])
{
  using namespace dealii;

  try
  {

    Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
    WaveEquation<2> wave_problem;
    wave_problem.run();
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