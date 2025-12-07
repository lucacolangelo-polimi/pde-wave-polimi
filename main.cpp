#include "include/WaveEquation.h"
#include <iostream>

// ----------------------------------------------------------------------------
// BLOCK: Main
// ----------------------------------------------------------------------------
int main()
{
    try
    {
        std::cout << "Starting Wave Equation Solver (deal.II)..." << std::endl;
        
        // Instantiate the problem in 2 Dimensions
        WaveEquation<2> wave_problem;
        
        // Lancia la simulazione
        wave_problem.run();
    }
    catch (std::exception &exc)
    {
        std::cerr << std::endl
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
                  << "----------------------------------------------------"
                  << std::endl
                  << "Unknown exception!" << std::endl
                  << "Aborting!" << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        return 1;
    }
    return 0;
}