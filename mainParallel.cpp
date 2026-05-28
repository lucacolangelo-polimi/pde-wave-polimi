#include "WaveEquationParallel.hpp"

int main(int argc, char *argv[])
{
    Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);
    MPI_Comm mpi_comm = MPI_COMM_WORLD;

    try
    {
        // Parsing args (simple manual parsing, can be improved with a library)
        std::string mode_str = "run";
        unsigned int ref = 6;
        for (int i = 1; i < argc; ++i)
        {
            if (std::string(argv[i]) == "--mode" && i + 1 < argc) mode_str = argv[++i];
            if (std::string(argv[i]) == "--ref" && i + 1 < argc)  ref = std::stoi(argv[++i]);
        }

        // Setup 
        WaveEquation<2> wave(mpi_comm);
        wave.initial_refinement = ref;
        wave.time_step = 5.0e-4;
        wave.mode = SimulationMode::PEBBLE_IN_POND;
        wave.use_amr = false; // AMR off for a clean scaling test

        // Dispatcher based on mode
        if (mode_str == "strong" || mode_str == "weak")
        {
            wave.prepare_for_analysis(); 
            auto result = wave.run_scaling_benchmark(100);
            wave.print_scaling_table({result}, "benchmark_results.txt");
        }
        else if (mode_str == "dispersion")
        {
            wave.prepare_for_analysis(); 
            const double h = 1.0 / std::pow(2.0, wave.initial_refinement);
            std::vector<double> kvals;
            for (int i = 1; i <= 8; ++i) kvals.push_back(M_PI * i / (8.0 * h));
            wave.run_dispersion_analysis(kvals);
        }
        else
        {
            // standard run
            wave.run();
        }
    }
    catch (std::exception &exc)
    {
        std::cerr << "Errore: " << exc.what() << "\n";
        return 1;
    }

    return 0;
}
