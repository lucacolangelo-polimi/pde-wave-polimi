#include "WaveEquationParallel.hpp"             //last version

int main(int argc, char *argv[])
{
    Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);
    MPI_Comm mpi_comm = MPI_COMM_WORLD;

    try
    {
        // Parsing command-line parameters
        std::string mode_str = "mms";       // Default: MMS convergence study
        unsigned int ref = 4;               // Initial baseline refinement level
        
        for (int i = 1; i < argc; ++i)
        {
            if (std::string(argv[i]) == "--mode" && i + 1 < argc) mode_str = argv[++i];
            if (std::string(argv[i]) == "--ref" && i + 1 < argc)  ref = std::stoi(argv[++i]);
        }

        WaveEquation<2> wave(mpi_comm);
        wave.initial_refinement = ref;

        // =========================================================================
        // 1. STRONG / WEAK SCALING
        // =========================================================================
        if (mode_str == "strong" || mode_str == "weak")
        {
            wave.time_scheme = TimeScheme::NEWMARK; // Newmark used for scaling benchmarks;        
            wave.time_step = 5.0e-4;
            wave.mode = SimulationMode::PEBBLE_IN_POND;
            wave.use_amr = false; // AMR turned off for clean scaling tests

            wave.prepare_for_analysis(); 
            auto result = wave.run_scaling_benchmark(100);
            
            // Save a dedicated CSV file based on the chosen mode
            wave.print_scaling_table({result}, mode_str + "_scaling.csv");
        }
        // =========================================================================
        // 2. NUMERICAL DISPERSION ANALYSIS MODE
        // =========================================================================
        else if (mode_str == "dispersion")
        {
            wave.time_scheme = TimeScheme::LEAPFROG;
            wave.prepare_for_analysis(); 
            
            const double h = 1.0 / std::pow(2.0, wave.initial_refinement);
            std::vector<double> kvals;
            for (int i = 1; i <= 8; ++i) 
                kvals.push_back(M_PI * i / (8.0 * h));
                
            wave.run_dispersion_analysis(kvals);
        }
        // =========================================================================
        // 3. MATH VERIFICATION MODE (MMS - To be launched with only 1 process)
        // =========================================================================
        else if (mode_str == "mms")
        {
            wave.time_scheme = TimeScheme::NEWMARK;
            wave.run_convergence_study();
        }
        // =========================================================================
        // 4. PEBBLE IN POND (With AMR support)
        // =========================================================================
        else if (mode_str == "pebble")
        {
            wave.mode                = SimulationMode::PEBBLE_IN_POND;
            wave.time_scheme         = TimeScheme::LEAPFROG; // Very fast for physics scenarios
            wave.time_step           = 1.0e-3;
            wave.end_time            = 1.0;
            
            // AMR module configuration (p4est)
            wave.use_amr             = true;
            wave.amr_every_n_steps   = 5;
            wave.max_refinement_level = 7; 
            
            wave.output_every_n_steps = 10; // Generate a .pvtu file every 10 steps
            wave.track_energy         = true;

            wave.run();
        }
        // =========================================================================
        // 5. INTERFERENCE (With AMR support)
        // =========================================================================

        else if (mode_str == "interference")
        {
            wave.mode                = SimulationMode::INTERFERENCE;
            wave.time_scheme         = TimeScheme::LEAPFROG;
            wave.time_step           = 1.0e-3;
            wave.end_time            = 1.0;
            
            // AMR module configuration (p4est)
            wave.use_amr             = true;
            wave.amr_every_n_steps   = 5;
            wave.max_refinement_level = 7;
            
            wave.output_every_n_steps = 10;
            wave.track_energy         = true;

            wave.run();
        }
        // =========================================================================
        // 6. REFRACTION (Eterogeneous medium with AMR support)
        // =========================================================================
        else if (mode_str == "refraction")
        {
            wave.mode                = SimulationMode::REFRACTION;
            wave.time_scheme         = TimeScheme::LEAPFROG;
            // For refraction, the CFL is dictated by c_fast (which is larger than c).
            // We use a smaller dt to ensure we don't blow up the Leapfrog.
            wave.time_step           = 5.0e-4; 
            wave.end_time            = 1.0;
            
            wave.use_amr             = true;
            wave.amr_every_n_steps   = 5;
            wave.max_refinement_level = 7;
            
            wave.output_every_n_steps = 20; // Let's save a little less often so as not to fill up the disk
            wave.track_energy         = true;

            wave.run();
        }
        // =========================================================================
        // 7. DIFFRACTION  (Internal obstacle + slit, with AMR support)
        // =========================================================================
        else if (mode_str == "diffraction")
        {
            wave.mode                = SimulationMode::DIFFRACTION;
            wave.time_scheme         = TimeScheme::LEAPFROG;
            wave.time_step           = 5.0e-4; // dt caution for rigid crack edges
            wave.end_time            = 1.0;
            
            wave.use_amr             = true;
            wave.amr_every_n_steps   = 5;
            wave.max_refinement_level = 7;
            
            wave.output_every_n_steps = 20; 
            wave.track_energy         = true;

            wave.run();
        }
        // =========================================================================
        // 8. DAMPING (Sponge Layers / Absorbing BC)
        // =========================================================================
        else if (mode_str == "damping")
        {

            wave.mode                = SimulationMode::DAMPED_WAVE;
            wave.time_scheme         = TimeScheme::LEAPFROG;
            wave.time_step           = 1.0e-3;
            wave.end_time            = 1.5; // Longer to give the wave time to exit the screen
            
            wave.use_amr             = true;
            wave.amr_every_n_steps   = 5;
            wave.max_refinement_level = 7;
            
            wave.output_every_n_steps = 10; 
            wave.track_energy         = true; // we need it to demonstrate absorption

            wave.run();
        }

        /// =========================================================================
        // 9.  ABSORBING BC (Sommerfeld Boundary Conditions)
        // =========================================================================
        else if (mode_str == "absorbing")
        {
            wave.mode                = SimulationMode::ABSORBING_BC;
            wave.use_absorbing_bc    = true; 
            wave.time_scheme         = TimeScheme::LEAPFROG;
            wave.time_step           = 1.0e-3;
            wave.end_time            = 1.5; 
            
            wave.use_amr             = true;
            wave.amr_every_n_steps   = 5;
            wave.max_refinement_level = 7;
            
            wave.output_every_n_steps = 10; 
            wave.track_energy         = true;

            wave.run();
        }

        // =========================================================================
        // UNRECOGNIZED MODE
        // =========================================================================
        else
        {
            if (Utilities::MPI::this_mpi_process(mpi_comm) == 0)
            {
                std::cout << "Unrecognized mode!\n"
                          << "Choose from: strong, weak, dispersion, mms, pebble, interference, refraction, diffraction, damping, absorbing\n";
            }
        }
    }
    
    catch (std::exception &exc)
    {
        std::cerr << "Error caught in main: " << exc.what() << "\n";
        return 1;
    }

    return 0;
}