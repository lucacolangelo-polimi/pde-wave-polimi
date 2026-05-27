
#include "WaveEquationOptimized.hpp"

int main(int argc, char *argv[])
{
    using namespace dealii;

    try
    {
        Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

        // ============================================================
        // Select simulation mode by commenting/uncommenting:
        //
        //  PEBBLE_IN_POND  - Gaussian pulse, Dirichlet BC, AMR
        //  DAMPED_WAVE     - Gaussian pulse with viscous damping
        //  ABSORBING_BC    - Gaussian pulse with Sommerfeld ABC
        //  MMS_CONVERGENCE - Full convergence study (use run_convergence_study())
        // ============================================================

        WaveEquation<2> wave;

        // ---- Physical parameters ----
        wave.c       = 1.0;     // Wave speed
        wave.damping = 0.0;     // Damping: 0 = undamped, try 1.0 for damped mode

        // ---- Time parameters ----
        wave.time_step = 5.0e-4;  // Must satisfy CFL: dt < h/(c*sqrt(2))
        wave.end_time  = 0.5; //1.0;

        // ---- Mesh parameters ----
        wave.initial_refinement   = 5;   // 4^6 = 4096 cells, 4^5 = 1024 cells
        wave.max_refinement_level = 8;

        // ---- Feature toggles ----
        wave.use_amr          = true;   // Adaptive Mesh Refinement
        wave.use_absorbing_bc = false;  // true = Sommerfeld, false = Dirichlet
        wave.track_energy     = true;   // Write energy_log.csv

        // ---- AMR settings ----
        wave.amr_every_n_steps   = 20;
        wave.amr_refine_fraction  = 0.30;
        wave.amr_coarsen_fraction = 0.10;

        // ---- Output ----
        wave.output_every_n_steps = 10;

        // ---- Choose mode ----
        wave.mode = SimulationMode::PEBBLE_IN_POND;
        // wave.mode = SimulationMode::DAMPED_WAVE;
        // wave.mode = SimulationMode::ABSORBING_BC;

        // For PEBBLE_IN_POND, DAMPED_WAVE, ABSORBING_BC:
        wave.run();

        // ---- Convergence study (MMS) ----
        // Uncomment to run instead of wave.run():
        // WaveEquation<2> mms_wave;
        // mms_wave.c = 1.0;
        // mms_wave.run_convergence_study();
    }
    catch (std::exception &exc)
    {
        std::cerr << "\n----------------------------------------------------\n"
                  << "Exception: " << exc.what() << "\nAborting!\n"
                  << "----------------------------------------------------\n";
        return 1;
    }
    catch (...)
    {
        std::cerr << "\n----------------------------------------------------\n"
                  << "Unknown exception! Aborting!\n"
                  << "----------------------------------------------------\n";
        return 1;
    }

    return 0;
}