#include "WaveEquationParallel.hpp"

int main(int argc, char *argv[])
{
    Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);
    MPI_Comm mpi_comm = MPI_COMM_WORLD;

    try
    {
        // Parsing dei parametri da riga di comando
        std::string mode_str = "mms";       // Default: studio di convergenza MMS
        unsigned int ref = 4;               // Livello di raffinamento iniziale baseline
        
        for (int i = 1; i < argc; ++i)
        {
            if (std::string(argv[i]) == "--mode" && i + 1 < argc) mode_str = argv[++i];
            if (std::string(argv[i]) == "--ref" && i + 1 < argc)  ref = std::stoi(argv[++i]);
        }

        WaveEquation<2> wave(mpi_comm);
        wave.initial_refinement = ref;

        // =========================================================================
        // 1. MODALITÀ DI BENCHMARK: STRONG / WEAK SCALING
        // =========================================================================
        if (mode_str == "strong" || mode_str == "weak")
        {
            wave.time_scheme = TimeScheme::NEWMARK;
            wave.time_step = 5.0e-4;
            wave.mode = SimulationMode::PEBBLE_IN_POND;
            wave.use_amr = false; // AMR spento per avere test di scaling puliti

            wave.prepare_for_analysis(); 
            auto result = wave.run_scaling_benchmark(100);
            
            // Salva un file CSV dedicato in base alla modalità scelta
            wave.print_scaling_table({result}, mode_str + "_scaling.csv");
        }
        // =========================================================================
        // 2. MODALITÀ ANALISI DI DISPERSIONE NUMERICA
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
        // 3. MODALITÀ VERIFICA MATEMATICA (MMS - Da lanciare con 1 solo processo)
        // =========================================================================
        else if (mode_str == "mms")
        {
            wave.time_scheme = TimeScheme::NEWMARK;
            wave.run_convergence_study();
        }
        // =========================================================================
        // 4. SCENARIO FISICO: PEBBLE IN POND (Con supporto AMR)
        // =========================================================================
        else if (mode_str == "pebble")
        {
            wave.mode                = SimulationMode::PEBBLE_IN_POND;
            wave.time_scheme         = TimeScheme::LEAPFROG; // Molto veloce per gli scenari fisici
            wave.time_step           = 1.0e-3;
            wave.end_time            = 1.0;
            
            // Configurazione modulo AMR (p4est)
            wave.use_amr             = true;
            wave.amr_every_n_steps   = 5;
            wave.max_refinement_level = 7; 
            
            wave.output_every_n_steps = 10; // Genera un file .pvtu ogni 10 step
            wave.track_energy         = true;

            wave.run();
        }
        // =========================================================================
        // 5. SCENARIO FISICO: INTERFERENCE (Con supporto AMR)
        // =========================================================================

        else if (mode_str == "interference")
        {
            wave.mode                = SimulationMode::INTERFERENCE;
            wave.time_scheme         = TimeScheme::LEAPFROG;
            wave.time_step           = 1.0e-3;
            wave.end_time            = 1.0;
            
            // Configurazione modulo AMR (p4est)
            wave.use_amr             = true;
            wave.amr_every_n_steps   = 5;
            wave.max_refinement_level = 7;
            
            wave.output_every_n_steps = 10;
            wave.track_energy         = true;

            wave.run();
        }
        // =========================================================================
        // 6. SCENARIO FISICO: REFRACTION (Eterogeneità del mezzo)
        // =========================================================================
        else if (mode_str == "refraction")
        {
            wave.mode                = SimulationMode::REFRACTION;
            wave.time_scheme         = TimeScheme::LEAPFROG;
            // Per la rifrazione, la CFL è dettata dalla c_fast (che è più grande di c).
            // Usiamo un dt più piccolo per essere sicuri di non far esplodere il Leapfrog.
            wave.time_step           = 5.0e-4; 
            wave.end_time            = 1.0;
            
            wave.use_amr             = true;
            wave.amr_every_n_steps   = 5;
            wave.max_refinement_level = 7;
            
            wave.output_every_n_steps = 20; // Salviamo un po' meno spesso per non riempire il disco
            wave.track_energy         = true;

            wave.run();
        }
        // =========================================================================
        // 7. SCENARIO FISICO: DIFFRACTION (Ostacolo con fessura)
        // =========================================================================
        else if (mode_str == "diffraction")
        {
            wave.mode                = SimulationMode::DIFFRACTION;
            wave.time_scheme         = TimeScheme::LEAPFROG;
            wave.time_step           = 5.0e-4; // dt prudenziale per i bordi rigidi della fessura
            wave.end_time            = 1.0;
            
            wave.use_amr             = true;
            wave.amr_every_n_steps   = 5;
            wave.max_refinement_level = 7;
            
            wave.output_every_n_steps = 20; 
            wave.track_energy         = true;

            wave.run();
        }
        // =========================================================================
        // 8. SCENARIO FISICO: DAMPING (Sponge Layers / Absorbing BC)
        // =========================================================================
        else if (mode_str == "damping")
        {
            // Assicurati di avere DAMPING (o ABC) definito nel tuo enum SimulationMode
            wave.mode                = SimulationMode::DAMPED_WAVE;
            wave.time_scheme         = TimeScheme::LEAPFROG;
            wave.time_step           = 1.0e-3;
            wave.end_time            = 1.5; // Più lungo per dare tempo all'onda di uscire dallo schermo
            
            wave.use_amr             = true;
            wave.amr_every_n_steps   = 5;
            wave.max_refinement_level = 7;
            
            wave.output_every_n_steps = 10; 
            wave.track_energy         = true; // CRITICO: ci serve per dimostrare l'assorbimento

            wave.run();
        }

        /// =========================================================================
        // 9. SCENARIO FISICO: ABSORBING BC (Sommerfeld Boundary Conditions)
        // =========================================================================
        else if (mode_str == "absorbing")
        {
            wave.mode                = SimulationMode::ABSORBING_BC;
            wave.use_absorbing_bc    = true; // Usiamo il flag booleano che hai nel tuo .hpp!
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
        // MODALITÀ SCONOSCIUTA
        // =========================================================================
        else
        {
            if (Utilities::MPI::this_mpi_process(mpi_comm) == 0)
            {
                std::cout << "Modalita non riconosciuta!\n"
                          << "Scegli tra: strong, weak, dispersion, mms, pebble, interference, refraction, diffraction, damping, absorbing\n";
            }
        }
    }
    
    catch (std::exception &exc)
    {
        std::cerr << "Errore catturato nel main: " << exc.what() << "\n";
        return 1;
    }

    return 0;
}