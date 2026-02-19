#include <iostream>
#include <vector>
#include <cmath>
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <random>
#include <fstream>
#include <filesystem>
#include <sstream>

// Pagmo
#include <pagmo/algorithm.hpp>
#include <pagmo/algorithms/sga.hpp>
#include <pagmo/algorithms/nsga2.hpp>
#include <pagmo/archipelago.hpp>
#include <pagmo/problem.hpp>

// Bridge
#include "SdaProblem_NSGA2.hpp"
#include "rng.h"
#include "SDA/SDA.h"
#include "Graph/Graph.h"

using namespace std;
namespace fs = std::filesystem;

std::mt19937 rng;


// params identical to og code
const int NUM_NODES = 256;
const int NUM_STATES = 12;
const int POP_SIZE = 52;
const int GENERATIONS = 100;
const int RUNS = 3;
const int NUM_CHARS = 2;
const int MAX_RESP_LEN = 2;
const int RUN_SIM = 30;

const int RUN_SIM_Check = RUN_SIM * 2;

const int TOURNAMENT_SIZE = 7;
const double CROSSOVER_RATE = 1.0;
const int maxMuts = 2;


// Mutation:
// og code: maxMuts = 2 (at 73 genes) -> about 2.7%
// MUTATION_RATE now gets calculated by maxMuts divided by vector size
const double MUTATION_RATE = static_cast<double>(maxMuts) / (1 + (NUM_STATES * NUM_CHARS) + (NUM_STATES * NUM_CHARS * MAX_RESP_LEN));
//const double MUTATION_RATE = 0.04;



int main() {
    cout << "=== replication of OG code with PAGMO (SGA) ===" << endl;
    cout << "Nodes: " << NUM_NODES << " | States: " << NUM_STATES << endl;
    cout << "Pop: " << POP_SIZE << " | Gens: " << GENERATIONS << endl;
    cout << "Mutation Rate: " << MUTATION_RATE << " (per Gen)" << endl;

    stringstream folderName_ss;
    folderName_ss << "Output - " << POP_SIZE << "PS, "
   << GENERATIONS << "Mevs, "
   << TOURNAMENT_SIZE << "TS, "
   << MUTATION_RATE * 100 << "%MuR, "
   << CROSSOVER_RATE * 100 << "%CrR, "
    << RUN_SIM << "SEpis, "
   << NUM_STATES << "ST";

    string outputFolder = folderName_ss.str();

    fs::path parentFolder = "Output";
    fs::path finalPath = parentFolder / outputFolder;

    try {
        if (fs::create_directories(finalPath)) {
            cout << "New folder: " << finalPath << std::endl;
        } else {
            cout << "folder already exists: " << finalPath << std::endl;
        }
    } catch (const exception& e) {
        cerr << "Error while creating folder: " << e.what() << std::endl;
    }

    vector<SDA> bestSDAs;
    vector<double> sdaEpiAvg;
    sdaEpiAvg.reserve(RUNS);
    vector<double> sdaEpiAvgChamps;
    sdaEpiAvgChamps.reserve(RUNS);
    double epiLenAvg = 0.0;

    for (int run = 1; run <= RUNS; ++run) {
        stringstream runFile_ss;
        runFile_ss << "run" << run << ".csv";
        string runFile = runFile_ss.str();
        fs::path runPath = finalPath / runFile;



        if (fs::exists(runPath)) {
            continue;
        }
        int currentSeed = 42 + run - 1;
        pagmo::random_device::set_seed(currentSeed);     // set seed for random number generator
        rng.seed(currentSeed);

        std::ofstream runData(runPath, std::ios::app);

        runData << "Generation;BestFitness;WorstFitness;MeanFitness;StdDev;CI95\n";





        cout << "\n--- Run " << run << " of " << RUNS << " ---" << endl;

        // 1. we define a problem
        pagmo::problem p{sda_epi_length_problem(NUM_NODES, NUM_STATES, RUN_SIM)};

        // 2. we tell pagmo to use SGA (SGA - Simple Genetic Algorithm)
        // sga(gen_per_step, cr, m, param_m, param_c, tournament_size)
//        pagmo::algorithm algo{pagmo::sade(1)};
        pagmo::algorithm algo{pagmo::sga(
                    1,                          // gen (Generationen pro Step)
                    CROSSOVER_RATE,             // cr
                    2.0,                       // eta_c (Distribution Index Crossover)
                    MUTATION_RATE,              // m
                    1.0,                       // param_m (Distribution Index Mutation) - je kleiner, desto "wilder"
                    (unsigned)TOURNAMENT_SIZE,  // param_s (Tournament Size)
                    "sbx",              // crossover strategy (Standard)
                    "polynomial",               // mutation strategy (Standard)
                    "tournament",               // selection strategy (Standard)
                    1u                          // <--- n_elite: HIER IST DER SCHLÜSSEL! (1 Champion überlebt)
                )};



        // 3. group of island --> 1u means we want 1 island, algo means that we follow the rules of sga, p is problem we have
        pagmo::archipelago archi{1u, algo, p, (unsigned)POP_SIZE};



        int steps = 10; // 10 steps of 1000 gens
        int genPerStep = GENERATIONS / steps;

        for(int i=0; i<steps; ++i) {
            archi.evolve(genPerStep);
            archi.wait_check();

            // 1. get the population
            pagmo::population pop = (*archi.begin()).get_population();

            // 2. getting all the fitness values
            // Pagmo gives us vector<vector<double>> (since it could be mo)
            auto pagmo_fits = pop.get_f();

            // 3.transform into list (vector<double>) --> getting ready for stats
            std::vector<double> currentFits;
            currentFits.reserve(pagmo_fits.size());

            for(const auto& f : pagmo_fits) {
                currentFits.push_back( -f[0] );
            }

            FitnessStats stats = calcStats(currentFits);

            cout << "Gen " << (i+1)*genPerStep
                 << ": Best=" << stats.best
                 << ", Worst=" << stats.worst
                 << ", Mean=" << stats.mean
                 << ", StdDev=" << stats.stdDev
                 << ", CI95=" << stats.ci95 << endl;


            runData << (1 + i) * genPerStep << ";"
            << stats.best << ";"
            << stats.worst << ";"
            << stats.mean << ";"
            << stats.stdDev << ";"
            << stats.ci95 << "\n";

            runData.flush();

            //double bestNow = -(*archi.begin()).get_population().champion_f()[0];
            //cout << "Gen " << (i+1)*genPerStep << ": Best Epi-Length = " << bestNow << endl;
        }
        string bestRun = "best01.csv";
        fs::path bestRunPath = finalPath / bestRun;
        std::ofstream bestRunData(bestRunPath, std::ios::app);

        // *archi.begin gives us first island, get pop gives us 50 individs, champion gives us the one with best fitness at 0
        double finalBest = -(*archi.begin()).get_population().champion_f()[0];
        cout << ">>> RUN " << run << " FINISH. BEST: " << finalBest << endl;
        std::vector<double> bestRawSDA = (*archi.begin()).get_population().champion_x();

        bestRunData << run << ". Run -> Best Fitness: " << finalBest
        << "\n>>> Champion Genes: { ";
        for (size_t i = 0; i < bestRawSDA.size(); ++i) {
            bestRunData << bestRawSDA[i];
            if (i < bestRawSDA.size() - 1) bestRunData << ", ";
        }
        bestRunData << " }" << endl;

        cout << ">>> Champion Genes: { ";
        for (size_t i = 0; i < bestRawSDA.size(); ++i) {
            cout << bestRawSDA[i];
            if (i < bestRawSDA.size() - 1) cout << ", ";
        }
        cout << " }" << endl;

        SDA sdaBEST(NUM_STATES, NUM_CHARS, MAX_RESP_LEN, NUM_NODES * (NUM_NODES - 1)/2);
        sdaBEST.setGenes(bestRawSDA);
        sdaBEST.print(std::cout);
        bestRunData << "SDA:\n";
        sdaBEST.print(bestRunData);
        std::vector<int> weights(NUM_NODES * (NUM_NODES - 1) / 2);
        const_cast<SDA&>(sdaBEST).fillOutput(weights, false, std::cout);
        Graph g(NUM_NODES);
        g.fill(weights, true);
        bestRunData << "Graph:\n";
        g.print(bestRunData);
        bestRunData << "\n";

        if (finalBest > 0) {
            sdaEpiAvg.push_back(finalBest);
            bestSDAs.reserve(RUNS);
            bestSDAs.push_back(sdaBEST);
        } else {
            cout << ">>> SKIPPING RESULT (Penalty/Necrotic)" << endl;
        }

    }
    cout << "\n" << std::left
             << std::setw(4) << "Run" << " | "
             << std::setw(13) << "Best Length" << " | "
             << "Average after " << RUN_SIM_Check << " Simulations" << endl;
    cout << std::string(60, '-') << endl;

    for (int sdaNr = 0; sdaNr < bestSDAs.size(); ++sdaNr) {
        std::vector<int> weights(NUM_NODES * (NUM_NODES - 1) / 2);
        const_cast<SDA&>(bestSDAs[sdaNr]).fillOutput(weights, false, std::cout);
        Graph g(NUM_NODES);
        g.fill(weights, true);

        long epiLenSum = 0;
        // run the simulation RUN_SIM - times and resetting the graph before each simulation
        for (int sim = 0; sim < RUN_SIM_Check; sim++) {
            double alpha = 0.5;
            vector<int> epiProfile(2000, 0);
            int totInf = 0;
            int epiLenSingle = g.SIR(0, alpha, epiProfile, totInf);
            epiLenSum = epiLenSum + epiLenSingle;
        }
        epiLenAvg = static_cast<double>(epiLenSum) / RUN_SIM_Check;
        sdaEpiAvgChamps.push_back(epiLenAvg);
        cout << std::left
                     << std::setw(4) << (sdaNr + 1) << " | "
                     << std::setw(13) << std::fixed << std::setprecision(2) << sdaEpiAvg[sdaNr] << " | "
                     << std::fixed << epiLenAvg << endl;
    }
    cout << "\nValid Runs: " << bestSDAs.size() << " / " << RUNS << endl;

    if (!bestSDAs.empty()) {
        FitnessStats champStats = calcStats(sdaEpiAvgChamps);
        cout << "\nStats of best SDAs in " << bestSDAs.size() << " valid runs." << endl;
        cout << "\nAverage mean: " << champStats.mean;
        cout << "\nStdDev: " << champStats.stdDev;
        cout << "\nCI95: " << champStats.ci95;
    }
    return 0;
}
