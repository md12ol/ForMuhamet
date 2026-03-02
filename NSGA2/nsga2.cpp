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
#include <pagmo/utils/multi_objective.hpp>
#include <pagmo/utils/hypervolume.hpp>
#include <pagmo/topologies/ring.hpp>

// Bridge
#include "SdaProblem_NSGA2.hpp"
#include "../rng.h"
#include "../SDA/SDA.h"
#include "../Graph/Graph.h"

using namespace std;
namespace fs = std::filesystem;

std::mt19937 rng;


// params identical to og code
const int NUM_NODES = 256;
const int NUM_STATES = 12;
const int GENERATIONS = 200;
const int RUNS = 5;
const int NUM_CHARS = 2;
const int MAX_RESP_LEN = 2;
const int RUN_SIM = 30;

const int RUN_SIM_Check = RUN_SIM * 2;
const int maxMuts = 2;

int POP_SIZE = 52;
double MUTATION_RATE = static_cast<double>(maxMuts) / (1 + (NUM_STATES * NUM_CHARS) + (NUM_STATES * NUM_CHARS * MAX_RESP_LEN));
double CROSSOVER_RATE = 0.99;
double ETA = 2.0;
int NUM_ISLANDS = 1;




// Mutation:
// og code: maxMuts = 2 (at 73 genes) -> about 2.7%
// MUTATION_RATE now gets calculated by maxMuts divided by vector size
//double MUTATION_RATE = static_cast<double>(maxMuts) / (1 + (NUM_STATES * NUM_CHARS) + (NUM_STATES * NUM_CHARS * MAX_RESP_LEN));
//const double MUTATION_RATE = 0.04;



int main(int argc, char* argv[]) {

    if (argc >= 6) {
        POP_SIZE = std::stoi(argv[1]);
        MUTATION_RATE = std::stod(argv[2]);
        CROSSOVER_RATE = std::stod(argv[3]);
        ETA = std::stod(argv[4]);
        NUM_ISLANDS = std::stoi(argv[5]);
    } else {
        cout << "Warning: No parameters provided. Using default values." << endl;
        //cout << "Usage: ./your_program [PopSize] [MutRate] [CrossRate]" << endl;
        //cout << "Example: ./your_program 100 0.05 0.5" << endl;
        cout << "Usage: ./pagmo_nsga2 [PopSize] [MutRate] [CrossRate] [ETA] [NumIslands]" << endl;
        cout << "Example: ./pagmo_nsga2 52 0.05 0.5 2.0 1" << endl;
    }

    cout << "=== replication of OG code with PAGMO (SGA) ===" << endl;
    cout << "Nodes: " << NUM_NODES << " | States: " << NUM_STATES << endl;
    cout << "Pop: " << POP_SIZE << " | Gens: " << GENERATIONS << endl;
    cout << "Mutation Rate: " << MUTATION_RATE << " (per Gen)" << endl;
    cout << "Crossover Rate: " << CROSSOVER_RATE << "(per Gen)" << endl;

    stringstream folderName_ss;
    folderName_ss << "Output (NSGA2) - " << POP_SIZE << "PS, "
   << GENERATIONS << "Mevs, "
   << MUTATION_RATE * 100 << "%MuR, "
   << CROSSOVER_RATE * 100 << "%CrR, "
    << ETA << "ETA, "
    << NUM_ISLANDS << "Islands, "
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

    string bestParetoFile = "pareto_front_best.csv";
    fs::path bestParetoPath = finalPath / bestParetoFile;

    bool bestFileExists = fs::exists(bestParetoPath);

    std::ofstream bestParetoData(bestParetoPath, std::ios::app);

    if (!bestFileExists) {
        bestParetoData << "Run;Hypervolume;CrowdingDistance;BestParetoLength_Edges;SDA_Genes" << endl;
    }


    vector<SDA> bestSDAs;
    vector<double> sdaEpiAvg;
    sdaEpiAvg.reserve(RUNS);
    vector<double> sdaEpiAvgChamps;
    sdaEpiAvgChamps.reserve(RUNS);
    double epiLenAvg = 0.0;

    for (int run = 1; run <= RUNS; ++run) {

        int currentSeed = 43 + run;
        pagmo::random_device::set_seed(currentSeed);     // set seed for random number generator
        rng.seed(currentSeed);

        stringstream runFile_ss;
        runFile_ss << "run" << run << ".csv";
        string runFile = runFile_ss.str();
        fs::path runPath = finalPath / runFile;



        if (fs::exists(runPath)) {
            continue;
        }


        std::ofstream runData(runPath, std::ios::app);
        runData << "Generation;Hypervolume;CrowdingDistance;ParetoLength_Edges" << endl;

        cout << "\n--- Run " << run << " of " << RUNS << " ---" << endl;

        // 1. we define a problem
        pagmo::problem p{sda_problem_nsga2(NUM_NODES, NUM_STATES, RUN_SIM)};


        // 2. we tell pagmo to use NSGA2
        // nsga2(gen_per_step, cr, param_c, m, param_m)
        pagmo::algorithm algo{pagmo::nsga2(
            1,
            CROSSOVER_RATE,
            ETA,
            MUTATION_RATE,
            ETA
        )};

        unsigned pop_per_island = POP_SIZE / NUM_ISLANDS;
        if (POP_SIZE % NUM_ISLANDS != 0) {
            cout << "WARNING: PopSize (" << POP_SIZE << ") cannot be divided by " << NUM_ISLANDS << " islands!" << endl;
        }

        pagmo::archipelago archi;
        if (NUM_ISLANDS > 1) {
            archi = pagmo::archipelago(pagmo::ring(), (unsigned)NUM_ISLANDS, algo, p, pop_per_island);
        } else {
            // 3. group of island --> 1u means we want 1 island, algo means that we follow the rules of nsga2, p is problem we have
            archi = pagmo::archipelago(1u, algo, p, (unsigned)POP_SIZE);
        }


        int steps = 10; // 10 steps of 1000 gens
        int genPerStep = GENERATIONS / steps;

        // stringstream hyperconverge_ss;
        // hyperconverge_ss << "hypervolumeConvergence_run" << run << ".csv";
        // fs::path hyperPath = finalPath / hyperconverge_ss.str();
        // std::ofstream hyperConvergeData(hyperPath);
        //
        // hyperConvergeData << "Generation_X;Hypervolume_Y;Spread\n";

        for(int i=0; i<steps; ++i) {
            archi.evolve(genPerStep);
            archi.wait_check();
            bool lastStep = false;
            if (i == steps - 1) lastStep = true;

            pagmo::population pop = (*archi.begin()).get_population();

            // all fitness values of the pop and also the genes
            auto fits = pop.get_f();
            auto xs = pop.get_x();

            // let pagmo sort the individuals
            // output--> first element of tuple are all the different fronts (list of lists),
            // second: how many dominate the individual (list of size 52),
            // third: which one does the indivual dominate (list of lists),
            // fourth: rank of all the individuals (list of size 52)
            auto ndf_tuple = pagmo::fast_non_dominated_sorting(fits);

            // take the first element --> list of all the SDAs sorted by their ranking
            const auto &all_fronts = std::get<0>(ndf_tuple);

            // out of this list we take the champions (non dominated solutions)
            const auto &best_front_indices = all_fronts[0];

            cout << "\n>>> Generation " << genPerStep * (i + 1) << " FINISH. Found " << best_front_indices.size() << " optimal trade-offs (Pareto Front)!" << endl;

            runData << genPerStep * (i + 1) << ";";

            vector<vector<double>> champions_fitness;
            vector<vector<double>> genesAll;
            // go through the champs
            bool first = true;
            for (auto idx : best_front_indices) {
                double epiLen = -fits[idx][0];
                double edges = fits[idx][1];
                champions_fitness.push_back(fits[idx]);
                if (lastStep) genesAll.push_back(xs[idx]);


                cout << "Length: " << std::setw(8) << epiLen
                     << " | Edges: " << edges << endl;
                // if (!first) {
                //     runData << ",";
                // }
                // runData << "(" << epiLen << "," << edges << ")";
                // first = false;
            }
            //runData << "}" << endl;
            vector<double> cdValuesFiltered;
            if (champions_fitness.size() > 2) {
                int infinite = 0;
                vector<double> cdValues = pagmo::crowding_distance(champions_fitness);
                for (int j=0; j < cdValues.size(); j++) {
                    if (std::isinf(cdValues[j])) {
                        infinite++;
                    }else{
                        cdValuesFiltered.push_back(cdValues[j]);
                        if (j - infinite == 0) {
                            cout << cdValues[j];
                        }else {
                            cout << ", " << cdValues[j];
                        }
                    }
                }
            }

            vector<double> ref_point = {100000.0, 100000.0};
            pagmo::hypervolume hv(champions_fitness);
            double current_hv = hv.compute(ref_point);

            FitnessStats cdStats = calcStats(cdValuesFiltered);

            if (lastStep) bestParetoData << run << ";" << current_hv << ";" << cdStats.stdDev << ";" << "{";
            //hyperConvergeData << genPerStep * (i + 1) << ";" << current_hv << ";" << cdStats.stdDev << endl;
            runData << current_hv << ";" << cdStats.stdDev << ";" << "{";
            for (int fitPair = 0; fitPair < champions_fitness.size(); fitPair++) {
                if (lastStep) {
                    bestParetoData << "(" << champions_fitness[fitPair][0] << "," << champions_fitness[fitPair][1] << ")";
                    if (fitPair < champions_fitness.size() - 1) bestParetoData << ",";
                }
                runData << "(" << champions_fitness[fitPair][0] << "," << champions_fitness[fitPair][1] << ")";
                if (fitPair < champions_fitness.size() - 1) runData << ",";
            }
            runData << "}" << endl;
            if (lastStep) {
                bestParetoData << "};{(";
                for (int genesSDA = 0; genesSDA < genesAll.size(); genesSDA++) {
                    for (int geneId = 0; geneId < genesAll[genesSDA].size(); geneId++) {
                        bestParetoData << genesAll[genesSDA][geneId];
                        if (geneId < genesAll[genesSDA].size() - 1) bestParetoData << ",";
                    }
                    if (genesSDA < genesAll.size() - 1) bestParetoData << "),(";
                    else bestParetoData << ")";
                }
            }
            if (lastStep) bestParetoData << "}" << endl;


            // if (i == steps - 1){
            //     cout << "\nBest SDAs found so far: " << endl;
            //
            //     stringstream bestParetoFile_ss;
            //     bestParetoFile_ss << "pareto_front_best.csv";
            //     fs::path bestParetoPath = finalPath / bestParetoFile_ss.str();
            //     std::ofstream bestParetoData(bestParetoPath);
            //
            //     bestParetoData << "Run;Hypervolume;";
            //
            //     for (auto idx : best_front_indices) {
            //         double epiLen = fits[idx][0];
            //         double edges = fits[idx][1];
            //         vector<double> &genes = xs[idx];
            //
            //         cout << "Length: " << std::setw(8) << epiLen
            //              << " | Edges: " << std::setw(5) << edges << " | Champions Genes: {";
            //         for (int i = 0; i < genes.size(); ++i) {
            //             cout  << genes[i];
            //             if (i < genes.size() - 1) cout << ", ";
            //         }
            //         cout << "}" << endl;
            //         paretoData << epiLen << ";" << edges << endl;
            //     }
            //     paretoData.close();
            // }
        }
    }
    return 0;
}