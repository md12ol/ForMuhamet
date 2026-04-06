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
#include <pagmo/batch_evaluators/thread_bfe.hpp>
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

thread_local std::mt19937 rng;

// params identical to og code
const int NUM_NODES = 256;
const int NUM_STATES = 12;
const int GENERATIONS = 10000;
const int RUNS = 30;
const int NUM_CHARS = 2;
const int MAX_RESP_LEN = 2;
const int RUN_SIM = 30;
const double ETA = 2.0;
int NUM_ISLANDS = 1;
const bool ISTOTINF_COST = false;

const int RUN_SIM_Check = RUN_SIM * 2;
//const int maxMuts = 2;

const int POP_SIZE = 496;

// Mutation:
// og code: maxMuts = 2 (at 73 genes) -> about 2.7%
// MUTATION_RATE now gets calculated by maxMuts divided by vector size
//double MUTATION_RATE = static_cast<double>(maxMuts) / (1 + (NUM_STATES * NUM_CHARS) + (NUM_STATES * NUM_CHARS * MAX_RESP_LEN));
double MUTATION_RATE = 0.027;
//double MUTATION_RATE = static_cast<double>(maxMuts) / (1 + (NUM_STATES * NUM_CHARS) + (NUM_STATES * NUM_CHARS * MAX_RESP_LEN));
const double CROSSOVER_RATE = 0.99;


int main(int argc, char* argv[]) {

    if (argc >= 3) {
        NUM_ISLANDS = std::stoi(argv[1]);
        MUTATION_RATE = std::stod(argv[2]);

    } else {
        cout << "Warning: No parameters provided. Using default values." << endl;
        cout << "Usage: ./pagmo_nsga2 [PopSize] [MutRate]" << endl;
    }

    cout << "=== PAGMO (NSGA2) ===" << endl;
    cout << "Nodes: " << NUM_NODES << " | States: " << NUM_STATES << endl;
    cout << "Pop: " << POP_SIZE << " | Gens: " << GENERATIONS << endl;

    stringstream folderName_ss;
    if (ISTOTINF_COST) {
        folderName_ss << "Output (NSGA2 - spread vs cost) - " << POP_SIZE << "PS, " << GENERATIONS << "Mevs, "
                      << MUTATION_RATE * 100 << "%MuR, " << CROSSOVER_RATE * 100 << "%CrR, "
                      << ETA << "ETA, " << NUM_ISLANDS << "Islands, " << RUN_SIM << "SEpis, " << NUM_STATES << "ST";
    }else {
        folderName_ss << "Output (NSGA2 - length vs cost) - " << POP_SIZE << "PS, " << GENERATIONS << "Mevs, "
              << MUTATION_RATE * 100 << "%MuR, " << CROSSOVER_RATE * 100 << "%CrR, "
              << ETA << "ETA, " << NUM_ISLANDS << "Islands, " << RUN_SIM << "SEpis, " << NUM_STATES << "ST";
    }


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
        cerr << "Error while creating folder: " << e.what() << endl;
    }

    string bestParetoFile = "pareto_front_best.csv";
    fs::path bestParetoPath = finalPath / bestParetoFile;

    bool fileExists = fs::exists(bestParetoPath);

    std::ofstream bestParetoData(bestParetoPath, std::ios::app);

    if (!fileExists) {
        if (ISTOTINF_COST) {
            bestParetoData << "Run;Hypervolume;CrowdingDistance;BestParetoSpread_Edges;SDA_Genes" << endl;
        } else {
            bestParetoData << "Run;Hypervolume;CrowdingDistance;BestParetoLength_Edges;SDA_Genes" << endl;
        }
    }

    for (int run = 1; run <= RUNS; ++run) {
        int currentSeed = 43 + run;
        pagmo::random_device::set_seed(currentSeed); // set seed for random number generator
        rng.seed(currentSeed);

        stringstream runFile_ss;
        runFile_ss << "run" << run << ".csv";
        fs::path runPath = finalPath / runFile_ss.str();

        if (fs::exists(runPath)) continue;

        std::ofstream runData(runPath, std::ios::app);
        if (ISTOTINF_COST) {
            runData << "Generation;Hypervolume;CrowdingDistance;ParetoSpread_Edges" << endl;
        } else {
            runData << "Generation;Hypervolume;CrowdingDistance;ParetoLength_Edges" << endl;
        }

        cout << "\n--- Run " << run << " of " << RUNS << " ---" << endl;

        // 1. we define a problem
        pagmo::problem p{sda_problem_nsga2(NUM_NODES, NUM_STATES, RUN_SIM, ISTOTINF_COST)};

        // 2. we tell pagmo to use NSGA2
        // nsga2(gen_per_step, cr, param_c, m, param_m)
        pagmo::algorithm algo{pagmo::nsga2(1, CROSSOVER_RATE, ETA, MUTATION_RATE, ETA)};
        algo.extract<pagmo::nsga2>()->set_bfe(pagmo::bfe{});

        pagmo::archipelago archi;
        unsigned pop_per_island = POP_SIZE / NUM_ISLANDS;
        if (NUM_ISLANDS > 1) {
            archi = pagmo::archipelago(pagmo::ring(), (unsigned)NUM_ISLANDS, algo, p, pop_per_island);
        } else {
            // 3. group of island --> 1u means we want 1 island, algo means that we follow the rules of nsga2, p is problem we have
            archi = pagmo::archipelago(1u, algo, p, (unsigned)POP_SIZE);
        }

        int steps = 100; // 10 steps of 1000 gens (angepasst für besseres Logging)
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
            bool lastStep = (i == steps - 1);

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

            // --- HV & CD for rank 0 ---
            vector<vector<double>> champions_fitness;

            // go through the champs
            for (auto idx : best_front_indices) {
                champions_fitness.push_back(fits[idx]);
            }

            vector<double> cdValuesFiltered;
            if (champions_fitness.size() > 2) {
                vector<double> cdValues = pagmo::crowding_distance(champions_fitness);
                for (double val : cdValues) if (!std::isinf(val)) cdValuesFiltered.push_back(val);
            }
            FitnessStats cdStats = calcStats(cdValuesFiltered);

            pagmo::hypervolume hv(champions_fitness);
            double current_hv = hv.compute({100000.0, 100000.0});

            // --- original part ---
            /*
            runData << genPerStep * (i + 1) << ";" << current_hv << ";" << cdStats.stdDev << ";{";
            for (int fitPair = 0; fitPair < champions_fitness.size(); fitPair++) {
                runData << "(" << champions_fitness[fitPair][0] << "," << champions_fitness[fitPair][1] << ")";
                if (fitPair < champions_fitness.size() - 1) runData << ",";
            }
            runData << "}" << endl;
            */

            // ---- update: all ranks ----
            runData << genPerStep * (i + 1) << ";" << current_hv << ";" << cdStats.stdDev << ";{";
            if (lastStep) bestParetoData << run << ";" << current_hv << ";" << cdStats.stdDev << ";" << "{";

            bool firstGlobal = true;
            for (int r = 0; r < all_fronts.size(); ++r) {
                for (auto idx : all_fronts[r]) {
                    if (!firstGlobal) runData << ",";

                    // runX.csv (EpiLen, Edges, Rank)
                    runData << "(" << -fits[idx][0] << "," << fits[idx][1] << "," << r << ")";

                    // bestParetoData (only rank 0)
                    if (lastStep && r == 0) {
                        if (idx != best_front_indices[0]) bestParetoData << ",";
                        bestParetoData << "(" << fits[idx][0] << "," << fits[idx][1] << ")";
                    }
                    firstGlobal = false;
                }
            }
            runData << "}" << endl;
            // --------------------------------------------------

            if (lastStep) {
                bestParetoData << "};{(";
                for (int g = 0; g < best_front_indices.size(); g++) {
                    int idx = best_front_indices[g];
                    for (int j = 0; j < xs[idx].size(); j++) {
                        bestParetoData << xs[idx][j] << (j < xs[idx].size() - 1 ? "," : "");
                    }
                    bestParetoData << (g < best_front_indices.size() - 1 ? "),(" : ")");
                }
                bestParetoData << "}" << endl;
            }
        }
        cout << "   -> Run " << run << " completed! [Mutation: " << MUTATION_RATE << " | Islands: " << NUM_ISLANDS << "]" << endl;
    }
    return 0;
}