#include <iostream>
#include <vector>
#include <cmath>
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <random>

// Pagmo
#include <pagmo/algorithm.hpp>
#include <pagmo/algorithms/sga.hpp>
#include <pagmo/archipelago.hpp>
#include <pagmo/problem.hpp>

// Bridge
#include "SdaEpiProblem.hpp"
#include "rng.h"
#include "SDA/SDA.h"

using namespace std;

std::mt19937 rng;

// params identical to og code
const int NUM_NODES = 256;
const int NUM_STATES = 12;
const int POP_SIZE = 50;
const int GENERATIONS = 10000;
const int RUNS = 5;
const int NUM_CHARS = 2;
const int MAX_RESP_LEN = 2;

const int TOURNAMENT_SIZE = 7;
const double CROSSOVER_RATE = 1.0;
const int maxMuts = 2;

// Mutation:
// og code: maxMuts = 2 (at 73 genes) -> about 2.7%
// MUTATION_RATE now gets calculated by maxMuts divided by vector size
const double MUTATION_RATE = static_cast<double>(maxMuts) / (1 + (NUM_STATES * NUM_CHARS) + (NUM_STATES * NUM_CHARS * MAX_RESP_LEN));
//const double MUTATION_RATE = 0.04;
int main() {
    pagmo::random_device::set_seed(42);     // set seed for random number generator
    rng.seed((42));
    cout << "=== replication of OG code with PAGMO (SGA) ===" << endl;
    cout << "Nodes: " << NUM_NODES << " | States: " << NUM_STATES << endl;
    cout << "Pop: " << POP_SIZE << " | Gens: " << GENERATIONS << endl;
    cout << "Mutation Rate: " << MUTATION_RATE << " (per Gen)" << endl;


    for (int run = 1; run <= RUNS; ++run) {
        cout << "\n--- Run " << run << " of " << RUNS << " ---" << endl;

        // 1. we define a problem
        pagmo::problem p{sda_epi_length_problem(NUM_NODES, NUM_STATES)};

        // 2. we tell pagmo to use SGA (SGA - Simple Genetic Algorithm)
        // sga(gen_per_step, cr, m, param_m, param_c, tournament_size)
        pagmo::algorithm algo{pagmo::sga(
            1,                  // generation per step
            CROSSOVER_RATE,     // 1, we always do crossover
            10.0,               // crossover distribution with 10 as standard (how much we mix)
            MUTATION_RATE,
            10.0,       // mutation distribution also standard (how much do we change the values)
            (unsigned)TOURNAMENT_SIZE
        )};

        // 3. group of island --> 1u means we want 1 island, algo means that we follow the rules of sga, p is problem we have
        pagmo::archipelago archi{1u, algo, p, (unsigned)POP_SIZE};



        int steps = 100; // 10 steps of 1000 gens
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

            //double bestNow = -(*archi.begin()).get_population().champion_f()[0];
            //cout << "Gen " << (i+1)*genPerStep << ": Best Epi-Length = " << bestNow << endl;
        }

        // *archi.begin gives us first island, get pop gives us 50 individs, champion gives us the one with best fitness at 0
        double finalBest = -(*archi.begin()).get_population().champion_f()[0];
        cout << ">>> RUN " << run << " FINISH. BEST: " << finalBest << endl;
        std::vector<double> bestRawSDA = (*archi.begin()).get_population().champion_x();
        SDA sdaBEST(NUM_STATES, NUM_CHARS, MAX_RESP_LEN, NUM_STATES * (NUM_STATES - 1)/2);
        sdaBEST.setGenes(bestRawSDA);
        sdaBEST.print(std::cout);
    }

    return 0;
}
