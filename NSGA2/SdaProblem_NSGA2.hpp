#ifndef SDA_PROBLEM_NSGA2_HPP
#define SDA_PROBLEM_NSGA2_HPP

#include <pagmo/problem.hpp>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include "../SDA/SDA.h"
#include "../Graph/Graph.h"

struct FitnessStats {
    double mean = 0.0;
    double stdDev = 0.0;
    double ci95 = 0.0;
    double best = 0.0;
    double worst = 0.0;
};

// this function directly take the fitness values and gets us the stats we want to get
// OG code but without biggerbetter since we want to maximize now
FitnessStats calcStats(const std::vector<double> &fits) {
    if (fits.empty()) return {0.0, 0.0, 0.0, 0.0, 0.0};

    double sum = 0.0;
    double bestVal = fits[0];
    double worstVal = fits[0];

    for (double val : fits) {
        sum += val;
        if (val > bestVal) bestVal = val;
        if (val < worstVal) worstVal = val;
    }

    double count = static_cast<double>(fits.size());
    double mean = sum / count;
    double stdDev = 0.0;
    double CI95 = 0.0;

    if (count > 1.0) {
        double stdDevSum = 0.0;
        for (double val : fits) {
            stdDevSum += std::pow(val - mean, 2);
        }
        stdDev = std::sqrt(stdDevSum / (count - 1.0));
        CI95 = 1.96 * (stdDev / std::sqrt(count));
    }

    return {mean, stdDev, CI95, bestVal, worstVal};
}
// ---

// Pagmo
struct sda_epi_length_problem {
    int n_nodes;
    int n_states;
    int run_sim;
    // consts as in og code
    const int NUM_CHARS = 2;
    const int MAX_RESP_LEN = 2;
    const int RUN_SIM = 30;

    sda_epi_length_problem(int nodes = 256, int states = 12, int runs = 5)
        : n_nodes(nodes), n_states(states), run_sim(runs) {}


    std::vector<double> fitness(const std::vector<double> &x) const {
        // first buildSDA
        int outputLen = n_nodes * (n_nodes - 1) / 2;
        SDA sda(n_states, NUM_CHARS, MAX_RESP_LEN, outputLen);

        // add genes
        const_cast<sda_epi_length_problem*>(this)->applyGenesToSDA(sda, x);

        // generate output and count edges
        std::vector<int> weights(outputLen);
        const_cast<SDA&>(sda).fillOutput(weights, false, std::cout);

        long edgeCount = 0;
        for(int w : weights) {
            if(w == 1) edgeCount++;
        }

        // necrotic filter 2.0
        long minEdges = 1 * n_nodes;      // 256
        long maxEdges = 6 * n_nodes;      // 1536

        // if graph is nectrotic, we add a penalty and graph will not be built at all
        if (edgeCount < minEdges || edgeCount > maxEdges) {
            // instead of kicking out the necrotic graphs, we penalize them adding 1000 to the
            // absolut value of the difference between the edgecount and our border --> we give pagmo a direction so it knows that less edges can be better and vice versa
            double penalty = 1000.0 + std::abs(edgeCount - maxEdges);
            return { penalty };
        }

        // now we can start the simulation
        Graph g(n_nodes);
        g.fill(weights, true);

        long epiLenSum = 0;

        // run the simulation run_sim - times and resetting the graph before each simulation
        for (int sim = 0; sim < run_sim; sim++) {
            double alpha = 0.5;
            std::vector<int> epiProfile(2000, 0);
            int totInf = 0;
            int epiLenSingle = g.SIR(0, alpha, epiProfile, totInf);
            epiLenSum = epiLenSum + epiLenSingle;
        }
        double epiLenAvg = static_cast<double>(epiLenSum) / run_sim;
        return { -static_cast<double>(epiLenAvg) };
    }

    // number of genes
    std::vector<double>::size_type get_nx() const {
        // InitChar (1) + Transitions (States*Chars) + Responses (States*Chars*MaxRespLen)
        // 1 + (12*2) + (12*2*2) = 1 + 24 + 48 = 73
        return 1 + (n_states * NUM_CHARS) + (n_states * NUM_CHARS * MAX_RESP_LEN);
    }

    //bounds of genes
    std::pair<std::vector<double>, std::vector<double>> get_bounds() const {
        auto n = get_nx();
        return {std::vector<double>(n, 0.0), std::vector<double>(n, 1.0)};
    }

    // vector x to SDA
    void applyGenesToSDA(SDA &sda, const std::vector<double> &genes) const {
        sda.setGenes(genes);
    }
};

#endif