#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>

// Bridge zu deinen Modellen
#include "../SDA/SDA.h"
#include "../Graph/Graph.h"
#include "../rng.h" // Den Header mit einbinden

using namespace std;


std::mt19937 rng;

// params for the sda
const int NUM_NODES = 256;
const int NUM_STATES = 12;
const int NUM_CHARS = 2;
const int MAX_RESP_LEN = 2;
const int OUTPUT_LEN = NUM_NODES * (NUM_NODES - 1) / 2;

// --- strings to vector ---
vector<vector<double>> parseNestedTuples(const string& str) {
    vector<vector<double>> result;
    vector<double> current_vec;
    string current_num = "";

    for (char c : str) {
        if (c == '{' || c == '}' || c == ' ') continue;

        if (c == '(') {
            current_vec.clear();
        } else if (c == ')') {
            if (!current_num.empty()) {
                current_vec.push_back(stod(current_num));
                current_num = "";
            }
            if (!current_vec.empty()) {
                result.push_back(current_vec);
            }
        } else if (c == ',') {
            if (!current_num.empty()) {
                current_vec.push_back(stod(current_num));
                current_num = "";
            }
        } else {
            current_num += c;
        }
    }
    return result;
}

int main() {
    cout << "=== NETWORK EXTRACTOR ===" << endl;

    string targetFolder = "../cmake-build-release-wsl/Output/Output (NSGA2 - spread vs cost) - 96PS, 20000Mevs, 2.7%MuR, 99%CrR, 2ETA, 1Islands, 30SEpis, 12ST";
    string inputFile = targetFolder + "/pareto_front_best.csv";
    string outputFile = targetFolder + "/best_run_networks.csv";

    ifstream infile(inputFile);
    if (!infile.is_open()) {
        cerr << "error: " << inputFile << " was not found!" << endl;
        cerr << "Please check dir's name in code" << endl;
        return 1;
    }

    string line;

    getline(infile, line);

    double max_hv = -1.0;
    string best_fitness_str = "";
    string best_genes_str = "";
    int best_run = -1;

    // find run with best hv
    while (getline(infile, line)) {
        if (line.empty()) continue;

        stringstream ss(line);
        string token;
        vector<string> columns;


        while (getline(ss, token, ';')) {
            columns.push_back(token);
        }

        if (columns.size() >= 5) {
            double current_hv = stod(columns[1]);
            if (current_hv > max_hv) {
                max_hv = current_hv;
                best_run = stoi(columns[0]);
                best_fitness_str = columns[3]; // ParetoLength_Edges
                best_genes_str = columns[4];   // SDA_Genes
            }
        }
    }
    infile.close();

    if (best_run == -1) {
        cerr << "no data found..." << endl;
        return 1;
    }

    cout << ">>> Best Run found: Run " << best_run << " (Hypervolume: " << max_hv << ")" << endl;

    // 3. parse strings
    vector<vector<double>> best_fitness = parseNestedTuples(best_fitness_str);
    vector<vector<double>> best_genes = parseNestedTuples(best_genes_str);

    cout << ">>> " << best_genes.size() << " SDAs extracted on pareto front" << endl;
    cout << ">>> build network" << endl;

    // 4. build sdas and look at edges
    ofstream outfile(outputFile);
    outfile << "Spread;Edges;EdgeWeights\n";

    for (size_t i = 0; i < best_genes.size(); ++i) {
        // take fitness values
        double spread = abs(best_fitness[i][0]);
        double edges = best_fitness[i][1];

        // SDA init
        SDA sda(NUM_STATES, NUM_CHARS, MAX_RESP_LEN, OUTPUT_LEN);
        sda.setGenes(best_genes[i]);


        vector<int> weights(OUTPUT_LEN);
        stringstream dummy;
        sda.fillOutput(weights, false, dummy);

        outfile << spread << ";" << edges << ";[";
        for (size_t w = 0; w < weights.size(); ++w) {
            outfile << weights[w];
            if (w < weights.size() - 1) outfile << ",";
        }
        outfile << "]\n";

        if ((i + 1) % 10 == 0) cout << i + 1 << "/" << best_genes.size() << " SDAs done..." << endl;
    }

    outfile.close();
    cout << "=== Done! Datea saved in: best_run_networks.csv ===" << endl;

    return 0;
}