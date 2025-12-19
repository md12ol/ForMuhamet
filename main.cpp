#include "main.h"

std::mt19937 rng;

/**
 * Run the program: initialize a population, evolve the population for the specified number of mating events using
 * the specified fitness function, and output results.  This will be repeated runs times.
 * The command line arguments are explained below:
 * 1. Random Number Seed
 * 2. Fitness Function (0 -> Epidemic Length, 1 -> Profile Matching, 2 -> Epidemic Spread, 3 -> Epidemic Severity)
 * 3. Variants (0.0 -> No Variants, (0.0, 1.0] -> Variants with this Variant Probability)
 * 4. Initial Run Number
 * 5. Number of Runs
 * 6. Population Size
 * 7. Number of Generations/Mating Events
 * 8. Tournament Size
 * 9. Crossover Rate
 * 10. Mutation Rate
 * 11. Maximum Number of Mutations (MNM) (The number of mutations will be in the range [1, MNM])
 * 12. Number of Sample Epidemics
 * 13. Number of States in the SDAs
 * If Profile Matching (aka the first argument is set to 1)
 * 14. The Path to the Profile from Working Directory (likely ./)
 * 15. The Profile Number
 * If using Variants (aka the second argument is more than 0.0 and at most 1.0)
 * 14. The Number of 1s in the Initial Variant String
 * 15. The Minimum Number of Edits for Derived Variants, inclusive
 * 16. The Maximum Number of Edits for Derived Variants, inclusive
 * 17. Coupled Immunity and Infectivity (0 -> Uncoupled, 1 -> Coupled)
 * If using Variants and Uncoupled (aka the 17th argument is 0)
 * 18. The Maximum Absolute Difference Between a Parent and Child Variant's Alpha (for Uncoupled)
 * Note: The use of variants using profile matching fitness is not implemented/possible.
 *
 * @param numCmdLineArgs Number of Command Line Arguments
 * @param cmdLineArgs Those Arguments, Explained Above
 * @return Hopefully a Zero :)
 */

int main(int numCmdLineArgs, char *cmdLineArgs[]) {
    string filename;
    fstream runStats, expStats, readMe; // For File Output
    ostringstream oss;

    // Get the command line arguments and place them in their corresponding variables.
    getArgs(cmdLineArgs);

    // Create Directory for Output
    // This is the common part of the folder name:
    auto appendCommon = [&](std::ostringstream& oss) {
        oss << std::setw(3) << std::setfill('0') << popsize << "PS, "
            << std::setw(6) << generations << "Mevs, "
            << tournSize << "TS, "
            << std::setw(3) << int(crossoverRate * 100) << "%CrR, "
            << std::setw(3) << int(mutationRate * 100) << "%MuR, "
            << maxMuts << "MNM, "
            << std::setw(2) << numSampEpis << "SEpis, "
            << std::setw(2) << SDANumStates << "St";
    };
    // This gets the specific parts:
    if (ctrlFitnessFctn == 0) {
        if (newVarProb > 0.0) {
            if (varCoupled) {
                oss << outRoot << "Output - ELVarC w "
                    << std::fixed << std::setprecision(4) << (newVarProb * 100) << "%VarP, ";
                appendCommon(oss);
                oss << ", "
                    << std::setw(2) << initOneBits << "InitB, "
                    << std::setw(2) << minEdits << "-"
                    << std::setw(2) << maxEdits << "Edits/";
                pathToOut = oss.str();
            } else {
                oss << outRoot
                    << "Output - ELVarU w "
                    << std::fixed << std::setprecision(4) << (newVarProb * 100) << "%VarP, ";
                appendCommon(oss);
                oss << ", "
                    << std::setw(2) << initOneBits << "InitB, "
                    << std::setw(2) << minEdits << "-"
                    << std::setw(2) << maxEdits << "Edits, "
                    << std::fixed << std::setprecision(2) << varAlphaDelta << "%DA/";
                pathToOut = oss.str();
            }
        } else {
            oss << outRoot << "Output - EL w ";
            appendCommon(oss);
            oss << "/";
            pathToOut = oss.str();
        }
    } else if (ctrlFitnessFctn == 1) { // TODO: Not Yet Implemented.
        oss << outRoot
            << "Output - PM" << profileNum << " w ";
        appendCommon(oss);
        oss << "/";
        pathToOut = oss.str();
    } else if (ctrlFitnessFctn == 2) {
        if (newVarProb > 0.0) {
            if (varCoupled) {
                oss << outRoot
                    << "Output - ESprVarC w "
                    << std::fixed << std::setprecision(4) << (newVarProb * 100) << "%VarP, ";
                appendCommon(oss);
                oss << ", "
                    << std::setw(2) << initOneBits << "InitB, "
                    << std::setw(2) << minEdits << "-"
                    << std::setw(2) << maxEdits << "Edits/";
                pathToOut = oss.str();
            } else {
                oss << outRoot
                    << "Output - ESprVarU w "
                    << std::fixed << std::setprecision(4) << (newVarProb * 100) << "%VarP, ";
                appendCommon(oss);
                oss << ", "
                    << std::setw(2) << initOneBits << "InitB, "
                    << std::setw(2) << minEdits << "-"
                    << std::setw(2) << maxEdits << "Edits, "
                    << std::fixed << std::setprecision(2) << varAlphaDelta << "%DA/";
                pathToOut = oss.str();
            }
        } else {
            oss << outRoot << "Output - ESpr w ";
            appendCommon(oss);
            oss << "/";
            pathToOut = oss.str();
        }
    } else {
        if (varCoupled) {
            oss << outRoot
                << "Output - ESevVarC w "
                << std::fixed << std::setprecision(4) << (newVarProb * 100) << "%VarP, ";
            appendCommon(oss);
            oss << ", "
                << std::setw(2) << initOneBits << "InitB, "
                << std::setw(2) << minEdits << "-"
                << std::setw(2) << maxEdits << "Edits/";
            pathToOut = oss.str();
        } else {
            oss << outRoot
                << "Output - ESevVarU w "
                << std::fixed << std::setprecision(4) << (newVarProb * 100) << "%VarP, ";
            appendCommon(oss);
            oss << ", "
                << std::setw(2) << initOneBits << "InitB, "
                << std::setw(2) << minEdits << "-"
                << std::setw(2) << maxEdits << "Edits, "
                << std::fixed << std::setprecision(2) << varAlphaDelta << "%DA/";
            pathToOut = oss.str();
        }
    }
    std::filesystem::create_directories(pathToOut);

    // Determine the Location of the Different Output Files
    oss.str(""); // clears the buffer
    oss.clear(); // clears error flags
    oss << pathToOut << "best" << std::setw(2) << std::setfill('0') << initRunNum << ".dat";
    filename = oss.str();
    expStats.open(filename, ios::out);
    oss.str("");
    oss.clear();
    oss << pathToOut << "readme.dat";
    filename = oss.str();
    readMe.open(filename, ios::out);

    // Generate the Readme File
    createReadMe(readMe);
    readMe.close();

    // Okay, Let's Get Started!
    initAlg();
    cmdLineIntro(cout);

    // Actually start the EA
    for (int run = initRunNum; run < initRunNum + runs; run++) {
        oss.str("");
        oss.clear();
        oss << pathToOut << "run" << std::setw(2) << std::setfill('0') << run << ".dat";
        filename = oss.str();
        runStats.open(filename, ios::out);
        if (verbose) cmdLineRun(run, cout);
        initPop(); // Initialization
        report(runStats); // Initial Report
        for (int mev = 1; mev <= generations; mev++) { // Evolution
            matingEvent();
            if (mev % reportEvery == 0) { // Time to report
                if (verbose) {
                    cout << left << setw(5) << run;
                    cout << left << setw(4) << mev / reportEvery;
                }
                report(runStats);
            }
        }
        runStats.close();
        reportBest(expStats);
        cout << "Done run " << run << ".  " << runs - (run - initRunNum + 1) << " more to go. " << endl;
    }
    expStats.close();
    return (0);
}