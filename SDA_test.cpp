#include "SDA/SDA.h"

int main() {
    // Create an SDA object
    SDA sda(10, 5, 3, 40);

    // Create a vector called output that will store the output from an SDA
    vector<int> output(40);

    // Fill the vector output using output from the SDA
    sda.fillOutput(output, true, std::cout);

    return 0;
}

