#include "testing_cli.h"

#include <cstdlib>
#include <iostream>

// Purpose: Runs the program entry point.
// Inputs: argc (int); argv (char **).
// Output: The computed value or status code.
int main(int argc, char **argv) {
    int code = runCli(argc, argv);
    std::cout.flush();
    std::cerr.flush();
    std::_Exit(code);
}
