#include "testing_cli.h"

#include <cstdlib>
#include <iostream>

int main(int argc, char **argv) {
    int code = runCli(argc, argv);
    std::cout.flush();
    std::cerr.flush();
    std::_Exit(code);
}
