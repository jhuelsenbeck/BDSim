#include <iostream>
#include "Settings.h"



int main (int argc, char* argv[]) {

    // get the user settings
    Settings mySettings(argc, argv);

    std::cout << "test" << std::endl;
    return 0;
}