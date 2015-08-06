#include <iostream>
#include "MbRandom.h"
#include "Settings.h"
#include "Tree.h"



int main (int argc, char* argv[]) {

    // get the user settings
    Settings mySettings(argc, argv);
    
    // instantiate the random number generator
    MbRandom myRandom;
    
    // generate the tree
    Tree myTree(&myRandom, mySettings.getSpeciationRate(), mySettings.getExtinctionRate(), mySettings.getNumLivingTaxa(), 10.0);


    return 0;
}