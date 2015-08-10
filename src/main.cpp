#include <iostream>
#include "Data.h"
#include "MbRandom.h"
#include "Settings.h"
#include "Tree.h"



int main (int argc, char* argv[]) {

    // get the user settings
    Settings mySettings(argc, argv);
    
    // instantiate the random number generator
    MbRandom myRandom;
    
    // generate the tree
    Tree myTree(&myRandom, mySettings.getSpeciationRate(), mySettings.getExtinctionRate(), mySettings.getFossilizationRate(), mySettings.getNumLivingTaxa(), 10.0);
    myTree.printTree();
    std::cout << myTree.getNewick() << std::endl;
    
    // simulate data
    Data myMorphologicalData(myTree.getNumberOfNodes(), 50, &myRandom, &myTree, 1.0);
    std::vector<double> theta = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    std::vector<double> pi    = {0.25, 0.25, 0.25, 0.25};
    Data myMolecularData(myTree.getNumberOfNodes(), 50, &myRandom, &myTree, 1.0, theta, pi);
    
    return 0;
}