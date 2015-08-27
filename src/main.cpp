#include <iostream>
#include <fstream>
#include "Data.h"
#include "MbRandom.h"
#include "Settings.h"
#include "Tree.h"



int main (int argc, char* argv[]) {

    // get the user settings
    Settings mySettings(argc, argv);
    
    // instantiate the random number generator
    MbRandom myRandom;
    
    for (int i = 0; i<mySettings.getNumReplicates(); i++)
        {
    
        // generate the tree
        Tree myTree(&myRandom,
                    mySettings.getSpeciationRate(),
                    mySettings.getExtinctionRate(),
                    mySettings.getFossilizationRate(),
                    mySettings.getNumLivingTaxa(),
                    10.0);
            
        myTree.printTree();
            
        std::stringstream ss_tree;
        ss_tree << mySettings.getOutputFileName() << "." << i << ".tree";
        std::string treefile =  ss_tree.str();
        std::fstream treeFileStream;
        treeFileStream.open( treefile.c_str(), std::fstream::out);
        treeFileStream << myTree.getNewick() << std::endl;
        treeFileStream.close();
        
        // simulate data
            
        std::stringstream ss_morph;
        ss_morph << mySettings.getOutputFileName() << ".morph." << i << ".nex";
        std::string morphfile =  ss_morph.str();
        std::fstream morphFileStream;
        morphFileStream.open( morphfile.c_str(), std::fstream::out);
        Data myMorphologicalData(myTree.getNumberOfNodes(),
                                 mySettings.getNumMorphologicalCharacters(),
                                 &myRandom,
                                 &myTree,
                                 1.0);
        myMorphologicalData.printNexus(morphFileStream,true,true);
            
        std::stringstream ss_mol;
        ss_mol << mySettings.getOutputFileName() << ".mol." << i << ".nex";
        std::string molfile =  ss_mol.str();
        std::fstream molFileStream;
        molFileStream.open( molfile.c_str(), std::fstream::out);
        Data myMolecularData(myTree.getNumberOfNodes(),
                             mySettings.getNumMolecularCharacters(),
                             &myRandom,
                             &myTree,
                             1.0,
                             mySettings.getExchangeabilityParameters(),
                             mySettings.getStationaryFrequenciesParameters(),
                             mySettings.getGammaShapeParameter());
        myMolecularData.printNexus(molFileStream,false,true);
        
        }
    
    return 0;
}