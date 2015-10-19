#include <iostream>
#include <fstream>
#include "Data.h"
#include "FossilCalibration.h"
#include "MbRandom.h"
#include "Node.h"
#include "RevScript.h"
#include "Settings.h"
#include "Tree.h"



int main (int argc, char* argv[]) {

    // get the user settings
    Settings mySettings(argc, argv);
    
    // instantiate the random number generator
    MbRandom myRandom(11);
    
    for (int i = 0; i<mySettings.getNumReplicates(); i++)
        {
        // generate the tree
        Tree myTree(&myRandom,
                    mySettings.getSpeciationRate(),
                    mySettings.getExtinctionRate(),
                    mySettings.getFossilizationRate(),
                    mySettings.getNumLivingTaxa(),
                    mySettings.getSimDuration());
        int numFossils = myTree.getNumFossils();
        //myTree.printTree();
        //myTree.listCalibrations();
            
        std::stringstream ss_tree;
        ss_tree << mySettings.getFullFileName() << "." << i+1 << ".tree";
        std::string treefile =  ss_tree.str();
        std::fstream treeFileStream;
        treeFileStream.open( treefile.c_str(), std::fstream::out);
        treeFileStream << "#NEXUS" << std::endl << std::endl;

        treeFileStream << "begin taxa;" << std::endl;
        int n = 0;
        for (int j=0; j<myTree.getNumberOfDownPassNodes(); j++)
            {
            Node* p = myTree.getDownPassNode(j);
            if (myTree.isExtantTaxon(p) == true || myTree.isFossilTaxon(p) == true)
                n++;
            }
        treeFileStream << "   dimensions ntax=" << n << ";" << std::endl;
        treeFileStream << "   taxlabels" << std::endl;
        for (int j=0; j<myTree.getNumberOfDownPassNodes(); j++)
            {
            Node* p = myTree.getDownPassNode(j);
            if (myTree.isExtantTaxon(p) == true || myTree.isFossilTaxon(p) == true)
                {
                treeFileStream << "   " << p->getName() << std::endl;
                }
            }
        treeFileStream << "   ;" << std::endl;
        treeFileStream << "end;" << std::endl;
        
        treeFileStream << "begin trees;" << std::endl;
        treeFileStream << "   tree true_tree = [&R] " << myTree.getNewick() << ";" << std::endl;
        treeFileStream << "end;" << std::endl << std::endl;
        treeFileStream.close();

        std::stringstream snf_tree;
        snf_tree << mySettings.getFullFileName() << "." << i+1 << ".nfossils";
        std::string snffile =  snf_tree.str();
        std::fstream nfFileStream;
        nfFileStream.open( snf_tree.str().c_str(), std::fstream::out);
        nfFileStream << numFossils;
        nfFileStream << std::endl;
        nfFileStream.close();
        
        // simulate data
        std::stringstream ss_morph;
        ss_morph << mySettings.getFullFileName() << ".morph." << i+1 << ".nex";
        std::string morphfile =  ss_morph.str();
        std::fstream morphFileStream;
        morphFileStream.open( morphfile.c_str(), std::fstream::out);
        Data myMorphologicalData(myTree.getNumberOfNodes(),
                                 mySettings.getNumMorphologicalCharacters(),
                                 &myRandom,
                                 &mySettings,
                                 &myTree,
                                 mySettings.getMorphologicalRate());
        myMorphologicalData.printNexus(morphFileStream,true,true);
        morphFileStream.close();
            
        Data myMolecularData(myTree.getNumberOfNodes(),
                             mySettings.getNumMolecularCharacters(),
                             &myRandom,
                             &mySettings,
                             &myTree,
                             mySettings.getMolecularRate(),
                             mySettings.getExchangeabilityParameters(),
                             mySettings.getStationaryFrequenciesParameters(),
                             mySettings.getGammaShapeParameter());
            
        std::stringstream ss_mol;
        ss_mol << mySettings.getFullFileName() << ".mol." << i+1 << ".nex";
        std::string molfile =  ss_mol.str();
        std::fstream molFileStream;
        molFileStream.open( molfile.c_str(), std::fstream::out);
        myMolecularData.printNexus(molFileStream,false,true);
        molFileStream.close();
            
        std::stringstream ss_mol_extant;
        ss_mol_extant << mySettings.getFullFileName() << ".mol.extant.only." << i+1 << ".nex";
        std::string molfile_extant =  ss_mol_extant.str();
        std::fstream molFileStream_extant;
        molFileStream_extant.open( molfile_extant.c_str(), std::fstream::out);
        myMolecularData.printNexus(molFileStream_extant,false,false);
        molFileStream_extant.close();
        
        // print RevBayes control file
        RevScript script_fbdp = RevScript( mySettings.getFilePath(), mySettings.getFileName(), i+1, &myTree, RevScript::FBDP, RevScript::FOSSIL, false );
        script_fbdp.print();
            
        RevScript script_fbdp_mol_only = RevScript( mySettings.getFilePath(), mySettings.getFileName(), i+1, &myTree, RevScript::FBDP, RevScript::FOSSIL, true );
        script_fbdp_mol_only.print();
            
        RevScript script_bdp = RevScript( mySettings.getFilePath(), mySettings.getFileName(), i+1, &myTree, RevScript::BDP, RevScript::FOSSIL, false );
        script_bdp.print();
        
        RevScript script_bdp_perfect = RevScript( mySettings.getFilePath(), mySettings.getFileName(), i+1, &myTree, RevScript::BDP, RevScript::NODE, false );
        script_bdp_perfect.print();
            
//        std::stringstream rbCtlStrStrm;
//        rbCtlStrStrm << mySettings.getOutputFileName() << ".ctl." << i+1 << ".Rev";
//        std::string rbCtlfile =  rbCtlStrStrm.str();
//        std::fstream rbCtlStrm;
//        rbCtlStrm.open( rbCtlfile.c_str(), std::fstream::out);
//        printCtlFileSSBDP(rbCtlStrm, i, mySettings, myTree);
//        rbCtlStrm.close();
            
            
            
            
        myTree.removeFossils();
        std::stringstream ss_tree_extant;
        ss_tree_extant << mySettings.getFullFileName() << ".extant.only." << i+1 << ".tree";
        std::string treefile_extant =  ss_tree_extant.str();
        std::fstream treeFileStream_extant;
        treeFileStream_extant.open( treefile_extant.c_str(), std::fstream::out);
        treeFileStream_extant << "#NEXUS" << std::endl << std::endl;
            
        treeFileStream_extant << "begin taxa;" << std::endl;
        n = 0;
        for (int j=0; j<myTree.getNumberOfDownPassNodes(); j++)
            {
            Node* p = myTree.getDownPassNode(j);
            if (myTree.isExtantTaxon(p) == true || myTree.isFossilTaxon(p) == true)
                n++;
            }
        treeFileStream_extant << "   dimensions ntax=" << n << ";" << std::endl;
        treeFileStream_extant << "   taxlabels" << std::endl;
        for (int j=0; j<myTree.getNumberOfDownPassNodes(); j++)
            {
            Node* p = myTree.getDownPassNode(j);
            if (myTree.isExtantTaxon(p) == true || myTree.isFossilTaxon(p) == true)
                {
                treeFileStream_extant << "   " << p->getName() << std::endl;
                }
            }
        treeFileStream_extant << "   ;" << std::endl;
        treeFileStream_extant << "end;" << std::endl;
            
        treeFileStream_extant << "begin trees;" << std::endl;
        treeFileStream_extant << "   tree true_tree = [&R] " << myTree.getNewick() << ";" << std::endl;
        treeFileStream_extant << "end;" << std::endl << std::endl;
        treeFileStream_extant.close();
        
        }
    
    return 0;
}

