#include <iostream>
#include <fstream>
#include "Data.h"
#include "FossilCalibration.h"
#include "MbRandom.h"
#include "Settings.h"
#include "Tree.h"

void printCtlFile(std::fstream& ss, int rep, Settings& ms, Tree& t);



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
                    10.0);
        myTree.printTree();
        myTree.listCalibrations();
            
        std::stringstream ss_tree;
        ss_tree << mySettings.getOutputFileName() << "." << i+1 << ".tree";
        std::string treefile =  ss_tree.str();
        std::fstream treeFileStream;
        treeFileStream.open( treefile.c_str(), std::fstream::out);
        treeFileStream << myTree.getNewick() << std::endl;
        treeFileStream.close();
        
        // simulate data
        std::stringstream ss_morph;
        ss_morph << mySettings.getOutputFileName() << ".morph." << i+1 << ".nex";
        std::string morphfile =  ss_morph.str();
        std::fstream morphFileStream;
        morphFileStream.open( morphfile.c_str(), std::fstream::out);
        Data myMorphologicalData(myTree.getNumberOfNodes(),
                                 mySettings.getNumMorphologicalCharacters(),
                                 &myRandom,
                                 &myTree,
                                 1.0);
        myMorphologicalData.printNexus(morphFileStream,true,true);
        morphFileStream.close();
            
        std::stringstream ss_mol;
        ss_mol << mySettings.getOutputFileName() << ".mol." << i+1 << ".nex";
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
        molFileStream.close();
        
        // print RevBayes control file
        std::stringstream rbCtlStrStrm;
        rbCtlStrStrm << mySettings.getOutputFileName() << ".ctl." << i+1 << ".Rev";
        std::string rbCtlfile =  rbCtlStrStrm.str();
        std::fstream rbCtlStrm;
        rbCtlStrm.open( rbCtlfile.c_str(), std::fstream::out);
        printCtlFile(rbCtlStrm, i, mySettings, myTree);
        rbCtlStrm.close();
        }
    
    return 0;
}

void printCtlFile(std::fstream& ss, int rep, Settings& ms, Tree& t) {

    // get file names
    std::stringstream ss_mol;
    ss_mol << ms.getOutputFileName() << ".mol." << rep+1 << ".nex";
    std::string molfile = ss_mol.str();
    std::stringstream ss_morph;
    ss_morph << ms.getOutputFileName() << ".morph." << rep+1 << ".nex";
    std::string morphfile = ss_morph.str();
    std::stringstream ss_tree;
    ss_tree << ms.getOutputFileName() << "." << rep+1 << ".tree";
    std::string treefile = ss_tree.str();

    ss << std::endl;
    ss << "# read data" << std::endl;
    ss << "D1 <- readDiscreteCharacterData(file=\")" << molfile <<   "\");" << std::endl;
    ss << "D2 <- readDiscreteCharacterData(file=\")" << morphfile << "\");" << std::endl;
    ss << "n_mol_sites <- D1.nchar();" << std::endl;
    ss << "n_morph_sites <- D2.nchar();" << std::endl << std::endl;

    ss << "# initialize an iterator for the moves vector " << std::endl;
    ss << "mi = 1;" << std::endl << std::endl;

    ss << "# read the tree " << std::endl;
    ss << "T <- readTrees(\"" << treefile << "\")[1];" << std::endl;
    ss << "n_taxa <- T.ntips();" << std::endl;
    ss << "names <- T.names();" << std::endl << std::endl;

    ss << "# set up birth-death process model " << std::endl;
    ss << "diversification ~ dnExponential(10.0);" << std::endl;
    ss << "moves[mi++] = mvScale(diversification,lambda=1.0,tune=true,weight=3.0);" << std::endl;
    ss << "turnover ~ dnBeta(2.0, 2.0);" << std::endl;
    ss << "moves[mi++] = mvSlide(turnover,delta=1.0,tune=true,weight=3.0);" << std::endl;
    ss << "denom := abs(1.0 - turnover);" << std::endl;
    ss << "birth_rate := diversification / denom;" << std::endl;
    ss << "death_rate := (turnover * diversification) / denom;" << std::endl;
    ss << "rho <- 1.0;" << std::endl;
    ss << "tHesperocyon <- 38.0;" << std::endl;                                     // TEMP
    ss << "mean_ra <- 11.0;" << std::endl;                                          // TEMP
    ss << "stdv_ra <- 0.25;" << std::endl;                                          // TEMP
    ss << "mu_ra <- ln(mean_ra) - ((stdv_ra*stdv_ra) * 0.5);" << std::endl;         // TEMP
    ss << "root_time ~ dnLnorm(mu_ra, stdv_ra, offset=tHesperocyon);" << std::endl; // TEMP

    for (int i=0; i<t.getNumberOfCalibrations(); i++)
        {
        FossilCalibration* fc = t.getCalibrationIndexed(i);
        ss << "clade_" << i+1 << " <- clade(";
        for (int j=0; j<fc->getNumTaxaInClade(); j++)
            {
            ss << "\"" << fc->getCladeTaxonIndexed(j) << "\"";
            if (j+1 != fc->getNumTaxaInClade())
                ss << ",";
            }
        ss << ");" << std::endl;
        }
    ss << "constraints <- v(";
    for (int i=0; i<t.getNumberOfCalibrations(); i++)
        {
        ss << "clade_" << i+1;
        if (i+1 != t.getNumberOfCalibrations())
            ss << ",";
        }
    ss << ");" << std::endl;

    ss << "timetree ~ dnBDP(lambda=birth_rate, mu=death_rate, rho=rho, rootAge=root_time, samplingStrategy=\"uniform\", condition=\"nTaxa\", nTaxa=n_taxa, names=names,constraints=constraints);" << std::endl;
    ss << "timetree.setValue(T);" << std::endl;

#if 0
    ### first create a deterministic node for the age of the MRCA of all Ursidae
    tmrca_Ursidae := tmrca(timetree,clade_Ursidae)
    ### create an additional deterministic node that is negative, representing time in the past
    n_TMRCA_Ursidae := -(tmrca_Ursidae)
    ### the age of the node is a function of the fossil age, which as the observation time of -11.2 My from the present
    tKretzoiarctos <- -11.2
    ### the age of the fossil is a stochastic node that has a lognormal waiting time that is a function of the age of the calibration node
    ### we will model this waiting time so that the fossil is 10 My younger than its ancestor
    M <- 10
    sdv <- 0.25
    mu <- ln(M) - ((sdv * sdv) * 0.5)
    ### using these parameters we create a stochastic node for the fossil 
    ### Kretzoiarctos is a fossil panda, and thus a crown fossil
    crown_Ursid_fossil ~ dnLnorm(mu, sdv, offset=n_TMRCA_Ursidae)
    ### then we clamp this node with the observed fossil specimen, treating the fossil like data
    crown_Ursid_fossil.clamp(tKretzoiarctos)

    ### create a deterministic node for the age of the MRCA of all ursids and the gray seal
    tmrca_UrsidaePinn := tmrca(timetree,clade_UrsPinn)
    n_TMRCA_UrsidaePinn := -(tmrca_UrsidaePinn)
    ### the age of the node is a function of the fossil age, which as the observation time of -33.0 My from the present
    tParictis <- -33.9
    ### the age of the fossil is a stochastic node that has an exponential waiting time that is a function of the age of the calibration node
    ### we will model this waiting time so that the fossil is 30 My younger than its ancestor
    ### Parictis is the earliest known genus of bears, in the family Ursidae
    ### this fossil is a stem ursid fossil, and branched off the tree before the MRCA of living bears
    stem_Ursid_fossil ~ dnExponential(lambda=0.0333, offset=n_TMRCA_UrsidaePinn)
    ### then we clamp this node with the observed fossil specimen, treating the fossil like data
    stem_Ursid_fossil.clamp(tParictis)

    ### add moves on the tree node times, including the root time, which is outside of the timetree 
    moves[mi++] = mvNodeTimeSlideUniform(timetree, weight=30.0)
    moves[mi++] = mvSlide(root_time, delta=2.0, tune=true, weight=10.0)
    moves[mi++] = mvScale(root_time, lambda=2.0, tune=true, weight=10.0)
    moves[mi++] = mvTreeScale(tree=timetree, rootAge=root_time, delta=1.0, tune=true, weight=3.0)
#endif

    ss << "clock_rate ~ dnGamma(2.0,4.0);" << std::endl;
    ss << "moves[mi++] = mvScale(clock_rate,lambda=0.5,tune=true,weight=5.0);" << std::endl;
    ss << "Q <- fnJC(4);" << std::endl;
    ss << "phySeq ~ dnPhyloCTMC(tree=timetree, Q=Q, branchRates=clock_rate, nSites=n_sites, type=\"DNA\");" << std::endl;
    ss << "phySeq.clamp(D1);" << std::endl;
    
#if 0
    mymodel = model(er)

    ### set up the monitors that will output parameter values to file and screen 
    monitors[1] = mnFile(filename="output/GMC_bears_mcmc.log", printgen=10,
                                diversification, turnover, birth_rate, death_rate, 
                                root_time, tmrca_Ursidae, tmrca_UrsidaePinn,
                                er, sf,
                                clock_rate)
    monitors[2] = mnFile(filename="output/GMC_bears_mcmc.trees", printgen=10, timetree)
    monitors[3] = mnScreen(printgen=1000, clock_rate, root_time)

    ### workspace mcmc ###
    mymcmc = mcmc(mymodel, monitors, moves)

    ### pre-burnin to tune the proposals ###
    mymcmc.burnin(generations=1000,tuningInterval=250)

    ### run the MCMC ###
    mymcmc.run(generations=100000)

    ### display proposal acceptance rates and tuning ###
    mymcmc.operatorSummary()

    ### summarize the trees ###
    tt = readTreeTrace("output/GMC_bears_mcmc.trees", "clock")
    tt.summarize()

    ### write MAP tree to file
    mapTree(tt, "output/GMC_bears_mcmc_MAP.tre")

    ## quit ##
    #q()
#endif
}
