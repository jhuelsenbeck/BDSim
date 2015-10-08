#include <iostream>
#include <fstream>
#include "Data.h"
#include "FossilCalibration.h"
#include "MbRandom.h"
#include "Node.h"
#include "Settings.h"
#include "Tree.h"

void printCtlFile(std::fstream& ss, int rep, Settings& ms, Tree& t);
void printCtlFileSSBDP(std::fstream& ss, int rep, Settings& ms, Tree& t);



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
        printCtlFileSSBDP(rbCtlStrm, i, mySettings, myTree);
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
    ss << "D1 <- readDiscreteCharacterData(file=\"" << molfile <<   "\");" << std::endl;
    ss << "D2 <- readDiscreteCharacterData(file=\"" << morphfile << "\");" << std::endl;
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

void printCtlFileSSBDP(std::fstream& ss, int rep, Settings& ms, Tree& t) {

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

    ss << "mi = 0;" << std::endl;
    ss << "D_mol <- readDiscreteCharacterData(file=\"" << molfile <<   "\");" << std::endl;
    ss << "D_morph <- readDiscreteCharacterData(file=\"" << morphfile << "\");" << std::endl;
    ss << "starting_tree <- readTrees(\"" << treefile << "\")[1];" << std::endl;
    ss << "taxon_names <- D_mol.names();" << std::endl << std::endl;

    ss << "##############" << std::endl;
    ss << "# Tree model #" << std::endl;
    ss << "##############" << std::endl << std::endl;

    ss << "# Specify a prior on the diversification and turnover rate" << std::endl;
    ss << "diversification ~ dnGamma(5,1);" << std::endl;
    ss << "turnover ~ dnGamma(5,1);" << std::endl << std::endl;

    ss << "# now transform the diversification and turnover rates into speciation and extinction rates" << std::endl;
    ss << "speciation := diversification + turnover;" << std::endl;
    ss << "extinction := turnover;" << std::endl << std::endl;

    ss << "### root age ###" << std::endl;
    ss << "## Uniform prior on root age" << std::endl;
    ss << "root_age ~ dnUnif(0.0, 100.0);" << std::endl;
    ss << "root_age.setValue(20.0);" << std::endl << std::endl;

    ss << "## The sampling-through-time parameter" << std::endl;
    ss << "fossil_rate <- 0.1;" << std::endl << std::endl;

    ss << "## The probability of sampling at the present (rho)" << std::endl;
    ss << "## (this is still in a vector, with a single element as required by the model fxn)" << std::endl;
    ss << "## We use a strongly informative prior with an expected value of rho that is very small" << std::endl;
    ss << "###sampling_prob ~ dnBeta(1.0, 9999.0)" << std::endl;
    ss << "sampling_prob <- 1;" << std::endl << std::endl;

    ss << "## the Birth-death distribution ##" << std::endl;
    ss << "psi ~ dnBDPSerial(rootAge=root_age, lambda=speciation, mu=extinction, psi=fossil_rate, rho=sampling_prob, timeSinceLastSample=0, condition=\"time\", names=taxon_names);" << std::endl << std::endl;

    ss << "## Specify a stable starting tree" << std::endl;
    ss << "psi.setValue(starting_tree);" << std::endl << std::endl;

    ss << "# create some moves that change the stochastic variables" << std::endl;
    ss << "# all moves are sliding proposals but you could use scaling proposals for the rates too" << std::endl;
    ss << "moves[++mi] = mvSlide(diversification,delta=1,tune=true,weight=1);" << std::endl;
    ss << "moves[++mi] = mvSlide(turnover,delta=1,tune=true,weight=1);" << std::endl;
    ss << "moves[++mi] = mvSlide(root_age,delta=1,tune=true,weight=1);" << std::endl;

    ss << "moves[++mi] = mvNarrow(psi, weight=5.0);" << std::endl;
    ss << "moves[++mi] = mvNNI(psi, weight=1.0);" << std::endl;
    ss << "moves[++mi] = mvFNPR(psi, weight=3.0);" << std::endl;
    ss << "moves[++mi] = mvSubtreeScale(psi, weight=3.0);" << std::endl;
    ss << "moves[++mi] = mvNodeTimeSlideUniform(psi, weight=15.0);" << std::endl;

    ss << "################################" << std::endl;
    ss << "# Molecular Substitution Model #" << std::endl;
    ss << "################################" << std::endl << std::endl;

    ss << "# substitution model priors" << std::endl;
    ss << "bf <- v(1,1,1,1);" << std::endl;
    ss << "e <- v(1,1,1,1,1,1);" << std::endl;
    ss << "pi ~ dnDirichlet(bf);" << std::endl;
    ss << "er ~ dnDirichlet(e);" << std::endl << std::endl;

    ss << "# moves on the substitution process parameters" << std::endl;
    ss << "# first some moves only changing one value of the simplex" << std::endl;
    ss << "moves[++mi] = mvSimplexElementScale(pi, alpha=10.0, tune=true, weight=2.0);" << std::endl;
    ss << "moves[++mi] = mvSimplexElementScale(er, alpha=10.0, tune=true, weight=3.0);" << std::endl << std::endl;

    ss << "# the rate matrix" << std::endl;
    ss << "Q_mol := fnGTR(er,pi);" << std::endl << std::endl;

    ss << "####################################" << std::endl;
    ss << "# Morphological Substitution Model #" << std::endl;
    ss << "####################################" << std::endl << std::endl;

    ss << "# the rate matrix" << std::endl;
    ss << "Q_morph <- fnJC(2);" << std::endl << std::endl;

    ss << "#############################" << std::endl;
    ss << "# Among Site Rate Variation #" << std::endl;
    ss << "#############################" << std::endl << std::endl;

    ss << "alpha_prior <- 0.05;" << std::endl;
    ss << "alpha ~ dnExponential( alpha_prior );" << std::endl;
    ss << "gamma_rates := fnDiscretizeGamma( alpha, alpha, 4, false );" << std::endl << std::endl;

    ss << "# add moves for the stationary frequencies, exchangeability rates and the shape parameter" << std::endl;
    ss << "moves[++mi] = mvScale(alpha,weight=2);" << std::endl << std::endl;

    ss << "#############################" << std::endl;
    ss << "# Strict global clock model #" << std::endl;
    ss << "#############################" << std::endl << std::endl;

    ss << "log_mol_clock_rate ~ dnUniform(-10,2);" << std::endl;
    ss << "mol_clock_rate := 10^log_mol_clock_rate;" << std::endl;

    ss << "moves[++mi] = mvSlide(log_mol_clock_rate, weight=1.0);" << std::endl;

    ss << "log_morph_clock_rate ~ dnUniform(-10,2);" << std::endl;
    ss << "morph_clock_rate := 10^log_morph_clock_rate;" << std::endl;

    ss << "moves[++mi] = mvSlide(log_morph_clock_rate, weight=1.0);" << std::endl << std::endl;

    ss << "###################" << std::endl;
    ss << "# PhyloCTMC Model #" << std::endl;
    ss << "###################" << std::endl << std::endl;

    ss << "seq_mol ~ dnPhyloCTMC(tree=psi, Q=Q_mol, branchRates=mol_clock_rate, siteRates=gamma_rates, type=\"DNA\");" << std::endl;
    ss << "seq_morph ~ dnPhyloCTMC(tree=psi, Q=Q_morph, branchRates=morph_clock_rate, type=\"Standard\");" << std::endl << std::endl;

    ss << "# attach the data" << std::endl;
    ss << "seq_mol.clamp(D_mol);" << std::endl;
    ss << "seq_morph.clamp(D_morph);" << std::endl << std::endl;

    ss << "mymodel = model(bf);" << std::endl;

    ss << "monitor_index = 0;" << std::endl;
    ss << "monitors[++monitor_index] = mnModel(filename=\"output/SerialSampled_BDP.log\",printgen=10, separator = TAB);" << std::endl;
    ss << "monitors[++monitor_index] = mnFile(filename=\"output/SerialSampled_BDP.trees\",printgen=10, separator = TAB, psi);" << std::endl;
    ss << "monitors[++monitor_index] = mnScreen(printgen=1000, mol_clock_rate, morph_clock_rate, root_age);" << std::endl;

    ss << "mymcmc = mcmc(mymodel, monitors, moves);" << std::endl;

    ss << "mymcmc.burnin(generations=10000,250);" << std::endl;
    ss << "mymcmc.run(generations=30000);" << std::endl;


    ss << "# Now, we will analyze the tree output." << std::endl;
    ss << "# Let us start by reading in the tree trace" << std::endl;
    ss << "treetrace = readTreeTrace(\"output/SerialSampled_BDP.trees\", treetype=\"clock\");" << std::endl;
    ss << "# and get the summary of the tree trace" << std::endl;
    ss << "#treetrace.summarize()" << std::endl;

    ss << "map_tree = mapTree(treetrace,\"output/SerialSampled_BDP.tree\");" << std::endl << std::endl;

    ss << "# you may want to quit RevBayes now" << std::endl;
    ss << "q();" << std::endl;

}
