#include "FossilCalibration.h"
#include "RevScript.h"
#include "Tree.h"

#include <sstream>
#include <fstream>


RevScript::RevScript(const std::string &pa, const std::string &f, int r, Tree *t, TREE_PRIOR p, bool mo) :
    base_file_name( f ),
    base_file_path( pa ),
    rep( r ),
    tree( t ),
    prior( p ),
    clock( STRICT ),
    mol_only( mo )
{
    
    
}



void RevScript::print( void )
{
    
    std::stringstream rbCtlStrStrm;
    rbCtlStrStrm << base_file_path << "/" << base_file_name << ".";
    
    if ( prior == BDP )
    {
        rbCtlStrStrm << "BDP.";
    }
    else
    {
        rbCtlStrStrm << "FBDP.";
    }
    
    if ( clock == STRICT )
    {
        rbCtlStrStrm << "strict.";
    }
    else
    {
        rbCtlStrStrm << "UCLN.";
    }
    
    if ( mol_only == true )
    {
        rbCtlStrStrm << "mol_only.";
    }
    
    rbCtlStrStrm << rep << ".Rev";
    std::string rbCtlfile =  rbCtlStrStrm.str();
    
    std::fstream rbCtlStrm;
    rbCtlStrm.open( rbCtlfile.c_str(), std::fstream::out);

    printInitialSetting( rbCtlStrm );
    
    // read in the molecular data
    printReadMolecularData( rbCtlStrm );
    
    // read in the morphological data
    if ( mol_only == false )
    {
        printReadMorphologicalData( rbCtlStrm );
    }
    
    // read in the true tree
    printReadTree( rbCtlStrm );
    
    // print the tree prior
    if ( prior == BDP )
    {
        printBDP( rbCtlStrm );
    }
    else
    {
        printFBDP( rbCtlStrm );
    }
    
    // print the tree prior
    if ( clock == STRICT )
    {
        printStrictMolecularClock( rbCtlStrm );
    }
    else
    {
        printUclnMolecularClock( rbCtlStrm );
    }
    
    // molecular substitution model
    printMolecularSubstitutionModel( rbCtlStrm );
    
    // morphological substitution model
    if ( mol_only == false )
    {
        printRelativeMorphologicalClock( rbCtlStrm );
        printMorphologicalSubstitutionModel( rbCtlStrm );
    }
    
    // Model
    printModel( rbCtlStrm );

    // Monitors
    printMonitors( rbCtlStrm );
    
    // MCMC
    printMcmc( rbCtlStrm );
    
    // summarize trees
    //printSummarizeTrees( rbCtlStrm );
    
    // quit
    printQuit( rbCtlStrm );
    
    rbCtlStrm.close();

}


void RevScript::printBDP( std::ostream &out )
{
    
    
//    ss << "# read the tree " << std::endl;
//    ss << "T <- readTrees(\"" << treefile << "\")[1]" << std::endl;
//    ss << "n_taxa <- T.ntips()" << std::endl;
//    ss << "names <- T.names()" << std::endl << std::endl;
    
    out << "#######################" << std::endl;
    out << "# birth-death process #" << std::endl;
    out << "#######################" << std::endl;
    out << std::endl;
    out << "speciation ~ dnExponential(1.0)" << std::endl;
    out << "relative_extinction ~ dnBeta(2.0, 2.0)" << std::endl;
    out << "extinction := speciation * relative_extinction" << std::endl;
    out << "sampling_prob <- 1.0" << std::endl;
    out << std::endl;
    out << "moves[++mi] = mvScale(speciation,lambda=1.0,tune=true,weight=3.0)" << std::endl;
    out << "moves[++mi] = mvSlide(relative_extinction,delta=1.0,tune=true,weight=3.0)" << std::endl;
    out << std::endl;
    out << "root_age ~ dnUniform(0, 1000)" << std::endl;
    
    for (int i=0; i<tree->getNumberOfCalibrations(); i++)
    {
        FossilCalibration* fc = tree->getCalibrationIndexed(i);
        out << "clade_" << i+1 << " <- clade(";
        for (int j=0; j<fc->getNumTaxaInClade(); j++)
        {
            out << "\"" << fc->getCladeTaxonIndexed(j) << "\"";
            if (j+1 != fc->getNumTaxaInClade())
                out << ",";
        }
        out << ")" << std::endl;
    }
    out << "constraints <- v(";
    for (int i=0; i<tree->getNumberOfCalibrations(); i++)
    {
        out << "clade_" << i+1;
        if (i+1 != tree->getNumberOfCalibrations())
            out << ",";
    }
    out << ")" << std::endl;
    
    out << "psi ~ dnBDP(lambda=speciation, mu=extinction, rho=sampling_prob, rootAge=root_age, samplingStrategy=\"uniform\", condition=\"nTaxa\", nTaxa=n_taxa, names=taxon_names, constraints=constraints)" << std::endl;
    out << "psi.setValue(T)" << std::endl;
    out << std::endl;
    
}


void RevScript::printFBDP( std::ostream &out )
{
    
    out << "##############################" << std::endl;
    out << "# fossil-birth-death process #" << std::endl;
    out << "##############################" << std::endl;
    out << std::endl;
    out << "speciation ~ dnExponential(1.0)" << std::endl;
    out << "relative_extinction ~ dnBeta(2.0, 2.0)" << std::endl;
    out << "extinction := speciation * relative_extinction" << std::endl;
    out << "fossil_rate ~ dnExponential(1.0)" << std::endl;
    out << "sampling_prob <- 1.0" << std::endl;
    out << "process_time ~ dnUnif(0.0, 1000.0)" << std::endl;
    out << std::endl;
    out << "psi ~ dnFBDP(origin=process_time, lambda=speciation, mu=extinction, psi=fossil_rate, rho=sampling_prob, condition=\"survival\", names=taxon_names)" << std::endl;
    out << "psi.setValue(T)" << std::endl;
    out << std::endl;
    out << "root_age := psi.rootAge()" << std::endl;
    out << std::endl;
    out << "moves[++mi] = mvScale(speciation,lambda=1.0,tune=true,weight=3.0)" << std::endl;
    out << "moves[++mi] = mvSlide(relative_extinction,delta=1.0,tune=true,weight=3.0)" << std::endl;
    out << "moves[++mi] = mvSlide(process_time,delta=1.0,tune=true,weight=3.0)" << std::endl;
    out << "moves[++mi] = mvSubtreeScale(psi, weight=3.0)" << std::endl;
    out << "moves[++mi] = mvRootTimeSlideUniform(psi, origin=process_time, weight=3.0)" << std::endl;
    out << "moves[++mi] = mvNodeTimeSlideUniform(psi, weight=15.0)" << std::endl;
    out << "moves[++mi] = mvCollapseExpandFossilBranch(psi,process_time,weight=10.0)" << std::endl;
    out << std::endl;
    
}


void RevScript::printInitialSetting( std::ostream &out )
{
    
    out << std::endl;
    out << "####################" << std::endl;
    out << "# initial settings #" << std::endl;
    out << "####################" << std::endl;
    out << std::endl;
    out << "mi = 0" << std::endl;
    out << std::endl;
    
}



void RevScript::printModel( std::ostream &out )
{
    
    out << std::endl;
    out << "################" << std::endl;
    out << "# create model #" << std::endl;
    out << "################" << std::endl;
    out << std::endl;
    out << "my_model = model(psi)" << std::endl;
    out << std::endl;

}


void RevScript::printMonitors( std::ostream &out )
{
    std::stringstream mntr_strm;
    mntr_strm << base_file_name << ".";
    
    if ( prior == BDP )
    {
        mntr_strm << "BDP.";
    }
    else
    {
        mntr_strm << "FBDP.";
    }
    
    if ( mol_only == true )
    {
        mntr_strm << "mol_only.";
    }
    mntr_strm << rep;
    std::string mntr_file =  mntr_strm.str();

    out << std::endl;
    out << "###################" << std::endl;
    out << "# create monitors #" << std::endl;
    out << "###################" << std::endl;
    out << std::endl;
    out << "monitors[1] = mnModel(filename=\"output/" << mntr_file << ".log\", printgen=10)" << std::endl;
    out << "monitors[2] = mnFile(filename=\"output/" << mntr_file << ".trees\", printgen=10, psi)" << std::endl;
    out << "monitors[3] = mnScreen(printgen=1000, root_age)" << std::endl;
    out << std::endl;

}

void RevScript::printMcmc( std::ostream &out )
{
    
    out << std::endl;
    out << "###############" << std::endl;
    out << "# create MCMC #" << std::endl;
    out << "###############" << std::endl;
    out << std::endl;
    out << "my_mcmc = mcmc(my_model, monitors, moves)" << std::endl;
    out << std::endl;
    out << "### pre-burnin to tune the proposals ###" << std::endl;
    out << "my_mcmc.burnin(generations=10000,tuningInterval=250)" << std::endl;
    out << "" << std::endl;
    out << "### run the MCMC ###" << std::endl;
    out << "my_mcmc.run(generations=50000)" << std::endl;
    out << std::endl;

}


void RevScript::printMolecularSubstitutionModel( std::ostream &out )
{
    
    out << std::endl;
    out << "################################" << std::endl;
    out << "# molecular substitution model #" << std::endl;
    out << "################################" << std::endl;
    out << std::endl;
    out << "# substition model priors and parameter" << std::endl;
    out << "pi_prior <- rep(1,4)" << std::endl;
    out << "er_prior <- rep(1,6)" << std::endl;
    out << std::endl;
    out << "pi ~ dnDirichlet(pi_prior)" << std::endl;
    out << "er ~ dnDirichlet(er_prior)" << std::endl;
    out << std::endl;
    out << "# moves on the substitution process parameters" << std::endl;
    out << "moves[++mi] = mvSimplexElementScale(pi, alpha=10.0, tune=true, weight=2.0)" << std::endl;
    out << "moves[++mi] = mvSimplexElementScale(er, alpha=10.0, tune=true, weight=3.0)" << std::endl;
    out << std::endl;
    out << "# the rate matrix" << std::endl;
    out << "Q_mol := fnGTR(er,pi)" << std::endl;
    out << std::endl;
    out << "# the sequence evolution model" << std::endl;
    out << "seq_mol ~ dnPhyloCTMC(tree=psi, Q=Q_mol, branchRates=mol_clock_rate, type=\"DNA\")" << std::endl;
    out << std::endl;
    out << "# attach the data" << std::endl;
    out << "seq_mol.clamp(D_mol)" << std::endl;
    out << std::endl;

}


void RevScript::printMorphologicalSubstitutionModel( std::ostream &out )
{
    
    out << std::endl;
    out << "####################################" << std::endl;
    out << "# morphological substitution model #" << std::endl;
    out << "####################################" << std::endl;
    out << std::endl;
    out << "# the rate matrix" << std::endl;
    out << "Q_morph <- fnJC(2)" << std::endl;
    out << std::endl;
    out << "# the sequence evolution model" << std::endl;
    out << "seq_morph ~ dnPhyloCTMC(tree=psi, Q=Q_morph, branchRates=morph_clock_rate, type=\"Standard\")" << std::endl;
    out << std::endl;
    out << "# attach the data" << std::endl;
    out << "seq_morph.clamp(D_morph)" << std::endl;
    out << std::endl;
    
}


void RevScript::printQuit( std::ostream &out )
{
    
    out << std::endl;
    out << "########" << std::endl;
    out << "# quit #" << std::endl;
    out << "########" << std::endl;
    out << std::endl;
    out << "q()" << std::endl;
    out << std::endl;
    
}



void RevScript::printReadMolecularData( std::ostream &out )
{
    // get file names
    std::stringstream ss_mol;
    ss_mol << base_file_name << ".mol.";
    if ( prior == BDP ) ss_mol << "extant.only.";
    ss_mol << rep << ".nex";
    std::string molfile = ss_mol.str();
    
    out << std::endl;
    out << "#######################" << std::endl;
    out << "# read molecular data #" << std::endl;
    out << "#######################" << std::endl;
    out << std::endl;
    out << "D_mol <- readDiscreteCharacterData(file=\"" << molfile <<   "\")" << std::endl;
    out << std::endl;

}


void RevScript::printReadMorphologicalData( std::ostream &out )
{
    // get file names
    std::stringstream ss_morph;
    ss_morph << base_file_name << ".morph." << rep << ".nex";
    std::string morph_file = ss_morph.str();
    
    out << std::endl;
    out << "###########################" << std::endl;
    out << "# read morphological data #" << std::endl;
    out << "###########################" << std::endl;
    out << std::endl;
    out << "D_morph <- readDiscreteCharacterData(file=\"" << morph_file <<   "\")" << std::endl;
    out << std::endl;
    
}



void RevScript::printReadTree( std::ostream &out )
{
    
    // get file names
    std::stringstream ss_tree;
    ss_tree << base_file_name << ".";
    if ( prior == BDP ) ss_tree << "extant.only.";
    ss_tree << rep << ".tree";
    std::string tree_file = ss_tree.str();
    
    out << std::endl;
    out << "#################" << std::endl;
    out << "# read the tree #" << std::endl;
    out << "#################" << std::endl;
    out << std::endl;
    out << "T <- readTrees(\"" << tree_file << "\")[1]" << std::endl;
    out << "n_taxa <- T.ntips()" << std::endl;
    out << "taxon_names <- T.names()" << std::endl << std::endl;
    out << std::endl;
    
}



void RevScript::printRelativeMorphologicalClock( std::ostream &out )
{
    
    out << std::endl;
    out << "################################" << std::endl;
    out << "# relative morphological clock #" << std::endl;
    out << "################################" << std::endl;
    out << std::endl;
    out << "log_rel_morph_clock_rate ~ dnUniform(-10,2)" << std::endl;
    out << "rel_morph_clock_rate := 10^log_rel_morph_clock_rate" << std::endl;
    out << "morph_clock_rate := rel_morph_clock_rate * mol_clock_rate" << std::endl;
    out << std::endl;
    out << "moves[++mi] = mvSlide(log_rel_morph_clock_rate, weight=1.0)" << std::endl;
    out << std::endl;
    
}



void RevScript::printStrictMolecularClock( std::ostream &out )
{
    
    out << std::endl;
    out << "#######################" << std::endl;
    out << "# strict global clock #" << std::endl;
    out << "#######################" << std::endl;
    out << std::endl;
    out << "log_mol_clock_rate ~ dnUniform(-10,2)" << std::endl;
    out << "mol_clock_rate := 10^log_mol_clock_rate" << std::endl;
    out << std::endl;
    out << "moves[++mi] = mvSlide(log_mol_clock_rate, weight=1.0)" << std::endl;
    out << std::endl;
    
}


    
    
void RevScript::printSummarizeTrees( std::ostream &out )
{
    
    std::stringstream mntr_strm;
    mntr_strm << base_file_name << ".";
    
    if ( prior == BDP )
    {
        mntr_strm << "BDP.";
    }
    else
    {
        mntr_strm << "FBDP.";
    }
    
    if ( mol_only == true )
    {
        mntr_strm << "mol_only.";
    }
    mntr_strm << rep;
    std::string mntr_file =  mntr_strm.str();
        
    out << std::endl;
    out << "###################" << std::endl;
    out << "# summarize trees #" << std::endl;
    out << "###################" << std::endl;
    out << std::endl;
    out << "" << std::endl;
    out << "tt = readTreeTrace(\"output/" << mntr_file << ".trees\", \"clock\")" << std::endl;
    out << "map_tree = mapTree(tt, \"output/" << mntr_file << ".MAP.tre\")" << std::endl;
    out << std::endl;
    
}



void RevScript::printUclnMolecularClock( std::ostream &out )
{
    
    out << std::endl;
    out << "######################" << std::endl;
    out << "# UCLN relaxed clock #" << std::endl;
    out << "######################" << std::endl;
    out << std::endl;
    out << "log_mol_clock_rate ~ dnUniform(-10,2)" << std::endl;
    out << "mol_clock_rate := 10^log_mol_clock_rate" << std::endl;
    out << std::endl;
    out << "moves[++mi] = mvSlide(log_mol_clock_rate, weight=1.0)" << std::endl;
    out << std::endl;
    
}