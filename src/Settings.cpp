#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include "Settings.h"



Settings::Settings(int argc, char *argv[]) {

#	if 1
	/* set up fake command-line argument string */
	char *cmdString[20];
	cmdString[ 0] = (char*)"bdsm";
	cmdString[ 1] = (char*)"-lambda";
	cmdString[ 2] = (char*)"0.5";
	cmdString[ 3] = (char*)"-mu";
	cmdString[ 4] = (char*)"0.1";
	cmdString[ 5] = (char*)"-phi";
	cmdString[ 6] = (char*)"0.1";
	cmdString[ 7] = (char*)"-outFile";
    cmdString[ 8] = (char*)"test";
    cmdString[ 9] = (char*)"-outPath";
    cmdString[10] = (char*)"/Users/hoehna/Desktop/LondonSimulations/data";
    cmdString[11] = (char*)"-nr";
    cmdString[12] = (char*)"1";
	cmdString[13] = (char*)"-nt";
	cmdString[14] = (char*)"20";
    cmdString[15] = (char*)"-bf";
	cmdString[16] = (char*)"(0.25, 0.25, 0.25, 0.25)";
    cmdString[17] = (char*)"-exch";
	cmdString[18] = (char*)"(1.0, 1.0, 1.0, 1.0, 1.0, 1.0)";
	argc = 19;
	argv = cmdString;
#	endif

    // process the string looking for number lists
    preprocessStr(&argc, argv);
    
	enum Mode { OUTPUT_FILE, OUTPUT_PATH, NUM_REPS, LAMBDA, MU, PHI, NUM_TAXA, NUM_MORPH, NUM_MOL, GAMMA_SHAPE, EXC_PARM, BASE_FREQ, NONE };

	/* set default values for parameters */
    numReplicates                   = 1;
    outputFileName                  = "";
    outputFilePath                  = "";
	speciationRate                  = 1.0;
	extinctionRate                  = 0.5;
    fossilizationRate               = 0.5;
	numLivingTaxa                   = 10;
    numMorphologicalCharacters      = 50;
    numMolecularCharacters          = 50;
    gammaShapeParameter             = 1000.0;
    exchangeabilityParameters       = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
    stationaryFrequenciesParameters = { 0.25, 0.25, 0.25, 0.25 };
	
    int vecNum = 0;
	if (argc > 1)
		{
		if (argc % 2 == 0)
			{
			printUsage();
			}
			
		/* read the command-line arguments */
		int status = NONE;
		for (int i=1; i<cmds.size(); i++)
			{
			std::string cmd = cmds[i];
			std::cout << cmd << std::endl;
			if (status == NONE)
				{
				/* read the parameter specifier */
				if ( cmd == "-lambda" )
					status = LAMBDA;
				else if ( cmd == "-mu" )
					status = MU;
				else if ( cmd == "-phi" )
					status = PHI;
				else if ( cmd == "-outFile" )
                    status = OUTPUT_FILE;
                else if ( cmd == "-outPath" )
                    status = OUTPUT_PATH;
				else if ( cmd == "-nt" )
					status = NUM_TAXA;
                else if ( cmd == "-nr" )
                    status = NUM_REPS;
				else if ( cmd == "-nmorph" )
					status = NUM_MORPH;
				else if ( cmd == "-nmol" )
					status = NUM_MOL;
				else if ( cmd == "-gamma" )
					status = GAMMA_SHAPE;
				else if ( cmd == "-exch" )
					status = EXC_PARM;
				else if ( cmd == "-bf" )
					status = BASE_FREQ;
				else
					{
					std::cerr << "Could not interpret option \"" << cmd << "\"." << std::endl;
					exit(1);
					}
				}
			else
				{
				/* read the parameter */
				if ( status == LAMBDA )
					speciationRate = atof(argv[i]);
				else if ( status == MU )
					extinctionRate = atof(argv[i]);
				else if ( status == PHI )
					fossilizationRate = atof(argv[i]);
				else if ( status == OUTPUT_FILE )
                    outputFileName = argv[i];
                else if ( status == OUTPUT_PATH )
                    outputFilePath = argv[i];
				else if ( status == NUM_TAXA )
					numLivingTaxa = atoi(argv[i]);
                else if ( status == NUM_REPS )
                    numReplicates = atoi(argv[i]);
				else if ( status == NUM_MORPH )
					numMorphologicalCharacters = atoi(argv[i]);
				else if ( status == NUM_MOL )
					numMolecularCharacters = atoi(argv[i]);
				else if ( status == GAMMA_SHAPE )
					gammaShapeParameter = atof(argv[i]);
				else if ( status == EXC_PARM )
                    exchangeabilityParameters = tempVectors[vecNum++];
				else if ( status == BASE_FREQ )
					stationaryFrequenciesParameters = tempVectors[vecNum++];
				else
					{
					std::cerr << "Unknown status reading command line information" << std::endl;
					exit(1);
					}
				status = NONE;
				}
			}
		}
	else
		{
		printUsage();
		}	

    print();
}

void Settings::preprocessStr(int* argc, char *argv[]) {

//    std::vector<std::string> newCmdStr;
    
    bool readingVector = false;
    std::vector<double> v;
    for (int i=0; i<(*argc); i++)
        {
        std::string s = argv[i];
        std::string tempStr = "";
        bool isLineVector = false, isLineEndVector = false;
        for (int j=0; j<s.length(); j++)
            {
            if (s[j] == '(')
                {
                readingVector = true;
                isLineVector = true;
                tempStr = "";
                v.clear();
                }
            else if (s[j] == ')')
                {
				double x;
				std::istringstream buf(tempStr);
				buf >> x;
                v.push_back(x);
                tempVectors.push_back(v);
                v.clear();
                tempStr = "";
                readingVector = false;
                isLineVector = true;
                isLineEndVector = true;
                }
            else if (s[j] == ',' && readingVector == true)
                {
				double x;
				std::istringstream buf(tempStr);
				buf >> x;
                tempStr = "";
                v.push_back(x);
                isLineVector = true;
                }
            else if ( (isdigit(s[j]) == true || s[j] == '.') && readingVector == true )
                {
                isLineVector = true;
                tempStr += s[j];
                }
            }
            
        if (isLineVector == false)
            cmds.push_back(s);
        else
            {
            if (isLineEndVector == true)
                cmds.push_back("xxx");
            }
        
        }
//    
//    for (int i=0; i<newCmdStr.size(); i++)
//        argv[i] = (char*)newCmdStr[i].c_str();
}

void Settings::print(void) {

    std::cout << "Output file name                = \"" << outputFileName << "\"" << std::endl;
    std::cout << "Output file path                = \"" << outputFilePath << "\"" << std::endl;
    std::cout << "Number of simulation replicates = " << numReplicates << std::endl;
    std::cout << "Speciation rate                 = " << speciationRate << std::endl;
    std::cout << "Extinction rate                 = " << extinctionRate << std::endl;
    std::cout << "Fossilization rate              = " << fossilizationRate << std::endl;
    std::cout << "No. of extant taxa              = " << numLivingTaxa << std::endl;
    std::cout << "No. of morphological characters = " << numMorphologicalCharacters << std::endl;
    std::cout << "No. of molecular characters     = " << numMolecularCharacters << std::endl;
    std::cout << "Gamma shape parameter           = " << gammaShapeParameter << std::endl;
    std::cout << "Exchangeability parameters      = ( ";
    for (int i=0; i<exchangeabilityParameters.size(); i++)
        std::cout << exchangeabilityParameters[i] << " ";
    std::cout << ")" << std::endl;
    std::cout << "Base frequencies                = ( ";
    for (int i=0; i<stationaryFrequenciesParameters.size(); i++)
        std::cout << stationaryFrequenciesParameters[i] << " ";
    std::cout << ")" << std::endl;
}

void Settings::printUsage(void) {

	std::cout << "Usage:" << std::endl;
    std::cout << "   -outFile : Output file name" << std::endl;
    std::cout << "   -outPath : Output file path" << std::endl;
    std::cout << "   -nr      : Number of simulation replicates" << std::endl;
	std::cout << "   -nt      : Number of living taxa" << std::endl;
	std::cout << "   -lambda  : Speciation rate" << std::endl;
	std::cout << "   -mu      : Extinction rate" << std::endl;
	std::cout << "   -phi     : Fossilization rate" << std::endl;
	std::cout << "   -nmorph  : Number of morphological characters" << std::endl;
	std::cout << "   -nmol    : Number of molecular characters" << std::endl;
	std::cout << "   -gamma   : Gamma shape parameter" << std::endl;
	std::cout << "   -exch    : Exchangeability parameters" << std::endl;
	std::cout << "   -bf      : Base frequencies" << std::endl;
	std::cout << std::endl;
	exit(0);

}