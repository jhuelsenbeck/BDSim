#include <iostream>
#include <cmath>
#include <string>
#include "Settings.h"



Settings::Settings(int argc, char *argv[]) {

#	if 1
	/* set up fake command-line argument string */
	char *cmdString[15];
	cmdString[ 0] = (char*)"bdsm";
	cmdString[ 1] = (char*)"-lambda";
	cmdString[ 2] = (char*)"1.0";
	cmdString[ 3] = (char*)"-mu";
	cmdString[ 4] = (char*)"0.1";
	cmdString[ 5] = (char*)"-out";
	cmdString[ 6] = (char*)"/Users/johnh/Desktop/test";
	cmdString[ 7] = (char*)"-nt";
	cmdString[ 8] = (char*)"20";
	argc = 9;
	argv = cmdString;
#	endif

	enum Mode { OUTPUT_FILE, LAMBDA, MU, NUM_TAXA, NONE };

	/* set default values for parameters */
	outputFileName = "";
	speciationRate = 1000000;
	extinctionRate = 0;
	numLivingTaxa  = 10;
	
	if (argc > 1)
		{
		if (argc % 2 == 0)
			{
			printUsage();
			}
			
		/* read the command-line arguments */
		int status = NONE;
		for (int i=1; i<argc; i++)
			{
			std::string cmd = argv[i];
			//std::cout << cmd << std::endl;
			if (status == NONE)
				{
				/* read the parameter specifier */
				if ( cmd == "-lambda" )
					status = LAMBDA;
				else if ( cmd == "-mu" )
					status = MU;
				else if ( cmd == "-out" )
					status = OUTPUT_FILE;
				else if ( cmd == "-nt" )
					status = NUM_TAXA;
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
				else if ( status == OUTPUT_FILE )
					outputFileName = argv[i];
				else if ( status == NUM_TAXA )
					numLivingTaxa = atoi(argv[i]);
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

}

void Settings::printUsage(void) {

	std::cout << "Usage:" << std::endl;
	std::cout << "   -out    : Output file name" << std::endl;
	std::cout << "   -nt     : Number of living taxa" << std::endl;
	std::cout << "   -lambda : Speciation rate" << std::endl;
	std::cout << "   -mu     : Extinction rate" << std::endl;
	std::cout << std::endl;
	exit(0);

}