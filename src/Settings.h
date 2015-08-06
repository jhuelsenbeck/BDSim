#ifndef Settings_H
#define Settings_H

#include <string>


class Settings {

	public:
                        Settings(int argc, char *argv[]);
        double          getSpeciationRate(void) { return speciationRate; }
        double          getExtinctionRate(void) { return extinctionRate; }
        double          getFossilizationRate(void) { return fossilizationRate; }
        int             getNumLivingTaxa(void) { return numLivingTaxa; }
        std::string     getOutputFileName(void) { return outputFileName; }
        void            setSpeciationRate(double x) { speciationRate = x; }
        void            setExtinctionRate(double x) { extinctionRate = x; }
        void            setFossilizationRate(double x) { fossilizationRate = x; }
        void            setNumLivingTaxa(int x) { numLivingTaxa = x; }
        void            setOutputFileName(std::string s) { outputFileName = s; }

    private:
        void            printUsage(void);
        std::string     outputFileName;
        double          speciationRate;
        double          extinctionRate;
        double          fossilizationRate;
        int             numLivingTaxa;
};


#endif