#ifndef Settings_H
#define Settings_H

#include <string>
#include <vector>


class Settings {

	public:
                            Settings(int argc, char *argv[]);
        double              getSpeciationRate(void) { return speciationRate; }
        std::vector<double> getExchangeabilityParameters(void) { return exchangeabilityParameters; }
        double              getExtinctionRate(void) { return extinctionRate; }
        double              getFossilizationRate(void) { return fossilizationRate; }
        double              getGammaShapeParameter(void) { return gammaShapeParameter; }
        int                 getNumLivingTaxa(void) { return numLivingTaxa; }
        int                 getNumMolecularCharacters(void) { return numMolecularCharacters; }
        int                 getNumMorphologicalCharacters(void) { return numMorphologicalCharacters; }
        int                 getNumReplicates(void) { return numReplicates; }
        std::string         getOutputFileName(void) { return outputFileName; }
        std::vector<double> getStationaryFrequenciesParameters(void) { return stationaryFrequenciesParameters; }
        void                print(void);
        void                setSpeciationRate(double x) { speciationRate = x; }
        void                setExchangeabilityParameters(std::vector<double> x) { exchangeabilityParameters = x; }
        void                setExtinctionRate(double x) { extinctionRate = x; }
        void                setFossilizationRate(double x) { fossilizationRate = x; }
        void                setGammaShapeParameter(double x) { gammaShapeParameter = x; }
        void                setNumLivingTaxa(int x) { numLivingTaxa = x; }
        void                setNumMolecularCharacters(int x) { numMolecularCharacters = x; }
        void                setNumMorphologicalCharacters(int x) { numMorphologicalCharacters = x; }
        void                setOutputFileName(std::string s) { outputFileName = s; }
        void                setStationaryFrequenciesParameters(std::vector<double> x) { stationaryFrequenciesParameters = x; }

    private:
        void                preprocessStr(int* argc, char *argv[]);
        void                printUsage(void);
        std::string         outputFileName;
        double              speciationRate;
        double              extinctionRate;
        double              fossilizationRate;
        int                 numLivingTaxa;
        int                 numMorphologicalCharacters;
        int                 numMolecularCharacters;
        double              gammaShapeParameter;
        std::vector<double> exchangeabilityParameters;
        std::vector<double> stationaryFrequenciesParameters;
        std::vector< std::vector<double> >  tempVectors;
        int                 numReplicates;
};


#endif