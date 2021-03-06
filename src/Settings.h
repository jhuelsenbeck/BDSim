#ifndef Settings_H
#define Settings_H

#include <string>
#include <vector>


class Settings {

	public:
                                            Settings(int argc, char *argv[]);
        double                              getSpeciationRate(void) { return speciationRate; }
        std::vector<double>                 getExchangeabilityParameters(void) { return exchangeabilityParameters; }
        double                              getExtinctionRate(void) { return extinctionRate; }
        double                              getFossilizationRate(void) { return fossilizationRate; }
        double                              getGammaShapeParameter(void) { return gammaShapeParameter; }
        int                                 getNumLivingTaxa(void) { return numLivingTaxa; }
        int                                 getNumMolecularCharacters(void) { return numMolecularCharacters; }
        int                                 getNumMorphologicalCharacters(void) { return numMorphologicalCharacters; }
        int                                 getNumReplicates(void) { return numReplicates; }
        std::string                         getFileName(void) { return outputFileName; }
        std::string                         getFilePath(void) { return outputFilePath; }
        std::string                         getFullFileName(void) { return outputFilePath + "/" + outputFileName; }
        std::vector<double>                 getStationaryFrequenciesParameters(void) { return stationaryFrequenciesParameters; }
        double                              getMorphologicalRate(void) { return morphologicalRate; }
        double                              getMolecularRate(void) { return molecularRate; }
        double                              getExtinctionToSpeciationRate(void) { return extinctionToSpeciationRate; }
        double                              getSimDuration(void) { return simDuration; }
        double                              getMolEventRate(void) { return molEventRate; }
        double                              getMolEventMag(void) { return molEventMag; }
        double                              getMorphEventRate(void) { return morphEventRate; }
        double                              getMorphEventMag(void) { return morphEventMag; }
        double                              getPunctucatedFraction(void) { return punctucatedFraction; }
        void                                print(void);
        void                                setSpeciationRate(double x) { speciationRate = x; }
        void                                setExchangeabilityParameters(std::vector<double> x) { exchangeabilityParameters = x; }
        void                                setExtinctionRate(double x) { extinctionRate = x; }
        void                                setFossilizationRate(double x) { fossilizationRate = x; }
        void                                setGammaShapeParameter(double x) { gammaShapeParameter = x; }
        void                                setNumLivingTaxa(int x) { numLivingTaxa = x; }
        void                                setNumMolecularCharacters(int x) { numMolecularCharacters = x; }
        void                                setNumMorphologicalCharacters(int x) { numMorphologicalCharacters = x; }
        void                                setOutputFileName(std::string s) { outputFileName = s; }
        void                                setOutputFilePath(std::string s) { outputFilePath = s; }
        void                                setStationaryFrequenciesParameters(std::vector<double> x) { stationaryFrequenciesParameters = x; }
        void                                setMorphologicalRate(double x) { morphologicalRate = x; }
        void                                setMolecularRate(double x) { molecularRate = x; }
        void                                setExtinctionToSpeciationRate(double x) { extinctionToSpeciationRate = x; }
        void                                setSimDuration(double x) { simDuration = x; }
        void                                setMolEventRate(double x) { molEventRate = x; }
        void                                setMolEventMag(double x) { molEventMag = x; }
        void                                setMorphEventRate(double x) { morphEventRate = x; }
        void                                setMorphEventMag(double x) { morphEventMag = x; }
        void                                setPunctucatedFraction(double x) { punctucatedFraction = x; }
    
    private:
        void                                preprocessStr(int* argc, char *argv[]);
        void                                printUsage(void);
        std::string                         outputFileName;
        std::string                         outputFilePath;
        double                              speciationRate;
        double                              extinctionRate;
        double                              extinctionToSpeciationRate;
        double                              fossilizationRate;
        double                              morphologicalRate;
        double                              molecularRate;
        int                                 numLivingTaxa;
        int                                 numMorphologicalCharacters;
        int                                 numMolecularCharacters;
        double                              gammaShapeParameter;
        std::vector<double>                 exchangeabilityParameters;
        std::vector<double>                 stationaryFrequenciesParameters;
        std::vector< std::vector<double> >  tempVectors;
        int                                 numReplicates;
        double                              simDuration;
        double                              molEventRate;
        double                              molEventMag;
        double                              morphEventRate;
        double                              morphEventMag;
        double                              punctucatedFraction;
        std::vector<std::string>            cmds;
};


#endif