#ifndef Data_H
#define Data_H

#include <vector>

class MbRandom;
class Settings;
class Tree;

class Data {

	public:
                            Data(int nn, int nc, MbRandom* rp, Settings* sp, Tree* t, double r);
                            Data(int nn, int nc, MbRandom* rp, Settings* sp, Tree* t, double r, std::vector<double> theta, std::vector<double> pi, double alpha);
                           ~Data(void);
        void                print(bool isMorph);
        void                printExtant(bool isMorph);
        void                printNexus(std::iostream &str, bool isMorph, bool includeFossils);

    private:
                            Data(void) { }

        char                convertToDNA(unsigned val);
        void                simulateMorphologicalCharacters(void);
        void                simulateMolecularCharacters(std::vector<double> theta, std::vector<double> pi, double alpha);
        double              rate;
        MbRandom*           ranPtr;
        Tree*               treePtr;
        Settings*           settingsPtr;
        int                 numNodes;
        int                 numChar;
        unsigned**          dataMatrix;
};


#endif