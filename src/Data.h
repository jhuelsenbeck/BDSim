#ifndef Data_H
#define Data_H

#include <vector>

class MbRandom;
class Tree;

class Data {

	public:
                            Data(int nn, int nc, MbRandom* rp, Tree* t, double r);
                            Data(int nn, int nc, MbRandom* rp, Tree* t, double r, std::vector<double> theta, std::vector<double> pi);
                           ~Data(void);
        void                print(void);
        void                printExtant(void);

    private:
                            Data(void) { }
        void                simulateMorphologicalCharacters(void);
        void                simulateMolecularCharacters(std::vector<double> theta, std::vector<double> pi);
        double              rate;
        MbRandom*           ranPtr;
        Tree*               treePtr;
        int                 numNodes;
        int                 numChar;
        unsigned**          dataMatrix;
};


#endif