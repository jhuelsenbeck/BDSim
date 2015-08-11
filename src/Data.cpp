#include <iostream>
#include "Data.h"
#include "MbRandom.h"
#include "Node.h"
#include "Tree.h"



Data::Data(int nn, int nc, MbRandom* rp, Tree* t, double r) {

    numNodes = nn;
    numChar  = nc;
    ranPtr   = rp;
    treePtr  = t;
    rate     = r;
    
    dataMatrix = new unsigned*[numNodes];
    dataMatrix[0] = new unsigned[numNodes * numChar];
    for (int i = 1; i<numNodes; i++)
        dataMatrix[i] = dataMatrix[i-1] + numChar;
    for (int i=0; i<numNodes; i++)
        for (int j=0; j<numChar; j++)
            dataMatrix[i][j] = 0;
    
    simulateMorphologicalCharacters();
    printExtant();
}

Data::Data(int nn, int nc, MbRandom* rp, Tree* t, double r, std::vector<double> theta, std::vector<double> pi, double alpha) {

    numNodes = nn;
    numChar  = nc;
    ranPtr   = rp;
    treePtr  = t;
    rate     = r;
    
    dataMatrix = new unsigned*[numNodes];
    dataMatrix[0] = new unsigned[numNodes * numChar];
    for (int i = 1; i<numNodes; i++)
        dataMatrix[i] = dataMatrix[i-1] + numChar;
    for (int i=0; i<numNodes; i++)
        for (int j=0; j<numChar; j++)
            dataMatrix[i][j] = 0;
    
    simulateMolecularCharacters(theta, pi, alpha);
    printExtant();
}


Data::~Data(void) {

    delete [] dataMatrix[0];
    delete [] dataMatrix;
}

void Data::print(void) {

    int maxLen = treePtr->lengthOfLongestName();
    for (int n=0; n<treePtr->getNumberOfDownPassNodes(); n++)
        {
        Node* p = treePtr->getDownPassNode(n);
        if ( treePtr->isExtantTaxon(p) == true ||  treePtr->isFossilTaxon(p) == true )
            {
            std::string s = p->getName();
            std::cout << s << "   ";
            for (int i=0; i<maxLen-s.length(); i++)
                std::cout << " ";
            for (int c=0; c<numChar; c++)
                std::cout << dataMatrix[p->getIndex()][c];
            std::cout << std::endl;
            }
        }
}

void Data::printExtant(void) {

    int maxLen = treePtr->lengthOfLongestName();
    for (int n=0; n<treePtr->getNumberOfDownPassNodes(); n++)
        {
        Node* p = treePtr->getDownPassNode(n);
        if ( treePtr->isExtantTaxon(p) == true )
            {
            std::string s = p->getName();
            std::cout << s << "   ";
            for (int i=0; i<maxLen-s.length(); i++)
                std::cout << " ";
            for (int c=0; c<numChar; c++)
                std::cout << dataMatrix[p->getIndex()][c];
            std::cout << std::endl;
            }
        }
}

void Data::simulateMolecularCharacters(std::vector<double> theta, std::vector<double> pi, double alpha) {

    // set up rate matrix
    double q[4][4];
    int k = 0;
    for (int i=0; i<4; i++)
        {
        for (int j=i+1; j<4; j++)
            {
            q[i][j] = theta[k] * pi[j];
            q[j][i] = theta[k] * pi[j];
            k++;
            }
        }
    double averageRate = 0.0;
    for (int i=0; i<4; i++)
        {
        double sum = 0.0;
        for (int j=0; j<4; j++)
            {
            if (i != j)
                sum += q[i][j];
            }
        q[i][i] = -sum;
        averageRate += pi[i] * sum;
        }
    for (int i=0; i<4; i++)
        {
        for (int j=0; j<4; j++)
            q[i][j] /= averageRate;
        }
    
    // simulate
    for (int c=0; c<numChar; c++)
        {
        double r = 1.0;
        if (alpha < 100.0)
            r = ranPtr->gammaRv(alpha, alpha);
        for (int n=0; n<treePtr->getNumberOfDownPassNodes(); n++)
            {
            Node* p = treePtr->getDownPassNode(n);
            if ( treePtr->isRoot(p) == true )
                {
                // root node
                double u = ranPtr->uniformRv();
                double sum = 0.0;
                for (int i=0; i<4; i++)
                    {
                    sum += pi[i];
                    if (u < sum)
                        {
                        dataMatrix[p->getIndex()][c] = i;
                        break;
                        }
                    }
                }
            else
                {
                double t = (p->getTime() - p->getAncestor()->getTime()) * r;
                double curT = 0.0;
                int curState = dataMatrix[p->getAncestor()->getIndex()][c];
                while (curT < t)
                    {
                    curT += ranPtr->exponentialRv(-q[curState][curState]);
                    if (curT < t)
                        {
                        double u = ranPtr->uniformRv();
                        double sum = 0.0;
                        for (int i=0; i<4; i++)
                            {
                            if (i != curState)
                                {
                                sum += (-q[curState][i] / q[curState][curState]);
                                if (u < sum)
                                    {
                                    curState = i;
                                    break;
                                    }
                                }
                            }
                        }
                    }
                dataMatrix[p->getIndex()][c] = curState;
                }
            }
        }
            
}

void Data::simulateMorphologicalCharacters(void) {

    for (int c=0; c<numChar; c++)
        {
        bool isVariable = false;
        do
            {
            for (int n=0; n<treePtr->getNumberOfDownPassNodes(); n++)
                {
                Node* p = treePtr->getDownPassNode(n);
                if ( treePtr->isRoot(p) == true )
                    {
                    // root node
                    if (ranPtr->uniformRv() < 0.5)
                        dataMatrix[p->getIndex()][c] = 0;
                    else
                        dataMatrix[p->getIndex()][c] = 1;
                    }
                else
                    {
                    double t = (p->getTime() - p->getAncestor()->getTime());
                    double curT = 0.0;
                    int curState = dataMatrix[p->getAncestor()->getIndex()][c];
                    while (curT < t)
                        {
                        curT += ranPtr->exponentialRv(rate);
                        if (curT < t)
                            {
                            if (curState == 0)
                                curState = 1;
                            else
                                curState = 0;
                            }
                        }
                    dataMatrix[p->getIndex()][c] = curState;
                    }
                }
                
                
            std::vector<Node*> living = treePtr->getExtantTaxa();
            int numZeros = 0;
            for (int i=0; i<living.size(); i++)
                {
                if (dataMatrix[ living[i]->getIndex() ][c] == 0)
                    numZeros++;
                }
            if (numZeros == 0 || numZeros == living.size())
                isVariable = false;
            else
                isVariable = true;
            } while (isVariable == false);
            
        }
            

}



