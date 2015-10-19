#include <iostream>
#include "Data.h"
#include "EventHistory.h"
#include "MbRandom.h"
#include "Node.h"
#include "Settings.h"
#include "Tree.h"



Data::Data(int nn, int nc, MbRandom* rp, Settings* sp, Tree* t, double r) {

    numNodes = nn;
    numChar  = nc;
    ranPtr   = rp;
    treePtr  = t;
    rate     = r;
    settingsPtr = sp;
    
    dataMatrix = new unsigned*[numNodes];
    dataMatrix[0] = new unsigned[numNodes * numChar];
    for (int i = 1; i<numNodes; i++)
        dataMatrix[i] = dataMatrix[i-1] + numChar;
    for (int i=0; i<numNodes; i++)
        for (int j=0; j<numChar; j++)
            dataMatrix[i][j] = 0;
    
    simulateMorphologicalCharacters();
//    printExtant(true);
    //print(true);
}

Data::Data(int nn, int nc, MbRandom* rp, Settings* sp, Tree* t, double r, std::vector<double> theta, std::vector<double> pi, double alpha) {

    numNodes = nn;
    numChar  = nc;
    ranPtr   = rp;
    treePtr  = t;
    rate     = r;
    settingsPtr = sp;
    
    dataMatrix = new unsigned*[numNodes];
    dataMatrix[0] = new unsigned[numNodes * numChar];
    for (int i = 1; i<numNodes; i++)
        dataMatrix[i] = dataMatrix[i-1] + numChar;
    for (int i=0; i<numNodes; i++)
        for (int j=0; j<numChar; j++)
            dataMatrix[i][j] = 0;
    
    simulateMolecularCharacters(theta, pi, alpha);
    printExtant(false);
    //print(false);
}


Data::~Data(void) {

    delete [] dataMatrix[0];
    delete [] dataMatrix;
}


char Data::convertToDNA(unsigned int val)
{
    switch ( val )
    {
        case 0:
            return 'A';
        case 1:
            return 'C';
        case 2:
            return 'G';
        case 3:
            return 'T';
            
        default:
            return '?';
    }
    
}

void Data::print(bool isMorph) {

    int maxLen = treePtr->lengthOfLongestName();
    for (int n=0; n<treePtr->getNumberOfDownPassNodes(); n++)
        {
        Node* p = treePtr->getDownPassNode(n);
        if ( treePtr->isExtantTaxon(p) == true || treePtr->isFossilTaxon(p) == true )
            {
            std::string s = p->getName();
            std::cout << s << "   ";
            for (int i=0; i<maxLen-s.length(); i++)
                std::cout << " ";
            for (int c=0; c<numChar; c++)
                {
                if ( isMorph == true )
                    {
                    std::cout << dataMatrix[p->getIndex()][c];
                    }
                else
                    {
                    if (s[0] == 'F')
                        std::cout << "?";
                    else
                        std::cout << convertToDNA(dataMatrix[p->getIndex()][c]);
                    }
                }
            std::cout << std::endl;
            }
        }
}

void Data::printExtant(bool isMorph) {

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
                {
                if ( isMorph == true )
                    {
                    std::cout << dataMatrix[p->getIndex()][c];
                    }
                else
                    {
                    std::cout << convertToDNA(dataMatrix[p->getIndex()][c]);
                    }
                }
            std::cout << std::endl;
            }
        }
    
}

void Data::printNexus(std::iostream &outStream, bool isMorph, bool includeFossils) {
    
    int ntaxa = 0;
    for (int n=0; n<treePtr->getNumberOfDownPassNodes(); n++)
        {
        Node* p = treePtr->getDownPassNode(n);
        if ( treePtr->isExtantTaxon(p) == true || (includeFossils == true && treePtr->isFossilTaxon(p) == true ) )
            {
            ++ntaxa;
            }
        }
    
    outStream << "#NEXUS" << std::endl << std::endl;
    outStream << "begin data;" << std::endl;
    outStream << "   dimensions ntax=" << ntaxa << " nchar=" << numChar << ";" << std::endl;
    if ( isMorph == true )
        {
        outStream << "   format datatype=Standard" << " missing=?"<< ";" << std::endl;
        }
    else
        {
        outStream << "   format datatype=DNA" << " missing=?"<< ";" << std::endl;
        }
    
    
    outStream << "   matrix" << std::endl;
    
    int maxLen = treePtr->lengthOfLongestName();
    for (int n=0; n<treePtr->getNumberOfDownPassNodes(); n++)
        {
        Node* p = treePtr->getDownPassNode(n);
        if ( treePtr->isExtantTaxon(p) == true || (includeFossils == true && treePtr->isFossilTaxon(p) == true ) )
            {
            std::string s = p->getName();
            outStream << "   " << s << "   ";
            for (int i=0; i<maxLen-s.length(); i++)
                outStream << " ";
            for (int c=0; c<numChar; c++)
                {
                if ( isMorph == true )
                    {
                    outStream << dataMatrix[p->getIndex()][c];
                    }
                else
                    {
                    if ( treePtr->isFossilTaxon(p) == true )
                        {
                        outStream << "?";
                        }
                    else
                        {
                        outStream << convertToDNA(dataMatrix[p->getIndex()][c]);
                        }
                    }
                }
            outStream << std::endl;
            }
        }
    outStream << "   ;" << std::endl;
    outStream << "end;" << std::endl;
}


void Data::simulateMolecularCharacters(std::vector<double> theta, std::vector<double> pi, double alpha) {

    // get events of rate change
    EventHistory myEventHistory(treePtr, ranPtr, 50.0, settingsPtr->getMolEventRate());
    //myEventHistory.print();
    
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
    
    // normalize rates
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
    
    // get the basal rate for the characters
    double basalRate = rate;
    
    // simulate
    for (int c=0; c<numChar; c++)
        {
        // get character rate
        double r = 1.0;
        if (alpha < 100.0)
            r = ranPtr->gammaRv(alpha, alpha);
            
        for (int n=treePtr->getNumberOfDownPassNodes()-1; n>=0; n--)
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
                p->setVal(1.0);
                }
            else
                {
                // begin branch sim
                // initialize interval information for this branch
                std::vector<double> intervalStarts;
                std::vector<double> intervalEnds;
                std::vector<double> intervalBasalRate;
                double intStart = p->getAncestor()->getTime();
                double f = p->getAncestor()->getVal();
                for (int i=0; i<myEventHistory.numEventsForBranchWithIndex(p->getIndex()); i++)
                    {
                    Event* e = myEventHistory.getEvent(p->getIndex(), i);
                    intervalEnds.push_back( e->getBranchPos() );
                    intervalStarts.push_back( intStart );
                    intervalBasalRate.push_back( f );
                    f *= e->getRateFactor();
                    intStart = e->getBranchPos();
                    }
                intervalEnds.push_back( p->getTime() );
                intervalStarts.push_back( intStart );
                intervalBasalRate.push_back( f );
                p->setVal( f );
                /*std::cout << "simulating branch " << p->getIndex() << std::endl;
                for (int i=0; i<intervalEnds.size(); i++)
                    std::cout << i << " -- " << intervalStarts[i] << " - " << intervalEnds[i] << " " << intervalBasalRate[i] << std::endl; */
                    
                // simulate along branch
                int curState = dataMatrix[p->getAncestor()->getIndex()][c];
                for (int brI=0; brI<intervalEnds.size(); brI++)
                    {
                    double curT = intervalStarts[brI];
                    double endT = intervalEnds[brI];
                    double rateMultiplier = basalRate * r * intervalBasalRate[brI];

                    while (curT < endT)
                        {
                        double subRate = -q[curState][curState] * rateMultiplier;
                        curT += ranPtr->exponentialRv(subRate);
                        if (curT < endT)
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
                    }
                dataMatrix[p->getIndex()][c] = curState;
                // end branch sim
                }
            }
        }
            
}

void Data::simulateMorphologicalCharacters(void) {

    // get events of rate change
    EventHistory myEventHistory(treePtr, ranPtr, 50.0, settingsPtr->getMolEventRate());
    //myEventHistory.print();

    // get the basal rate for the characters
    double basalRate = rate;

    for (int c=0; c<numChar; c++)
        {
        for (int n=treePtr->getNumberOfDownPassNodes()-1; n>=0; n--)
            {
            Node* p = treePtr->getDownPassNode(n);
            if ( treePtr->isRoot(p) == true )
                {
                // root node
                if (ranPtr->uniformRv() < 0.5)
                    dataMatrix[p->getIndex()][c] = 0;
                else
                    dataMatrix[p->getIndex()][c] = 1;
                p->setVal(1.0);
                }
            else
                {

                // initialize interval information for this branch
                std::vector<double> intervalStarts;
                std::vector<double> intervalEnds;
                std::vector<double> intervalBasalRate;
                double intStart = p->getAncestor()->getTime();
                double f = p->getAncestor()->getVal();
                for (int i=0; i<myEventHistory.numEventsForBranchWithIndex(p->getIndex()); i++)
                    {
                    Event* e = myEventHistory.getEvent(p->getIndex(), i);
                    intervalEnds.push_back( e->getBranchPos() );
                    intervalStarts.push_back( intStart );
                    intervalBasalRate.push_back( f );
                    f *= e->getRateFactor();
                    intStart = e->getBranchPos();
                    }
                intervalEnds.push_back( p->getTime() );
                intervalStarts.push_back( intStart );
                intervalBasalRate.push_back( f );
                p->setVal( f );
                /*std::cout << "simulating branch " << p->getIndex() << std::endl;
                for (int i=0; i<intervalEnds.size(); i++)
                    std::cout << i << " -- " << intervalStarts[i] << " - " << intervalEnds[i] << " " << intervalBasalRate[i] << std::endl; */
                    
                // simulate along branch
                int curState = dataMatrix[p->getAncestor()->getIndex()][c];
                for (int brI=0; brI<intervalEnds.size(); brI++)
                    {
                    double curT = intervalStarts[brI];
                    double endT = intervalEnds[brI];
                    double rateMultiplier = basalRate * 1.0 * intervalBasalRate[brI];
                    while (curT < endT)
                        {
                        double subRate = 1.0 * rateMultiplier;
                        curT += ranPtr->exponentialRv(subRate);
                        if (curT < endT)
                            {
                            if (curState == 0)
                                curState = 1;
                            else
                                curState = 0;
                            }
                        }
                    }
                dataMatrix[p->getIndex()][c] = curState;
                }
            }
        }
            

}



