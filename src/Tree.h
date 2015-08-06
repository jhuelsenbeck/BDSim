#ifndef Tree_H
#define Tree_H

#include <set>
#include <sstream>
#include <string>
#include <vector>
class MbRandom;
class Node;


class Tree {
	public:
                            Tree(MbRandom* rp, double l, double m, double p, int n, double t);
                            Tree(Tree& t);
                           ~Tree(void);
        std::string         getNewick(void);
        void                printTree(void);
    
    private:
        Node*               addNode(void);
        void                copyTree(Tree& t);
        int                 dex(Node* p);
        void                clearTree(void);
        void                buildBirthDeathTree(void);
        void                pruneToReconstructedProcess(void);
        Node*               randomlyChosenNodeFromSet(std::set<Node*>& s);
        void                showTree(Node* p, int indent);
        void                writeTree(Node* p, std::stringstream& ss);

        MbRandom*           ranPtr;
        double              lambda;
        double              mu;
        double              phi;
        double              duration;
        int                 numLivingTaxa;
        std::vector<Node*>  nodes;
        Node*               root;
};


#endif