#ifndef Tree_H
#define Tree_H

#include <set>
#include <vector>
class MbRandom;
class Node;


class Tree {
	public:
                            Tree(MbRandom* rp, double l, double m, int n, double t);
                           ~Tree(void);
    private:
        Node*               addNode(void);
        void                clearTree(void);
        void                buildBirthDeathTree(void);
        Node*               randomlyChosenNodeFromSet(std::set<Node*>& s);
        MbRandom*           ranPtr;
        double              lambda;
        double              mu;
        double              duration;
        int                 numLivingTaxa;
        std::vector<Node*>  nodes;
        Node*               root;
};


#endif