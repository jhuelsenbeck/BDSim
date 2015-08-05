#ifndef Tree_H
#define Tree_H

#include <vector>
class Node;


class Tree {

	public:
                            Tree(int argc, char *argv[]);

    private:
        int                 numLivingTaxa;
        std::vector<Node*>  nodes;
};


#endif