#ifndef Node_H
#define Node_H

#include <vector>


class Node {

	public:
                            Node(void);
        int                 getIndex(void) { return index; }
        void                setIndex(int x) { index = x; }

    private:
        int                 index;
        Node*               ancestor;
        std::vector<Node*>  descendants;
        double              time;
};


#endif