#ifndef Node_H
#define Node_H

#include <vector>


class Node {

	public:
                            Node(void);
        void                addDescendant(Node* p);
        Node*               getAncestor(void) { return ancestor; }
        int                 getIndex(void) { return index; }
        double              getTime(void) { return time; }
        void                setAncestor(Node* p) { ancestor = p; }
        void                setIndex(int x) { index = x; }
        void                setTime(double x) { time = x; }

    private:
        int                 index;
        Node*               ancestor;
        std::vector<Node*>  descendants;
        double              time;
};


#endif