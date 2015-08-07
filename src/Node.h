#ifndef Node_H
#define Node_H

#include <string>
#include <vector>


class Node {

	public:
                            Node(void);
        void                addDescendant(Node* p);
        Node*               getAncestor(void) { return ancestor; }
        std::vector<Node*>& getDescendants(void) { return descendants; }
        Node*               getDescendant(int idx);
        int                 getIndex(void) { return index; }
        std::string         getName(void) { return name; }
        int                 getNumDescendants(void);
        double              getTime(void) { return time; }
        void                removeDescendants(void);
        void                setAncestor(Node* p) { ancestor = p; }
        void                setIndex(int x) { index = x; }
        void                setName(std::string s) { name = s; }
        void                setTime(double x) { time = x; }

    private:
        int                 index;
        Node*               ancestor;
        std::vector<Node*>  descendants;
        std::string         name;
        double              time;
};


#endif