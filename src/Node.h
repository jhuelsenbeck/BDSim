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
        bool                isExtinct(void) const { return extinct; }
        bool                isFossil(void) const { return fossil; }
        void                removeDescendant(Node* p);
        void                removeDescendants(void);
        void                setAncestor(Node* p) { ancestor = p; }
        void                setIndex(int x) { index = x; }
        void                setIsExtinct(bool tf) { extinct = tf; }
        void                setIsFossil(bool tf) { fossil = tf; }
        void                setName(std::string s) { name = s; }
        void                setTime(double x) { time = x; }

    private:
        int                 index;
        Node*               ancestor;
        std::vector<Node*>  descendants;
        std::string         name;
        double              time;
        bool                extinct;
        bool                fossil;
};


#endif