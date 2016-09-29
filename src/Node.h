#ifndef Node_H
#define Node_H

#include <string>
#include <vector>


class Node {

	public:
                                    Node(void);
        void                        addDescendant(Node* p);
        void                        addTaxonToBipartition(std::string s) { taxonBipartition.push_back(s); }
        Node*                       getAncestor(void) { return ancestor; }
        std::vector<Node*>&         getDescendants(void) { return descendants; }
        Node*                       getDescendant(int idx);
        int                         getIndex(void) { return index; }
        std::string                 getName(void) { return name; }
        int                         getNumDescendants(void);
        double                      getReconstructedTime();
        double                      getTime(void) { return time; }
        std::vector<std::string>&   getTaxonBipartition(void) { return taxonBipartition; }
        bool                        isExtinct(void) const { return extinct; }
        bool                        isFossil(void) const { return fossil; }
        bool                        getFlag(void) { return flag; }
        bool                        getIsInReconstructedTree(void) { return isInReconstructedTree; }
        double                      getVal(void) { return val; }
        bool                        getIsPunctuatedBranch(void) { return isPunctuatedBranch; }
        void                        removeDescendant(Node* p);
        void                        removeDescendants(void);
        void                        setAncestor(Node* p) { ancestor = p; }
        void                        setIndex(int x) { index = x; }
        void                        setIsExtinct(bool tf) { extinct = tf; }
        void                        setIsFossil(bool tf) { fossil = tf; }
        void                        setName(std::string s) { name = s; }
        void                        setTime(double x) { time = x; }
        void                        setFlag(bool tf) { flag = tf; }
        void                        setIsInReconstructedTree(bool tf) { isInReconstructedTree = tf; }
        void                        setVal(double x) { val = x; }
        void                        setIsPunctuatedBranch(bool tf) { isPunctuatedBranch = tf; }
    
    private:
        int                         index;
        Node*                       ancestor;
        std::vector<Node*>          descendants;
        std::string                 name;
        double                      time;
        bool                        extinct;
        bool                        fossil;
        bool                        flag;
        bool                        isInReconstructedTree;
        std::vector<std::string>    taxonBipartition;
        double                      val;
        bool                        isPunctuatedBranch;
};


#endif