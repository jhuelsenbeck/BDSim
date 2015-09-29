#ifndef Tree_H
#define Tree_H

#include <set>
#include <sstream>
#include <string>
#include <vector>
class FossilCalibration;
class MbRandom;
class Node;


class Tree {
	public:
                                        Tree(MbRandom* rp, double l, double m, double p, int n, double t, bool reconst=true);
                                        Tree(Tree& t);
                                       ~Tree(void);
        Node*                           getDownPassNode(size_t idx) { return downPassSequence[idx]; }
        std::vector<Node*>&             getExtantTaxa(void) { return extantTaxa; }
        std::string                     getNewick(void);
        int                             getNumberOfDownPassNodes(void) { return (int)downPassSequence.size(); }
        int                             getNumberOfNodes(void) { return (int)nodes.size(); }
        int                             getNumberOfCalibrations(void) { return (int)calibrations.size(); }
        FossilCalibration*              getCalibrationIndexed(int idx) { return calibrations[idx]; }
        bool                            isExtantTaxon(Node* p);
        bool                            isFossilTaxon(Node* p);
        bool                            isRoot(Node* p) { return (p == root); }
        int                             lengthOfLongestName(void);
        void                            listCalibrations(void);
        void                            printTree(void);
        void                            markReconstructedTree(void);
    
    private:
        Node*                           addNode(void);
        void                            copyTree(Tree& t);
        int                             dex(Node* p);
        void                            clearTree(void);
        void                            buildBirthDeathTree(void);
        void                            initializeCalibrations(void);
        void                            initializeDownPassSequence(void);
        Node*                           mrcaOfExtantTaxa(void);
        void                            passDown(Node* p);
        void                            pruneToReconstructedProcess(void);
        Node*                           randomlyChosenNodeFromSet(std::set<Node*>& s);
        void                            showTree(Node* p, int indent);
        void                            writeTree(Node* p, std::stringstream& ss);

        MbRandom*                       ranPtr;
        double                          lambda;
        double                          mu;
        double                          phi;
        double                          duration;
        int                             numLivingTaxa;
        std::vector<Node*>              nodes;
        std::vector<Node*>              extantTaxa;
        std::vector<Node*>              fossilTaxa;
        std::vector<Node*>              downPassSequence;
        Node*                           root;
        std::vector<FossilCalibration*> calibrations;
};


#endif