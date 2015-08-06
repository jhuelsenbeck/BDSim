#include <iostream>
#include <set>
#include "MbRandom.h"
#include "Node.h"
#include "Tree.h"



Tree::Tree(MbRandom* rp, double l, double m, int n, double t) {

    ranPtr        = rp;
    lambda        = l;
    mu            = m;
    numLivingTaxa = n;
    duration      = t;
    root          = NULL;

    buildBirthDeathTree();
}

Tree::~Tree(void) {

    clearTree();
}

Node* Tree::addNode(void) {

    Node* newNode = new Node;
    newNode->setIndex((int)nodes.size());
    nodes.push_back(newNode);
    return newNode;
}

void Tree::buildBirthDeathTree(void) {

    bool goodTree = false;
    do
        {
        // this is the list of active nodes
        std::set<Node*> activeNodes;
        
        // we start with one lineage
        root = addNode();
        root->setTime(0.0);
        Node* p = addNode();
        p->setAncestor(root);
        root->addDescendant(p);
        activeNodes.insert(p);
        
        double curT = 0.0;
        while (curT < duration && activeNodes.size() > 0)
            {
            double effectiveRate = (lambda + mu) * activeNodes.size();
            curT += ranPtr->exponentialRv(effectiveRate);
            if (curT < duration)
                {
                // pick a branch at random from the active nodes
                p = randomlyChosenNodeFromSet(activeNodes);
                p->setTime(curT);
                if ( ranPtr->uniformRv() < (mu / (lambda+mu)) )
                    {
                    // extinction event
                    activeNodes.erase(p);
                    }
                else
                    {
                    // speciation event
                    Node* n1 = addNode();
                    Node* n2 = addNode();
                    p->addDescendant(n1);
                    p->addDescendant(n2);
                    n1->setAncestor(p);
                    n2->setAncestor(p);
                    activeNodes.erase(p);
                    activeNodes.insert(n1);
                    activeNodes.insert(n2);
                    }
                }
            }
        for (std::set<Node*>::iterator it = activeNodes.begin(); it != activeNodes.end(); it++)
            (*it)->setTime(duration);
        std::cout << "Number of nodes = " << nodes.size() << " Number of extant taxa = " << activeNodes.size() << " Duration = " << duration << std::endl;
        if (activeNodes.size() != numLivingTaxa)
            {
            clearTree();
            activeNodes.clear();
            }
        else
            {
            goodTree = true;
            }
        } while(goodTree == false);
}

void Tree::clearTree(void) {

    for (int i=0; i<nodes.size(); i++)
        delete nodes[i];
    nodes.clear();
}

Node* Tree::randomlyChosenNodeFromSet(std::set<Node*>& s) {

    int whichElement = (int)(ranPtr->uniformRv() * s.size());
    int k = 0;
    for (std::set<Node*>::iterator it = s.begin(); it != s.end(); it++)
        {
        if (k == whichElement)
            return (*it);
        k++;
        }
    return NULL;
}





