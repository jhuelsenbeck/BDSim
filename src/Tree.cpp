#include <iomanip>
#include <iostream>
#include <set>
#include "FossilCalibration.h"
#include "MbRandom.h"
#include "Node.h"
#include "Tree.h"



Tree::Tree(MbRandom* rp, double l, double m, double p, int n, double t, bool r) {

    ranPtr        = rp;
    lambda        = l;
    mu            = m;
    phi           = p;
    numLivingTaxa = n;
    duration      = t;
    root          = NULL;

    buildBirthDeathTree();
    
    if ( r == true )
        {
        pruneToReconstructedProcess();
        initializeDownPassSequence();
        }
    markReconstructedTree();
    initializeCalibrations();
}

Tree::Tree(Tree& t) {

    copyTree(t);
}

Tree::~Tree(void) {

    clearTree();
}

Node* Tree::addNode(void) {

    Node* newNode = new Node;
    newNode->setIndex((int)nodes.size());
    char cStr[32];
    sprintf(cStr, "Node_%d", newNode->getIndex());
    std::string cppStr = cStr;
    newNode->setName(cppStr);
    nodes.push_back(newNode);
    return newNode;
}

void Tree::buildBirthDeathTree(void) {

    bool goodTree = false;
    do
        {
        // this is the list of active nodes
        std::set<Node*> activeNodes;
        
        // we start with two lineages
        root = addNode();
        root->setTime(0.0);
        Node* p = addNode();
        p->setAncestor(root);
        root->addDescendant(p);
        activeNodes.insert(p);
        p = addNode();
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
                double u = ranPtr->uniformRv();
                if ( u <= (mu / (lambda+mu+phi)) )
                    {
                    // extinction event
                    activeNodes.erase(p);
                    p->setIsExtinct(true);
                    }
                else if ( u > (mu / (lambda+mu+phi)) && u <= ((mu+phi) / (lambda+mu+phi)) )
                    {
                    // fossilization event
                    Node* n1 = addNode();
                    p->addDescendant(n1);
                    p->setIsFossil(true);
                    n1->setAncestor(p);
                    activeNodes.erase(p);
                    activeNodes.insert(n1);
                    fossilTaxa.push_back(p);
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
        initializeDownPassSequence();
        for (std::set<Node*>::iterator it = activeNodes.begin(); it != activeNodes.end(); it++)
            {
            (*it)->setTime(duration);
            extantTaxa.push_back(*it);
            }
        std::cout << "Number of nodes = " << nodes.size() << " Number of extant taxa = " << activeNodes.size() << " Duration = " << duration << std::endl;
        if ( (activeNodes.size() != numLivingTaxa) || (mrcaOfExtantTaxa() != root) || fossilTaxa.size() == 0 )
            {
            clearTree();
            activeNodes.clear();
            }
        else
            {
            goodTree = true;
            }
        } while(goodTree == false);
    
    for (size_t i=0; i<extantTaxa.size(); i++)
        {
        char cStr[20];
        sprintf(cStr, "Taxon_%lu", i+1);
        std::string cppStr = cStr;
        extantTaxa[i]->setName(cppStr);
        }
    for (size_t i=0; i<fossilTaxa.size(); i++)
        {
        char cStr[20];
        sprintf(cStr, "Fossil_%lu", i+1);
        std::string cppStr = cStr;
        fossilTaxa[i]->setName(cppStr);
        }
}

void Tree::clearTree(void) {

    for (int i=0; i<nodes.size(); i++)
        delete nodes[i];
    nodes.clear();
    extantTaxa.clear();
    fossilTaxa.clear();
    downPassSequence.clear();
}

void Tree::copyTree(Tree& t) {

    ranPtr        = t.ranPtr;
    lambda        = t.lambda;
    mu            = t.mu;
    phi           = t.phi;
    duration      = t.duration;
    numLivingTaxa = t.numLivingTaxa;
    
    clearTree();
    for (int i=0; i<t.nodes.size(); i++)
        addNode();
    
    for (int i=0; i<t.nodes.size(); i++)
        {
        Node* copyFrom = t.nodes[i];
        Node* copyTo   = nodes[i];
        
        if (copyFrom == t.root)
            root = copyTo;
        
        copyTo->setTime(copyFrom->getTime());
        if (copyFrom->getAncestor() != NULL)
            copyTo->setAncestor( nodes[copyFrom->getAncestor()->getIndex()] );
        else
            copyTo->setAncestor(NULL);
            
        std::vector<Node*> fromDesc = copyFrom->getDescendants();
        for (std::vector<Node*>::iterator it = fromDesc.begin(); it != fromDesc.end(); it++)
            {
            copyTo->addDescendant( nodes[(*it)->getIndex()] );
            }
        }
    
    extantTaxa.clear();
    for (int i=0; i<t.extantTaxa.size(); i++)
        extantTaxa.push_back( nodes[t.extantTaxa[i]->getIndex()] );

    fossilTaxa.clear();
    for (int i=0; i<t.fossilTaxa.size(); i++)
        fossilTaxa.push_back( nodes[t.fossilTaxa[i]->getIndex()] );
    
    downPassSequence.clear();
    for (int i=0; i<t.downPassSequence.size(); i++)
        downPassSequence.push_back( nodes[t.downPassSequence[i]->getIndex()] );
}

int Tree::dex(Node* p) {

	if (p == NULL)
		return -1;
	else 
		return p->getIndex();
}

std::string Tree::getNewick(void) {

	std::stringstream ss;
	writeTree(root, ss);
	std::string newick = ss.str();
	return newick;
}

void Tree::initializeCalibrations(void) {

    for (std::vector<Node*>::iterator it=nodes.begin(); it != nodes.end(); it++)
        {
        if ((*it)->isFossil() == true)
            {
            FossilCalibration* fc = new FossilCalibration;
            calibrations.push_back(fc);
            
            Node* p = (*it);
            while (p->getIsInReconstructedTree() == false && p != root)
                p = p->getAncestor();
            
            fc->setFossilName( (*it)->getName() );
            fc->setTime( (*it)->getTime() );
            std::vector<std::string> desNames = p->getTaxonBipartition();
            for (std::vector<std::string>::iterator it2 = desNames.begin(); it2 != desNames.end(); it2++)
                fc->addTaxonToClade(*it2);
            }
        }
}

void Tree::initializeDownPassSequence(void) {

    downPassSequence.clear();
    passDown(root);
}

bool Tree::isExtantTaxon(Node* p) {

    for (int i=0; i<extantTaxa.size(); i++)
        {
        if (p == extantTaxa[i])
            return true;
        }
    return false;
}

bool Tree::isFossilTaxon(Node* p) {

    for (int i=0; i<fossilTaxa.size(); i++)
        {
        if (p == fossilTaxa[i])
            return true;
        }
    return false;
}

int Tree::lengthOfLongestName(void) {

    int len = 0;
    for (int i=0; i<nodes.size(); i++)
        {
        std::string s = nodes[i]->getName();
        if (s.length() > len)
            len = (int)s.length();
        }
    return len;
}

void Tree::listCalibrations(void) {

    for (std::vector<Node*>::iterator it=nodes.begin(); it != nodes.end(); it++)
        {
        if ((*it)->isFossil() == true)
            {
            Node* p = (*it);
            while (p->getIsInReconstructedTree() == false && p != root)
                p = p->getAncestor();
            
            
            std::cout << (*it)->getName() << " ";
            std::cout << (*it)->getIndex() << " ";
            std::cout << (*it)->getTime() << " ";
            std::cout << "Calibrates node " << p->getIndex() << " ";
            std::vector<std::string> desNames = p->getTaxonBipartition();
            std::cout << "( ";
            for (std::vector<std::string>::iterator it2 = desNames.begin(); it2 != desNames.end(); it2++)
                std::cout << *it2 << " ";
            std::cout << ") ";
            std::cout << std::endl;
            }
        }    for (std::vector<Node*>::iterator it=nodes.begin(); it != nodes.end(); it++)
        {
        if ((*it)->isFossil() == true)
            {
            Node* p = (*it);
            while (p->getIsInReconstructedTree() == false && p != root)
                p = p->getAncestor();
            
            
            std::cout << (*it)->getName() << " ";
            std::cout << (*it)->getIndex() << " ";
            std::cout << (*it)->getTime() << " ";
            std::cout << "Calibrates node " << p->getIndex() << " ";
            std::vector<std::string> desNames = p->getTaxonBipartition();
            std::cout << "( ";
            for (std::vector<std::string>::iterator it2 = desNames.begin(); it2 != desNames.end(); it2++)
                std::cout << *it2 << " ";
            std::cout << ") ";
            std::cout << std::endl;
            }
        }

}

void Tree::markReconstructedTree(void) {

    for (int i=0; i<nodes.size(); i++)
        {
        nodes[i]->setFlag(false);
        nodes[i]->setIsInReconstructedTree(false);
        }
    
    for (int i=0; i<extantTaxa.size(); i++)
        {
        Node* p = extantTaxa[i];
        while (true)
            {
            if (p->getFlag() == false)
                {
                p->setFlag(true);
                if (p == root)
                    break;
                else
                    p = p->getAncestor();
                }
            else
                {
                p->setIsInReconstructedTree(true);
                break;
                }
            }
        }
    
    for (int i=0; i<downPassSequence.size(); i++)
        {
         Node* p = downPassSequence[i];
         
        std::vector<Node*> des = p->getDescendants();
        if (des.size() == 0)
            {
            if (isExtantTaxon(p) == true)
                p->addTaxonToBipartition(p->getName());
            }
        else
            {
            for (std::vector<Node*>::iterator it = des.begin(); it != des.end(); it++)
                {
                std::vector<std::string> desNames = (*it)->getTaxonBipartition();
                for (std::vector<std::string>::iterator it2 = desNames.begin(); it2 != desNames.end(); it2++)
                    p->addTaxonToBipartition(*it2);
                }
            }
        }
}

void Tree::passDown(Node* p) {

    if (p != NULL)
        {
        std::vector<Node*> ndeDescendants = p->getDescendants();
        for (std::vector<Node*>::iterator it = ndeDescendants.begin(); it != ndeDescendants.end(); it++)
            passDown(*it);
            
            p->setIndex((int)downPassSequence.size());
        downPassSequence.push_back(p);
        }
}

Node* Tree::mrcaOfExtantTaxa(void) {

    std::vector<bool> flags(nodes.size());
    for (size_t i=0; i<flags.size(); i++)
        flags[i] = false;
    
    for (size_t i=0; i<extantTaxa.size(); i++)
        {
        Node* p = extantTaxa[i];
        do
            {
            flags[p->getIndex()] = true;
            p = p->getAncestor();
            } while (p != NULL);
        }
    
    Node* mrca = NULL;
    for (size_t i=0; i<getNumberOfDownPassNodes(); i++)
        {
        Node* p = getDownPassNode(i);
        std::vector<Node*> ndeDescendants = p->getDescendants();
        int numFlagged = 0;
        for (std::vector<Node*>::iterator it = ndeDescendants.begin(); it != ndeDescendants.end(); it++)
            {
            if (flags[(*it)->getIndex()] == true)
                numFlagged++;
            }
        if (numFlagged == 2)
            mrca = p;
        }

    return mrca;
}

void Tree::printTree(void) {

    showTree(root, 3);
}

void Tree::pruneToReconstructedProcess(void) {
    
    for (size_t i=0; i<nodes.size(); ++i)
        {
        Node *curr_node = nodes[i];
        if ( curr_node->isExtinct() )
            {
            // let's prune this one
            
            // get the parent node
            Node *anc = curr_node->getAncestor();
            Node *sibling = NULL;
            
            // the ancestor might be
            if ( anc->getNumDescendants() > 1 )
                {
                sibling = anc->getDescendant( 0 );
                if ( sibling == curr_node )
                    {
                    sibling = anc->getDescendant( 1 );
                    }
                
                Node *great_anc = anc->getAncestor();
                great_anc->removeDescendant(anc);
                
                great_anc->addDescendant( sibling );
                sibling->setAncestor( great_anc );
                
                // remove the nodes from the list
                nodes.erase( nodes.begin() + i );
                nodes.erase(std::remove(nodes.begin(), nodes.end(), anc), nodes.end());
                }
            else
                {
                // just remove this node
                anc->removeDescendant(curr_node);
                
                // remove the nodes from the list
                nodes.erase( nodes.begin() + i );
                }
            
            // this is not efficient but for safety
            i = 0;
            }
        }
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

void Tree::showTree(Node* p, int indent) {

	if (p != NULL)
		{
        std::vector<Node*> ndeDescendants = p->getDescendants();
		for (int i=0; i<indent; i++)
			std::cout << " ";
        if (p->getIsInReconstructedTree() == true)
            std::cout << "R";
        else
            std::cout << " ";
        std::cout << p->getIndex();
        std::cout << " ( ";
        for (std::vector<Node*>::iterator it = ndeDescendants.begin(); it != ndeDescendants.end(); it++)
            std::cout << (*it)->getIndex() << " ";
        std::cout << ") " << p->getTime();
        if (p->getNumDescendants() == 0 || p->getNumDescendants() == 1)
            std::cout << " (" << p->getName() << ") ";
		if (p == root)
			std::cout << " <- Root ";
            
        std::cout << std::endl;
        for (std::vector<Node*>::iterator it = ndeDescendants.begin(); it != ndeDescendants.end(); it++)
            showTree(*it, indent+2);
		}
}

void Tree::writeTree(Node* p, std::stringstream& ss) {

	if (p != NULL)
		{
		if (p->getNumDescendants() == 0)
            {
            double bl = (p->getAncestor() == NULL ? 0 : p->getTime() - p->getAncestor()->getTime());
			ss << p->getName() << ":" << std::fixed << std::setprecision(2) << bl;
			}
        else if (p->getNumDescendants() == 1)
            {
            ss << "(";
            writeTree(p->getDescendants()[0], ss);
            double bl = (p->getAncestor() == NULL ? 0 : p->getTime() - p->getAncestor()->getTime());
            ss << ")" << p->getName() << ":" << std::fixed << std::setprecision(2) << bl;
            }
		else
			{
            ss << "(";
            std::vector<Node*> myDescendants = p->getDescendants();
            for (int i=0; i<myDescendants.size(); i++)
                {
                writeTree(myDescendants[i], ss);
                if (i + 1 != myDescendants.size())
                    ss << ",";
                }
            double bl = (p->getAncestor() == NULL ? 0 : p->getTime() - p->getAncestor()->getTime());
            ss << "):" << std::fixed << std::setprecision(2) << bl;
            }
		}
}





