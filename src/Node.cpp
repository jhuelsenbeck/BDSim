#include "Node.h"


Node::Node(void) {

    index    = 0;
    ancestor = NULL;
    time     = 0.0;
    name     = "";
    fossil   = false;
    extinct  = false;
    
}

void Node::addDescendant(Node* p) {

    descendants.push_back(p);
}

Node* Node::getDescendant(int idx) {

    if (idx < 0 || idx >= descendants.size())
        return NULL;
    return descendants[idx];
}

int Node::getNumDescendants(void) {

    return (int)descendants.size();
}

void Node::removeDescendant(Node* d)
{
    
    for (std::vector<Node*>::iterator it = descendants.begin(); it!=descendants.end(); ++it)
    {
        if ( *it == d )
        {
            descendants.erase(it);
            break;
        }
        
    }
     
}

void Node::removeDescendants(void) {

    descendants.clear();
}