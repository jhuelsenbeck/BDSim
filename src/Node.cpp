#include "Node.h"


Node::Node(void) {

    index    = 0;
    ancestor = NULL;
    time     = 0.0;
    name     = "";
}

void Node::addDescendant(Node* p) {

    descendants.push_back(p);
}

int Node::getNumDescendants(void) {

    return (int)descendants.size();
}

void Node::removeDescendants(void) {

    descendants.clear();
}