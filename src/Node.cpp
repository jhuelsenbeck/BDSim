#include "Node.h"


Node::Node(void) {

    index    = 0;
    ancestor = NULL;
    time     = 0.0;
}

void Node::addDescendant(Node* p) {

    descendants.push_back(p);
}
