#include <iomanip>
#include <iostream>
#include "Event.h"
#include "Node.h"



Event::Event(Node* n, double t, double f) {

    myNode = n;
    branchPos = t;
    rateFactor = f;
}

bool Event::operator<(const Event& e) const {

    if ( myNode->getIndex() < e.myNode->getIndex() )
        {
        return true;
        }
    else if (myNode->getIndex() == e.myNode->getIndex())
        {
        if (branchPos < e.branchPos)
            return true;
        }
    return false;
}

void Event::print(void) {

    std::cout << "(" << myNode->getIndex() << "," << branchPos << "," << rateFactor << ")";
}