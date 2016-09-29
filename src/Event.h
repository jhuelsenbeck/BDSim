#ifndef Event_H
#define Event_H

class Node;

class Event {

	public:
                            Event(Node* n, double t, double f);
        bool                operator<(const Event& e) const;
        Node*               getNode(void) { return myNode; }
        double              getBranchPos(void) { return branchPos; }
        double              getRateFactor(void) { return rateFactor; }
        void                print(void);
        void                setNode(Node* p) { myNode = p; }
        void                setBranchPos(double x) { branchPos = x; }
        void                setRateFactor(double x) { rateFactor = x; }

    private:
                            Event(void) { }
        Node*               myNode;
        double              branchPos;
        double              rateFactor;
};

class CompEvents {
    
    public:
        bool                operator()(Event* e1, Event* e2) const { return (*e1 < *e2); }
};

#endif