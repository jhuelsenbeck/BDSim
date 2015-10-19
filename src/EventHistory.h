#ifndef EventHistory_H
#define EventHistory_H

#include "Event.h"
#include <set>

class MbRandom;
class Node;
class Tree;

class EventHistory {

	public:
                                    EventHistory(Tree* t, MbRandom* rp, double a0, double er);
                                   ~EventHistory(void);
        void                        addEvent(void);
        int                         numEvents(void) { return (int)events.size(); }
        int                         numEventsForBranchWithIndex(int idx);
        Event*                      getEvent(int branchIdx, int idx);
        void                        print(void);

    private:
                                    EventHistory(void) { }
        double                      RndGamma(double s);
        double                      RndGamma1(double s);
        double                      RndGamma2(double s);
        double                      psi(double x);
        Tree*                       myTree;
        MbRandom*                   ranPtr;
        double                      alpha0;
        double                      eventRate;
        double                      treeLength;
        std::set<Event*,CompEvents> events;
};

#endif