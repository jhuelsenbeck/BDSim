#include <iostream>
#include "Event.h"
#include "EventHistory.h"
#include "MbRandom.h"
#include "Node.h"
#include "Tree.h"



EventHistory::EventHistory(Tree* t, MbRandom* rp, double a0, double er) {

    // set values
    myTree    = t;
    ranPtr    = rp;
    alpha0    = a0;
    eventRate = er;
    treeLength = 0.0;
    for (int n=0; n<myTree->getNumberOfDownPassNodes(); n++)
        {
        Node * p = myTree->getDownPassNode(n);
        if (myTree->isRoot(p) == false)
            treeLength += p->getTime() - p->getAncestor()->getTime();
        }
    
    // initialize events
    int numEvents = ranPtr->poissonRv(treeLength*eventRate);
    for (int i=0; i<numEvents; i++)
        addEvent();
}

EventHistory::~EventHistory(void) {

    for (std::set<Event*,CompEvents>::iterator it = events.begin(); it != events.end(); it++)
        delete *it;
}

void EventHistory::addEvent(void) {

    // pick a branch at random (and proportional to length)
    double u = ranPtr->uniformRv() * treeLength;
    double sum = 0.0, pos = 0.0;
    Node* whichBranch = NULL;
    for (int n=0; n<myTree->getNumberOfDownPassNodes(); n++)
        {
        Node * p = myTree->getDownPassNode(n);
        if (myTree->isRoot(p) == false)
            {
            sum += p->getTime() - p->getAncestor()->getTime();
            if (u < sum)
                {
                whichBranch = p;
                pos = (sum - u) + p->getAncestor()->getTime();
                break;
                }
            }
        }
    
    if (whichBranch != NULL)
        {
        //double r = ranPtr->gammaRv(alpha0, alpha0) / exp(psi(alpha0));
        double r = RndGamma(alpha0) / exp(psi(alpha0));
        Event* newE = new Event(whichBranch, pos, r);
        events.insert(newE);
        }
    else
        {
        std::cout << "ERROR: Can't find a branch for the event" << std::endl;
        exit(1);
        }
}

Event* EventHistory::getEvent(int branchIdx, int idx) {

    Event* evnt = NULL;
    int num = 0;
    for (std::set<Event*,CompEvents>::iterator it = events.begin(); it != events.end(); it++)
        {
        if ( (*it)->getNode()->getIndex() == branchIdx )
            {
            if (num == idx)
                return (*it);
            num++;
            }
        }
    
    return evnt;
}

int EventHistory::numEventsForBranchWithIndex(int idx) {

    int num = 0;
    for (std::set<Event*,CompEvents>::iterator it = events.begin(); it != events.end(); it++)
        {
        if ( (*it)->getNode()->getIndex() == idx )
            num++;
        }
    return num;
}

void EventHistory::print(void) {

    int i = 0;
    for (std::set<Event*,CompEvents>::iterator it = events.begin(); it != events.end(); it++)
        {
        std::cout << ++i << " -- ";
        (*it)->print();
        std::cout << std::endl;
        }
}

double EventHistory::psi(double x) {

	/* Algorithm AS 103 in Applied Statistics (1976) Vol. 25, No. 3

	Fortran code from netlib (http://www.netlib.no)
	Also at StatLib (http://www.stat.cmu.edu)
	Translated to C by Bret Larget on August 6, 1999

	psi(x) = d/dx [Gamma(x)] = Gamma'(x) / Gamma(x)

	x must be positive.
	psi does not check and assumes it is only called properly
	*/

	const double half = 0.5, one = 1.0;
	const double s = 1.0e-05;
	const double c = 8.5;
	const double s3 = 8.3333333333e-02;  /* 1/12 = 3rd Stirling coefficient */
	const double s4 = 8.3333333333e-03;  /* 1/120 = 4th Stirling coefficient */
	const double s5 = 3.9682539683e-03;  /* 1/252 = 5th Stirling coefficient */
	const double d1 = -0.5772156649;     /* - Euler's constant = psi(1) */
	double ans = 0.0;
	double r;

	/* approximation if x <= s */

	if (x <= s)
		return (d1 - one/x);

	/* reduce to case where psi(x+n) where (x+n) >= c */

	while (x < c)
		{
		ans -= one/x;
		x += one;
		}

	r = one/x;
	ans += log(x) - half*r;
	r *= r;
	ans -= r*(s3-r*(s4-r*s5));
	return (ans);
}

double EventHistory::RndGamma(double s) {

	double  r=0.0;
	
	if (s <= 0.0)      
		puts ("jgl gamma..");
	else if (s < 1.0)  
		r = RndGamma1(s);
	else if (s > 1.0)  
		r = RndGamma2(s);
	else           
		r =- log(ranPtr->uniformRv());
	return (r);
   
}

double EventHistory::RndGamma1(double s) {

	double			r, x=0.0, small=1e-37, w;
	static double   a, p, uf, ss=10.0, d;
	
	if (s!=ss) 
		{
		a  = 1.0-s;
		p  = a/(a+s*exp(-a));
		uf = p*pow(small/a,s);
		d  = a*log(a);
		ss = s;
		}
	for (;;) 
		{
		r = ranPtr->uniformRv();
		if (r > p)        
			x = a-log((1.0-r)/(1.0-p)), w=a*log(x)-d;
		else if (r>uf)  
			x = a*pow(r/p,1/s), w=x;
		else            
			return (0.0);
		r = ranPtr->uniformRv();
		if (1.0-r <= w && r > 0.0)
		if (r*(w+1.0) >= 1.0 || -log(r) <= w)  
			continue;
		break;
		}
	return (x);
   
}

double EventHistory::RndGamma2 (double s) {

	double			r ,d, f, g, x;
	static double	b, h, ss=0;
	
	if (s!=ss) 
		{
		b  = s-1.0;
		h  = sqrt(3.0*s-0.75);
		ss = s;
		}
	for (;;) 
		{
		r = ranPtr->uniformRv();
		g = r-r*r;
		f = (r-0.5)*h/sqrt(g);
		x = b+f;
		if (x <= 0.0) 
			continue;
		r = ranPtr->uniformRv();
		d = 64*r*r*g*g*g;
		if (d*x < x-2.0*f*f || log(d) < 2*(b*log(x/b)-f))  
			break;
		}
	return (x);
   
}
