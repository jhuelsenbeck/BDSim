//
//  RevScript.hpp
//  BDSim
//
//  Created by Sebastian Hoehna on 10/9/15.
//  Copyright © 2015 Sebastian Hoehna. All rights reserved.
//

#ifndef RevScript_H
#define RevScript_H

#include <iostream>
#include <string>


class Tree;

class RevScript {
    
    public:

        enum TREE_PRIOR { BDP, FBDP };
        enum CLOCK_TYPE { STRICT, UCLN };

    
        RevScript(const std::string &pa, const std::string &f, int r, Tree *t, TREE_PRIOR p, bool mo);
    
    
        void            print( void );
    


    private:
    
//        void            printQuit(std::ostream &o);
        void            printBDP(std::ostream &o);
        void            printFBDP(std::ostream &o);
        void            printInitialSetting(std::ostream &o);
        void            printMcmc(std::ostream &o);
        void            printModel(std::ostream &o);
        void            printMolecularSubstitutionModel(std::ostream &o);
        void            printMorphologicalSubstitutionModel(std::ostream &o);
        void            printMonitors(std::ostream &o);
        void            printQuit(std::ostream &o);
        void            printReadMolecularData(std::ostream &o);
        void            printReadMorphologicalData(std::ostream &o);
        void            printReadTree(std::ostream &o);
        void            printRelativeMorphologicalClock(std::ostream &o);
        void            printSetting(std::ostream &o);
        void            printStrictMolecularClock(std::ostream &o);
        void            printSummarizeTrees(std::ostream &o);
        void            printUclnMolecularClock(std::ostream &o);
    
    
        std::string     base_file_name;
        std::string     base_file_path;
        int             rep;
        Tree*           tree;
        TREE_PRIOR      prior;
        CLOCK_TYPE      clock;
        bool            mol_only;

};


#endif /* RevScript_hpp */
