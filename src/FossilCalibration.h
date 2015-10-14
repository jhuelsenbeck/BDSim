#ifndef FossilCalibration_H
#define FossilCalibration_H

#include <string>
#include <vector>

class FossilCalibration {

	public:
                                    FossilCalibration(void);
        void                        addTaxonToClade(std::string s) { cladeNames.push_back(s); }
        std::string                 getFossilName(void) { return fossilName; }
        int                         getNumTaxaInClade(void) { return (int)cladeNames.size(); }
        std::string                 getCladeTaxonIndexed(int idx) { return cladeNames[idx]; }
        double                      getFossilTime(void) { return fossilTime; }
        double                      getNodeTime(void) { return nodeTime; }
        void                        setFossilName(std::string s) { fossilName = s; }
        void                        setFossilTime(double x) { fossilTime = x; }
        void                        setNodeTime(double x) { nodeTime = x; }

    private:
        std::string                 fossilName;
        double                      fossilTime;
        double                      nodeTime;
        std::vector<std::string>    cladeNames;
};


#endif