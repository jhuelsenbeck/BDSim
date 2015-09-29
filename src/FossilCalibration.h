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
        double                      getTime(void) { return time; }
        void                        setFossilName(std::string s) { fossilName = s; }
        void                        setTime(double x) { time = x; }

    private:
        std::string                 fossilName;
        double                      time;
        std::vector<std::string>    cladeNames;
};


#endif