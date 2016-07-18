#ifndef HGCSSPUenergy_h
#define HGCSSPUenergy_h

#include <string>
#include <fstream>
#include <vector>
#include "TF1.h"
#include "TMath.h"

class HGCSSPUenergy {

public:
	HGCSSPUenergy() {
	}
	;
	HGCSSPUenergy(const std::string filePath);
	~HGCSSPUenergy();
	double getDensity(const double & eta, const unsigned layer,
			const double & cellSize, const unsigned PU) const;

private:
	std::vector<double> p0_;
	std::vector<double> p1_;
};

#endif 
