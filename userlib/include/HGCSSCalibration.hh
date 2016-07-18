#ifndef HGCSSCalibration_h
#define HGCSSCalibration_h

#include <string>
#include <vector>
#include "TH2D.h"

#include "HGCSSDetector.hh"

class HGCSSCalibration {

public:
	HGCSSCalibration() {
		vtx_x_ = 0;
		vtx_y_ = 0;
		vtx_z_ = 0;
	}
	;

	HGCSSCalibration(std::string filePath, const bool bypassR = false,
			const unsigned nSi = 2);

	~HGCSSCalibration();

	inline void setVertex(const double & x, const double & y,
			const double & z) {

		vtx_x_ = x;
		vtx_y_ = y;
		vtx_z_ = z;

	}
	;

	double addTimeOfFlight(const double & aTime, const double & posx,
			const double & posy, const double & posz);

	double addTimeOfFlight(const double & aTime, const double & posx,
			const double & posy, const double & posz, const double & vtxx,
			const double & vtxy, const double & vtxz);

	double correctTime(const double & aTime, const double & posx,
			const double & posy, const double & posz);

	double correctTime(const double & aTime, const double & posx,
			const double & posy, const double & posz, const double & vtxx,
			const double & vtxy, const double & vtxz);

	double MeVToMip(const unsigned layer, const bool absWeight = false) const;

	//double MeVToMip(const unsigned layer, const double aEta,
	//const bool absWeight=false) const;

	double MeVToMip(const unsigned layer, const double aRadius,
			const bool absWeight = false) const;

private:

	double vtx_x_;
	double vtx_y_;
	double vtx_z_;
	bool bypassRadius_;
	unsigned nSiLayers_;
};

#endif

