#ifndef _hgcssevent_hh_
#define _hgcssevent_hh_
#include <iomanip>
#include <vector>
#include "Rtypes.h"
#include <sstream>
#include <map>
#include "TVector3.h"
class HGCSSEvent {

public:
	HGCSSEvent() :
			event_(0) {

	}
	;

	~HGCSSEvent() {
	}
	;

	inline unsigned eventNumber() const {
		return event_;
	}
	;

	inline void eventNumber(const unsigned aNum) {
		event_ = aNum;
	}
	;

	inline double vtx_x() const {
		return xvtx_;
	}
	;
	inline double steelThick() const {
		return steelThick_;
	}
	;
	inline void vtx_x(const double x) {
		xvtx_ = x;
	}
	;

	inline double vtx_y() const {
		return yvtx_;
	}
	;

	inline TVector3 seeds() const {
		return seeds_;
	}
	;
	inline TVector3 status() const {
		return status_;
	}
	;
	inline double dep() const {
		return dep_;
	}
	;
	inline double wgtDep() {
		return wgtDep_;
	}
	;
	inline void vtx_y(const double y) {
		yvtx_ = y;
	}
	;

	inline double vtx_z() const {
		return zvtx_;
	}
	;

	inline void vtx_z(const double z) {
		zvtx_ = z;
	}
	;
	inline void steelThick(const double thick) {
		steelThick_ = thick;
	}
	;
	inline void seeds(TVector3 seeds) {
		seeds_ = seeds;
	}
	;
	inline void status(TVector3 status) {
		status_ = status;
	}
	;
	inline void dep(double dep) {
		dep_ = dep;
	}
	;
	inline void wgtDep(double wgtDep) {
		wgtDep_ = wgtDep;
	}
	;
private:

	unsigned event_;
	double xvtx_;
	double yvtx_;
	double zvtx_;
	double steelThick_;
	double dep_;
	double wgtDep_;
	TVector3 seeds_,status_;

ClassDef(HGCSSEvent,1)
	;

};

#endif
