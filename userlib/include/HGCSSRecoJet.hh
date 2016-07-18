#ifndef _hgcssrecojet_hh_
#define _hgcssrecojet_hh_

#include <iomanip>
#include <vector>
#include "Rtypes.h"
#include <sstream>

class HGCSSRecoJet {

public:
	HGCSSRecoJet() :
			energy_(0), px_(0), py_(0), pz_(0), nConstituents_(0), area_(0), area_error_(
					0) {
	}
	;

	~HGCSSRecoJet() {
	}
	;

	HGCSSRecoJet(const double & px, const double & py, const double & pz,
			const double & E);

	inline double energy() const {
		return energy_;
	}
	;

	inline double E() const {
		return energy_;
	}
	;

	inline void energy(const double & energy) {
		energy_ = energy;
	}
	;

	inline double px() const {
		return px_;
	}
	;

	inline double py() const {
		return py_;
	}
	;

	inline double pz() const {
		return pz_;
	}
	;

	inline double area() const {
		return area_;
	}
	;

	inline void area(const double & aArea) {
		area_ = aArea;
	}
	;

	inline double area_error() const {
		return area_error_;
	}
	;

	inline void area_error(const double & aErr) {
		area_error_ = aErr;
	}
	;

	inline unsigned nConstituents() const {
		return nConstituents_;
	}
	;

	inline void nConstituents(const double n) {
		nConstituents_ = n;
	}
	;

	void Print(std::ostream & aOs) const;

private:

	double energy_;
	double px_;
	double py_;
	double pz_;
	unsigned nConstituents_;
	double area_;
	double area_error_;

ClassDef(HGCSSRecoJet,1)
	;

};

typedef std::vector<HGCSSRecoJet> HGCSSRecoJetVec;

#endif
