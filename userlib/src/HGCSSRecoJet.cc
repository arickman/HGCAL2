#include "HGCSSRecoJet.hh"

#include <iomanip>
#include <cmath>
#include <stdlib.h>

HGCSSRecoJet::HGCSSRecoJet(const double & px, const double & py,
		const double & pz, const double & E) {
	energy_ = E;
	px_ = px;
	py_ = py;
	pz_ = pz;
}

void HGCSSRecoJet::Print(std::ostream & aOs) const {
	aOs << " -----------------------------------" << std::endl
			<< " --------- Printing jet: -----------" << std::endl << " -- E = "
			<< energy_ << std::endl << " -- px = " << px_ << std::endl
			<< " -- py = " << py_ << std::endl << " -- pz = " << pz_
			<< std::endl << " -- nConstituents = " << nConstituents_
			<< std::endl << " -- area = " << area_ << " +/- " << area_error_
			<< std::endl << " -----------------------------------" << std::endl;
}
