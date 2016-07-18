#include "G4VPhysicalVolume.hh"

#include "SamplingSection.hh"

//
std::pair<G4bool,G4bool> SamplingSection::add(G4double depositRawE,  G4VPhysicalVolume* vol) {
	std::string lstr = vol->GetName();
	bool breakSwitch = false;
	bool isSens = false;
	for (unsigned ie(0); ie < n_elements * n_sectors; ++ie) {
			if (breakSwitch) return std::make_pair(breakSwitch,isSens);
		if (sublayer_vol[ie] && lstr == sublayer_vol[ie]->GetName()) {
			breakSwitch = true;
			unsigned idx = getSensitiveLayerIndex(lstr);
			unsigned eleidx = ie % n_elements;
			isSens = isSensitiveElement(eleidx);
			sublayer_RawDep[eleidx] += depositRawE;

			} //if in right material
		} //loop on available materials
	return std::make_pair(breakSwitch,isSens);
	}

//
void SamplingSection::report(bool header) {
	if (header)
		G4cout
				<< "E/[MeV]\t  Si\tAbsorber\tTotal\tSi g frac\tSi e frac\tSi mu frac\tSi had frac\tSi <t> \t nG4SiHits"
				<< G4endl;
		G4cout << std::setprecision(3) << "\t  " << getTotalSensE()
		<< G4endl;
	}


G4double SamplingSection::getTotalSensE() {
	double etot = 0;
	for (unsigned ie(0); ie < n_elements; ++ie) {
		if (isSensitiveElement(ie))
			etot += sublayer_RawDep[ie];
	}
	return etot;
}
