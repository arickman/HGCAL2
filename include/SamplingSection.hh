#ifndef _samplingsection_hh_
#define _samplingsection_hh_

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Colour.hh"

#include <iomanip>
#include <vector>

#include "G4SiHit.hh"

class SamplingSection {
public:
	//CTOR
	SamplingSection(std::vector<std::pair <G4double,std::string>> iEle) {

			std::vector<G4double> aThicknessVec;std::vector<std::string> aMaterialVec;
				for (unsigned i = 0; i < iEle.size(); i++)
				{
					aThicknessVec.push_back(iEle.at(i).first);
					aMaterialVec.push_back(iEle.at(i).second);

				}
				Total_thick = 0;
				n_sens_elements=0;
				n_elements=0;
				n_sectors=0;
				sublayer_thick.clear();
				ele_name.clear();
				sublayer_X0.clear();
				sublayer_L0.clear();
				sublayer_vol.clear();
				hasScintillator = false;
				for (unsigned ie(0); ie<aThicknessVec.size(); ++ie) {
					//consider only material with some non-0 width...
					if (aThicknessVec[ie]>0) {
						sublayer_thick.push_back(aThicknessVec[ie]);
						ele_name.push_back(aMaterialVec[ie]);
						if (aMaterialVec[ie]== "Scintillator") hasScintillator = true;
						sublayer_X0.push_back(0);
						sublayer_dEdx.push_back(0);
						sublayer_L0.push_back(0);
						sublayer_vol.push_back(0);
						Total_thick+=aThicknessVec[ie];
						++n_elements;
						//the following method check the total size...
						//so incrementing first.
						if (isSensitiveElement(n_elements-1)) {
							++n_sens_elements;
						}
					}
				}

				resetCounters();

				std::cout << " -- End of sampling section initialisation. Input " << aThicknessVec.size() << " elements, constructing " << n_elements << " elements with " << n_sens_elements << " sensitive elements." << std::endl;

			};

			//DTOR
			~SamplingSection() {};

			inline void setNumberOfSectors(const unsigned nSec) {
				n_sectors = nSec;
				sublayer_vol.clear();
				for (unsigned ie(0); ie<n_elements*n_sectors; ++ie) {
					sublayer_vol.push_back(0);
				}
			};

			//
			std::pair<G4bool,G4bool> add(G4double depositRawE,G4VPhysicalVolume* vol);

			inline bool isSensitiveElement(const unsigned & aEle) {
				if (aEle < n_elements &&
						(ele_name[aEle] == "Si" || ele_name[aEle] == "Scintillator")
				) return true;
				return false;
			};

			inline unsigned getSensitiveLayerIndex(std::string astr) {
				if (astr.find("_")== astr.npos) return 0;
				size_t pos = astr.find("phys");
				if (pos != astr.npos && pos>1) {
					unsigned idx = 0;
					std::istringstream(astr.substr(pos-1,1))>>idx;
					return idx;
				}
				return 0;
			};

			inline G4Colour g4Colour(const unsigned & aEle) {
				if (isSensitiveElement(aEle)) return G4Colour::Red();
				if (ele_name[aEle] == "Cu") return G4Colour::Black();
				if (isAbsorberElement(aEle)) return G4Colour::Gray();
				if (ele_name[aEle] == "PCB") return G4Colour::Blue();
				if (ele_name[aEle] == "Air") return G4Colour::Cyan();
				return G4Colour::Yellow();
			};

			inline bool isAbsorberElement(const unsigned & aEle) {
				if (aEle < n_elements &&
						(
								ele_name[aEle] == "Pb" || ele_name[aEle] == "Cu" ||
								ele_name[aEle] == "W" || ele_name[aEle] == "Brass" ||
								ele_name[aEle] == "Fe" || ele_name[aEle] == "Steel" ||
								ele_name[aEle] == "SSteel" || ele_name[aEle] == "Al" ||
								ele_name[aEle] == "WCu" || ele_name[aEle] == "NeutMod"
						)
				) return true;
				return false;
			};

			//reset
			inline void resetCounters()
			{
				sublayer_RawDep.clear();
				sublayer_RawDep.resize(n_elements,0);

			}



			G4double getTotalSensE();

			//
			void report(bool header=false);

			//members
			unsigned n_elements;
			unsigned n_sectors;
			unsigned n_sens_elements;
			std::vector<G4double> sublayer_thick;
			std::vector<std::string> ele_name;
			std::vector<G4double> sublayer_X0;
			std::vector<G4double> sublayer_dEdx;
			std::vector<G4double> sublayer_L0;
			std::vector<G4double> sublayer_RawDep;
			std::vector<G4double> sublayer_PrimaryDep;
			std::vector<G4double> sublayer_dl;
			std::vector<G4VPhysicalVolume*> sublayer_vol;
			G4double Total_thick;
			bool hasScintillator;

		};

#endif
