#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include <boost/algorithm/string.hpp>

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Polyhedra.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4PhysicalConstants.hh"

using namespace std;

//

DetectorConstruction::DetectorConstruction(G4int ver, G4int mod, double steelThick) :
		version_(ver), model_(mod),steelThick_(steelThick) {

		G4cout << "[DetectorConstruction] starting v_HGCALSYM_v1" << G4endl;
		buildTracker();
		buildECal();
		buildHCal(steelThick);
		//Add the target



	DefineMaterials();
	SetMagField(0);
	m_detectorMessenger = new DetectorMessenger(this);
	UpdateCalorSize();
}
void DetectorConstruction::buildTracker(){

		std::vector<std::pair <G4double,std::string>> iEleL;
		initLayer(1);

		if (version_ == HGCAL_E26_TH || version_ == HGCAL_E26_T || version_ == HGCAL_E40_TH || version_ == HGCAL_E40_T || version_ == T){
			for (int i = 0; i < 6; i ++){
				iEleL.push_back(make_pair(.7*mm,"Si"));
				iEleL.push_back(make_pair(99.3*mm,"G4_Galactic"));
			}

			iEleL.push_back(make_pair(.7*mm,"Si"));
			iEleL.push_back(make_pair(6.8*mm,"G4_Galactic"));

			iEleL.push_back(make_pair(.3504*mm,"W"));
			iEleL.push_back(make_pair(6.8*mm,"G4_Galactic"));


			for (int i = 0; i < 3; i ++){
				iEleL.push_back(make_pair(.7*mm,"Si"));
				if (i < 3){
					iEleL.push_back(make_pair(15*mm,"G4_Galactic"));

				}
			}

			iEleL.push_back(make_pair(.7*mm,"Si"));
			iEleL.push_back(make_pair(44.8*mm,"G4_Galactic"));
			iEleL.push_back(make_pair(.7*mm,"Si"));
			iEleL.push_back(make_pair(84.3*mm,"G4_Galactic"));
			iEleL.push_back(make_pair(.7*mm,"Si"));
			iEleL.push_back(make_pair(14.3*mm,"G4_Galactic"));


			m_caloStruct.push_back( SamplingSection(iEleL) );
		}
		else if (version_ == HGCAL_E26 || version_ == HGCAL_E40){
			iEleL.push_back(make_pair(.3504*mm,"W"));
			iEleL.push_back(make_pair(20*cm,"G4_Galactic"));
		}

}

void DetectorConstruction::buildECal(){

		std::vector<std::pair <G4double,std::string>> iEleL,iEleR;
		initLayer(1);
		if (version_ == HGCAL_E26_TH || version_ == HGCAL_E26_T || version_ == HGCAL_E26_H || version_ == HGCAL_E26){
			G4double airThick = 2*mm,
					pcbThick = 2*mm,
					wThick = 2.*mm,
					wcuThick = 0.6*mm;

			iEleL.clear();
			iEleL.push_back(make_pair(0.5*mm,"Cu"));
			iEleL.push_back(make_pair(0.5*mm,"CFMix"));
			iEleL.push_back(make_pair(wThick,"W"));

			iEleL.push_back(make_pair(0.5*mm,"CFMix"));
			iEleL.push_back(make_pair(0.5*mm,"Cu"));
			iEleL.push_back(make_pair(airThick,"Air"));
			iEleL.push_back(make_pair(pcbThick,"PCB"));
			iEleL.push_back(make_pair(0.1*mm,"Si"));
			iEleL.push_back(make_pair(0.1*mm,"Si"));
			iEleL.push_back(make_pair(0.1*mm,"Si"));

			std::vector<std::pair <G4double,std::string>> iEleR;
			iEleR.push_back(make_pair(wcuThick,"WCu"));
			iEleR.push_back(make_pair(6*mm,"Cu"));
			iEleR.push_back(make_pair(wcuThick,"WCu"));
			iEleR.push_back(make_pair(0.1*mm,"Si"));
			iEleR.push_back(make_pair(0.1*mm,"Si"));
			iEleR.push_back(make_pair(0.1*mm,"Si"));
			iEleR.push_back(make_pair(pcbThick,"PCB"));
			iEleR.push_back(make_pair(airThick,"Air"));

				unsigned Nmodule=4;
				for(unsigned i=0; i<Nmodule; i++) {
					m_caloStruct.push_back( SamplingSection(iEleL) );
					m_caloStruct.push_back( SamplingSection(iEleR) );
				}

				Nmodule=5;
				iEleL[2].first = 2.8*mm;
				iEleR[0].first = 1.2*mm;
				iEleR[2].first = 1.2*mm;
				for(unsigned i=0; i<Nmodule; i++) {
					m_caloStruct.push_back( SamplingSection(iEleL) );
					m_caloStruct.push_back( SamplingSection(iEleR) );
				}

				Nmodule=4;
				iEleL[2].first = 4.2*mm;
				iEleR[0].first = 2.2*mm;
				iEleR[2].first = 2.2*mm;
				for(unsigned i=0; i<Nmodule; i++) {
					m_caloStruct.push_back( SamplingSection(iEleL) );
					m_caloStruct.push_back( SamplingSection(iEleR) );
				}
		}
		else if (version_ == HGCAL_E40_TH || version_ == HGCAL_E40_T || version_ == HGCAL_E40_H || version_ == HGCAL_E40 || version_ == DUMMY ){
			G4double airThick = 2*mm,
					pcbThick = 2*mm,
					wThick = 2.5*mm,
					wcuThick = 0.6*mm;

			iEleL.clear();
			iEleL.push_back(make_pair(0.5*mm,"Cu"));
			iEleL.push_back(make_pair(0.5*mm,"CFMix"));
			iEleL.push_back(make_pair(wThick,"W"));
			iEleL.push_back(make_pair(0.5*mm,"CFMix"));
			iEleL.push_back(make_pair(0.5*mm,"Cu"));
			iEleL.push_back(make_pair(airThick,"Air"));
			iEleL.push_back(make_pair(pcbThick,"PCB"));
			iEleL.push_back(make_pair(0.1*mm,"Si"));
			iEleL.push_back(make_pair(0.1*mm,"Si"));
			iEleL.push_back(make_pair(0.1*mm,"Si"));


			unsigned Nmodule=10;
			for(unsigned i=0; i<Nmodule; i++) {
				m_caloStruct.push_back( SamplingSection(iEleL) );
			}
			iEleL[2].first = 3.2*mm;
			Nmodule=32;
			for(unsigned i=0; i<Nmodule; i++) {
				m_caloStruct.push_back( SamplingSection(iEleL) );
			}
		}

}

void DetectorConstruction::buildHCal(double steelThick){
	std::vector<std::pair <G4double,std::string>> iEleL;

	if (version_ == H){
		iEleL.push_back(make_pair(80*cm,"G4_Galactic"));
		m_caloStruct.push_back( SamplingSection(iEleL));
		iEleL.clear();
	}

	if (version_ == HGCAL_E26_TH || version_ == HGCAL_E26_H || version_ == HGCAL_E40_TH || version_ == HGCAL_E40_H || version_ == H){



		iEleL.push_back(make_pair(3*mm,"Cu"));
		iEleL.push_back(make_pair(1*mm,"Pb"));
		iEleL.push_back(make_pair(steelThick*mm,"SSteel"));
		iEleL.push_back(make_pair(0.5*mm,"Cu"));
		iEleL.push_back(make_pair(9*mm,"Scintillator"));
		unsigned Nmodule=15;

		for(unsigned i=0; i<Nmodule; i++) {
			m_caloStruct.push_back( SamplingSection(iEleL) );
		}
	}
}
DetectorConstruction::~DetectorConstruction() {
	delete m_detectorMessenger;
}

//
void DetectorConstruction::DefineMaterials() {
	G4NistManager* nistManager = G4NistManager::Instance();

	m_materials["Al"] = nistManager->FindOrBuildMaterial("G4_Al", false);
	m_dEdx["Al"] = 0.4358;
	m_materials["W"] = nistManager->FindOrBuildMaterial("G4_W", false);
	m_dEdx["W"] = 2.210;
	m_materials["Pb"] = nistManager->FindOrBuildMaterial("G4_Pb", false);
	m_dEdx["Pb"] = 1.274;
	m_materials["Cu"] = nistManager->FindOrBuildMaterial("G4_Cu", false);
	m_dEdx["Cu"] = 1.257;
	m_materials["Si"] = nistManager->FindOrBuildMaterial("G4_Si", false);
	m_dEdx["Si"] = 0.3876;
	m_materials["Zn"] = nistManager->FindOrBuildMaterial("G4_Zn", false);
	m_dEdx["Zn"] = 1.007;
	m_materials["Air"] = nistManager->FindOrBuildMaterial("G4_AIR", false);
	m_dEdx["Air"] = 0;
	m_materials["G4_Galactic"] = nistManager->FindOrBuildMaterial("G4_Galactic",
			false);
	m_dEdx["G4_Galactic"] = 0;
	m_materials["Fe"] = nistManager->FindOrBuildMaterial("G4_Fe", false);
	m_dEdx["Fe"] = 1.143;
	m_materials["Mn"] = nistManager->FindOrBuildMaterial("G4_Mn", false);
	m_dEdx["Mn"] = 1.062;
	m_materials["C"] = nistManager->FindOrBuildMaterial("G4_C", false);
	m_dEdx["C"] = 0.3952;
	m_materials["H"] = nistManager->FindOrBuildMaterial("G4_H", false);
	m_dEdx["H"] = 0;
	m_materials["Cl"] = nistManager->FindOrBuildMaterial("G4_Cl", false);
	m_dEdx["Cl"] = 0;
	m_materials["Cr"] = nistManager->FindOrBuildMaterial("G4_Cr", false);
	m_dEdx["Cr"] = 1.046;
	m_materials["Ni"] = nistManager->FindOrBuildMaterial("G4_Ni", false);
	m_dEdx["Ni"] = 1.307;
	m_materials["O"] = nistManager->FindOrBuildMaterial("G4_O", false);
	m_materials["Br"] = nistManager->FindOrBuildMaterial("G4_Br", false);

	m_materials["PCB"] = new G4Material("FR4", 1.700 * g / cm3, 5);
	m_materials["PCB"]->AddMaterial(m_materials["Si"], 0.18077359);
	m_materials["PCB"]->AddMaterial(m_materials["O"], 0.4056325);
	m_materials["PCB"]->AddMaterial(m_materials["C"], 0.27804208);
	m_materials["PCB"]->AddMaterial(m_materials["H"], 0.068442752);
	m_materials["PCB"]->AddMaterial(m_materials["Br"], 0.067109079);
	m_dEdx["PCB"] = 0;

	m_materials["Brass"] = new G4Material("Brass", 8.53 * g / cm3, 2);
	m_materials["Brass"]->AddMaterial(m_materials["Cu"], 70 * perCent);
	m_materials["Brass"]->AddMaterial(m_materials["Zn"], 30 * perCent);
	m_dEdx["Brass"] = 0.7 * m_dEdx["Cu"] + 0.3 * m_dEdx["Zn"];
	m_materials["Steel"] = new G4Material("Steel", 7.87 * g / cm3, 3);
	m_materials["Steel"]->AddMaterial(m_materials["Fe"], 0.9843);
	m_materials["Steel"]->AddMaterial(m_materials["Mn"], 0.014);
	m_materials["Steel"]->AddMaterial(m_materials["C"], 0.0017);
	m_dEdx["Steel"] = 0.9843 * m_dEdx["Fe"] + 0.014 * m_dEdx["Mn"]
			+ 0.0017 * m_dEdx["C"];
	m_materials["SSteel"] = new G4Material("SSteel", 8.02 * g / cm3, 4);
	m_materials["SSteel"]->AddMaterial(m_materials["Fe"], 0.70);
	m_materials["SSteel"]->AddMaterial(m_materials["Mn"], 0.01);
	m_materials["SSteel"]->AddMaterial(m_materials["Cr"], 0.19);
	m_materials["SSteel"]->AddMaterial(m_materials["Ni"], 0.10);
	m_dEdx["SSteel"] = 0.7 * m_dEdx["Fe"] + 0.01 * m_dEdx["Mn"]
			+ 0.19 * m_dEdx["Cr"] + 0.1 * m_dEdx["Ni"];

	m_materials["Scintillator"] = new G4Material("Scintillator",
			1.032 * g / cm3, 2);
	m_materials["Scintillator"]->AddMaterial(m_materials["C"],
			91.512109 * perCent);
	m_materials["Scintillator"]->AddMaterial(m_materials["H"],
			8.4878906 * perCent);
	m_dEdx["Scintillator"] = m_dEdx["C"];

	G4cout << m_materials["Scintillator"] << G4endl;
	m_materials["Polystyrole"] = new G4Material("Polystyrole", 1.065 * g / cm3,
			2);
	m_materials["Polystyrole"]->AddMaterial(m_materials["H"], 50 * perCent);
	m_materials["Polystyrole"]->AddMaterial(m_materials["C"], 50 * perCent);
	m_dEdx["Polystyrole"] = 0.5 * m_dEdx["C"];

	m_materials["PVC"] = new G4Material("PVC", 1.350 * g / cm3, 3);
	m_materials["PVC"]->AddMaterial(m_materials["H"], 50 * perCent);
	m_materials["PVC"]->AddMaterial(m_materials["C"], 33.33 * perCent);
	m_materials["PVC"]->AddMaterial(m_materials["Cl"], 16.67 * perCent);
	m_dEdx["PVC"] = 0.33 * m_dEdx["C"];

	m_materials["CFMix"] = new G4Material("CFMix", 0.845 * g / cm3, 3);
	m_materials["CFMix"]->AddMaterial(m_materials["C"], 0.84491305);
	m_materials["CFMix"]->AddMaterial(m_materials["H"], 0.042542086);
	m_materials["CFMix"]->AddMaterial(m_materials["O"], 0.11254487);
	m_dEdx["CFMix"] = 0;

	m_materials["Foam"] = new G4Material("Foam", 0.0999 * g / cm3, 2);
	m_materials["Foam"]->AddMaterial(m_materials["C"], 0.856);
	m_materials["Foam"]->AddMaterial(m_materials["H"], 0.144);
	m_dEdx["Foam"] = 1.749 * 0.856 * 0.0999 / 10.;

	m_materials["WCu"] = new G4Material("WCu", 14.979 * g / cm3, 2);
	m_materials["WCu"]->AddMaterial(m_materials["W"], 75 * perCent);
	m_materials["WCu"]->AddMaterial(m_materials["Cu"], 25 * perCent);
	m_dEdx["WCu"] = 0.75 * m_dEdx["W"] + 0.25 * m_dEdx["Cu"];

	m_materials["NeutMod"] = new G4Material("NeutMod", 0.950 * g / cm3, 2);
	m_materials["NeutMod"]->AddMaterial(m_materials["C"], 0.85628);
	m_materials["NeutMod"]->AddMaterial(m_materials["H"], 0.14372);
	m_dEdx["NeutMod"] = 1.749 * 0.86 * 0.950 / 10.;

}

//
void DetectorConstruction::UpdateCalorSize() {

	m_CalorSizeZ = 0;
	G4double maxRadLen = 0;
	for (size_t i = 0; i < m_caloStruct.size(); i++) {
		m_CalorSizeZ = m_CalorSizeZ + m_caloStruct[i].Total_thick;
		if (i != 0)
			maxRadLen = maxRadLen + m_caloStruct[i].Total_thick;
	}
	maxRadLen = maxRadLen + 200;
	m_nSectors = 1;
	if (model_ == DetectorConstruction::m_SIMPLE_20) {
		m_CalorSizeXY = 200;
		m_sectorWidth = m_CalorSizeXY;
	} else if (model_ == DetectorConstruction::m_SIMPLE_50) {
		m_CalorSizeXY = 500;
		m_sectorWidth = m_CalorSizeXY;
	} else if (model_ == DetectorConstruction::m_SIMPLE_100) {
		m_CalorSizeXY = 1000;
		m_sectorWidth = m_CalorSizeXY;
	}

	else if (model_ == DetectorConstruction::m_FULLSECTION) {
		m_maxTheta = pi / 6.0;
		cout << "Constructing the model = " << model_ << "Calo." << endl;

		cout << "The m_CalorSizeZ " << m_CalorSizeZ << endl;
		cout << "The maxRadLen " << maxRadLen << endl;

		m_maxRadius = (maxRadLen ) * tan(m_maxTheta);
		cout << "The maximum radius is " << m_maxRadius << endl;
		m_CalorSizeXY = m_maxRadius * 2; //use full length for making hexagon map
		m_sectorWidth = m_CalorSizeXY;
	}

	for (size_t i = 0; i < m_caloStruct.size(); i++)
		m_caloStruct[i].setNumberOfSectors(m_nSectors);

	m_WorldSizeZ = m_CalorSizeZ * 1.1;
	if (m_nSectors > 1)
		m_WorldSizeXY = (m_CalorSizeXY + 2 * m_sectorWidth) * 1.1;
	else
		m_WorldSizeXY = m_CalorSizeXY * 1.1;

	G4cout << "[DetectorConstruction][UpdateCalorSize] Z x XY = "
			<< m_CalorSizeZ << " x " << m_CalorSizeXY << " mm "
			<< ", nsectors = " << m_nSectors << G4endl;

}

//
G4VPhysicalVolume* DetectorConstruction::Construct() {

	//clean old geometry
	G4GeometryManager::GetInstance()->OpenGeometry();
	G4PhysicalVolumeStore::GetInstance()->Clean();
	G4LogicalVolumeStore::GetInstance()->Clean();
	G4SolidStore::GetInstance()->Clean();

	//world
	G4double expHall_z = 6 * m;
	G4double expHall_x = 3 * m;
	G4double expHall_y = 3 * m;

	G4Box* experimentalHall_box = new G4Box("expHall_box", expHall_x, expHall_y,
			expHall_z);

	G4LogicalVolume* experimentalHall_log = new G4LogicalVolume(
			experimentalHall_box, m_materials["Air"], "expHall_log");
	G4VPhysicalVolume* experimentalHall_phys = new G4PVPlacement(0, // no rotation
			G4ThreeVector(0., 0., 0.), // translation position
			experimentalHall_log,    // its logical volume
			"expHall",               // its name
			0,                       // its mother volume
			false,                   // no boolean operations
			0);                      // its copy number

	//detector's World
	G4double pos_x = 0.;
	G4double pos_y = 0.;
	G4double pos_z = 0.;

	m_solidWorld = new G4Box("Wbox", m_WorldSizeXY / 2, m_WorldSizeXY / 2,
			m_WorldSizeZ / 2);

	cout << "m_WorldSizeXY = " << m_WorldSizeXY << ", m_WorldSizeZ = "
			<< m_WorldSizeZ << endl;

	m_logicWorld = new G4LogicalVolume(m_solidWorld, m_materials["Air"],
			"Wlog");
	m_physWorld = new G4PVPlacement(0, G4ThreeVector(pos_x, pos_y, pos_z),
			m_logicWorld, "Wphys", experimentalHall_log, false, 0);

	for (unsigned iS(0); iS < m_nSectors; ++iS) {
		G4double minL = m_sectorWidth * iS;
		buildSectorStack(iS, minL, m_sectorWidth);
		if (m_nSectors > 1)
			fillInterSectorSpace(iS, minL + m_sectorWidth, m_interSectorWidth);
	}
	// Visualization attributes
	//
	m_logicWorld->SetVisAttributes(G4VisAttributes::Invisible);

	//return m_physWorld;
	return experimentalHall_phys;
}

void DetectorConstruction::buildSectorStack(const unsigned sectorNum,
		const G4double & minL, const G4double & width) {

	//build the stack
	G4double zOffset(-m_CalorSizeZ / 2), zOverburden(0.);
	char nameBuf[10];
	G4VSolid *solid;

	G4double totalLengthX0 = 0;
	G4double totalLengthL0 = 0;

	for (size_t i = 0; i < m_caloStruct.size(); i++) {
		const unsigned nEle = m_caloStruct[i].n_elements;
		//index for counting Si sensitive layers
		unsigned idx = 0;

		for (unsigned ie(0); ie < nEle; ++ie) {
			std::string eleName = m_caloStruct[i].ele_name[ie];

			if (m_nSectors == 1)
				sprintf(nameBuf, "%s%d", eleName.c_str(), int(i + 1));

			else
				sprintf(nameBuf, "%s%d_%d", eleName.c_str(), int(sectorNum),
						int(i + 1));

			if (eleName == "Si") {
				if (m_nSectors == 1)
					sprintf(nameBuf, "Si%d_%d", int(i + 1), idx);

				else
					sprintf(nameBuf, "Si%d_%d_%d", int(sectorNum), int(i + 1),
							idx);
				idx++;
			}
			std::string baseName(nameBuf);
			G4double thick = m_caloStruct[i].sublayer_thick[ie];
			//

			if (thick > 0) {
#if 0
				cout << "solid = constructSolid("<<baseName
				<<",thick="<<thick
				<<",zOffset+zOverburden="<<zOffset+zOverburden
				<<",width="<<width<<");"<<endl;
#endif
				solid = constructSolid(baseName, thick, zOffset + zOverburden,
						minL, width, i);
				G4LogicalVolume *logi = new G4LogicalVolume(solid,
						m_materials[eleName], baseName + "log");
				m_caloStruct[i].sublayer_X0[ie] = m_materials[eleName]->GetRadlen();
				m_caloStruct[i].sublayer_dEdx[ie] = m_dEdx[eleName];
				m_caloStruct[i].sublayer_L0[ie] =
						m_materials[eleName]->GetNuclearInterLength();
				if (sectorNum == 0 || sectorNum == m_nSectors - 1) {
					G4cout << "************ " << eleName;
					if (m_nSectors > 1)
						G4cout << " sector " << sectorNum;
					G4cout << " layer " << i << " dEdx="
							<< m_caloStruct[i].sublayer_dEdx[ie] << " X0="
							<< m_caloStruct[i].sublayer_X0[ie] << " L0="
							<< m_caloStruct[i].sublayer_L0[ie] << " zpos="
							<< zOffset + zOverburden << "mm w="
							<< m_caloStruct[i].sublayer_thick[ie] << "mm";

					totalLengthX0 += m_caloStruct[i].sublayer_thick[ie]
							/ m_caloStruct[i].sublayer_X0[ie];
					G4cout << " TotX0=" << totalLengthX0;
					totalLengthL0 += m_caloStruct[i].sublayer_thick[ie]
							/ m_caloStruct[i].sublayer_L0[ie];
					G4cout << " TotLambda=" << totalLengthL0 << G4endl;
				}

				if (m_caloStruct[i].isSensitiveElement(ie))
					m_logicSi.push_back(logi);

				G4double xpvpos = -m_CalorSizeXY / 2. + minL + width / 2;
#if 0
				cout << "m_caloStruct[i].ele_vol[nEle*sectorNum+ie]=new G4PVPlacement(0, G4ThreeVector(xpvpos="<<xpvpos
				<< ",0.,zOffset+zOverburden+thick/2="<<zOffset+zOverburden+thick/2
				<< "), logi,"
				<< baseName+"phys, m_logicWorld, false, 0);" << endl;
#endif
				m_caloStruct[i].sublayer_vol[nEle * sectorNum + ie] =
						new G4PVPlacement(0,
								G4ThreeVector(xpvpos, 0.,
										zOffset + zOverburden + thick / 2),
								logi, baseName + "phys", m_logicWorld, false,
								0);

				G4VisAttributes *simpleBoxVisAtt = new G4VisAttributes(
						m_caloStruct[i].g4Colour(ie));
				simpleBoxVisAtt->SetVisibility(true);
				logi->SetVisAttributes(simpleBoxVisAtt);
				zOverburden = zOverburden + thick;
				//for sensitive volumes
				//add region to be able to set specific cuts for it
				//just for Si
				if (eleName == "Si") {
					unsigned nlogicsi = m_logicSi.size();
					G4Region* aRegion = new G4Region(baseName + "Reg");
					m_logicSi[nlogicsi - 1]->SetRegion(aRegion);
					aRegion->AddRootLogicalVolume(m_logicSi[nlogicsi - 1]);
				}
			}

		}	  //loop on elements
	}	  //loop on layers

}	  //buildstack

void DetectorConstruction::fillInterSectorSpace(const unsigned sectorNum,
		const G4double & minL, const G4double & width) {
	//build the stack
	G4double zOffset(-m_CalorSizeZ / 2), zOverburden(0.);
	char nameBuf[10];
	G4VSolid *solid;

	for (size_t i = 0; i < m_caloStruct.size(); i++) {
		const unsigned nEle = m_caloStruct[i].n_elements;
		for (unsigned ie(0); ie < nEle; ++ie) {

			std::string eleName = m_caloStruct[i].ele_name[ie];
			G4double thick = m_caloStruct[i].sublayer_thick[ie];

			eleName = "CFMix";
			sprintf(nameBuf, "%s%d_%d", eleName.c_str(), int(sectorNum),
					int(i + 1));
			std::string baseName(nameBuf);
			if (thick > 0) {
				solid = constructSolid(baseName, thick, zOffset + zOverburden,
						minL, width, i);
				G4LogicalVolume *logi = new G4LogicalVolume(solid,
						m_materials[eleName], baseName + "log");
				G4double xpvpos = -m_CalorSizeXY / 2. + minL + width / 2;
				G4PVPlacement *tmp = new G4PVPlacement(0,
						G4ThreeVector(xpvpos, 0.,
								zOffset + zOverburden + thick / 2), logi,
						baseName + "phys", m_logicWorld, false, 0);

				G4VisAttributes *simpleBoxVisAtt = new G4VisAttributes(
						G4Colour::Magenta);
				simpleBoxVisAtt->SetVisibility(true);
				simpleBoxVisAtt->SetForceSolid(true);
				logi->SetVisAttributes(simpleBoxVisAtt);
				zOverburden = zOverburden + thick;
			}
		}	  //loop on elements
	}	  //loop on layers
}	  //fill intersector space

//
void DetectorConstruction::SetMagField(G4double fieldValue) {

	if (fieldValue <= 0)
		return;

	//apply a global uniform magnetic field along Z axis
	G4FieldManager* fieldMgr =
			G4TransportationManager::GetTransportationManager()->GetFieldManager();
	if (m_magField)
		delete m_magField;                //delete the existing magn field
	m_magField = new G4UniformMagField(G4ThreeVector(0., 0., fieldValue));
	fieldMgr->SetDetectorField(m_magField);
	fieldMgr->CreateChordFinder(m_magField);
	fieldMgr->SetDetectorField(m_magField);
}

void DetectorConstruction::SetDetModel(G4int model) {
	if (model <= 0)
		return;
	std::cout << " -- Setting detector model to " << model << std::endl;
	model_ = model;
}

G4VSolid *DetectorConstruction::constructSolid(std::string baseName,
		G4double thick, G4double zpos, const G4double & minL,
		const G4double & width, size_t which_ele) {
	G4VSolid *solid;

	if (which_ele == 0) {
		solid = new G4Box(baseName + "box", width / 2, m_CalorSizeXY / 2,
				thick / 2);
		//set the offset!
		if (baseName == "W1"){
			m_z0pos = zpos;
		}
	} else {
		if (model_ == DetectorConstruction::m_FULLSECTION) {
			double layerR = tan(m_maxTheta) * (zpos - m_z0pos);
			std::cout << "The radius for this layer is " << layerR << std::endl;
			solid = new G4Tubs(baseName + "box", 0, layerR, thick / 2, minL,
					width);
		}

		else {
			solid = new G4Box(baseName + "box", width / 2, m_CalorSizeXY / 2,
					thick / 2);
		}
	}
	return solid;
}
