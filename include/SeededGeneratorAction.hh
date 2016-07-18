//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id$
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef SeededGeneratorAction_h
#define SeededGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include <map>
#include "EventAction.hh"

class G4Event;
class DetectorConstruction;
class StackingAction;
class SeededGeneratorMessenger;
class G4VPrimaryGenerator;
class G4ParticleGun;
class HepMCG4AsciiReader;
class HepMCG4PythiaInterface;
class TFile;
class TTree;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class SeededGeneratorAction: public G4VUserPrimaryGeneratorAction {
public:
	SeededGeneratorAction(G4int mod = 0,
			std::string data = "");
	virtual ~SeededGeneratorAction();

	void GeneratePrimaries(G4Event*);
	void SetRndmFlag(G4String val) {
		rndmFlag = val;
	}
	void PipeData();

	void SetGenerator(G4VPrimaryGenerator* gen);
	void SetGenerator(G4String genname);

	G4VPrimaryGenerator* GetGenerator() const;
	G4String GetGeneratorName() const;

private:
	EventAction *eventAction_;
	int model_;

	TFile* file_;
	TTree * tree_;
	HGCSSGenParticleVec* inc_ ;
	HGCSSEvent* evt_;
	SeededGeneratorMessenger* gunMessenger; //messenger of this class

	std::string data_;
	G4ParticleGun* particleGun;
	HepMCG4AsciiReader* hepmcAscii;
	HepMCG4PythiaInterface* pythiaGen;

	G4VPrimaryGenerator* currentGenerator;
	G4String currentGeneratorName;
	std::map<G4String, G4VPrimaryGenerator*> gentypeMap;
	//StackingAction* stacker_;
	DetectorConstruction* Detector;     //pointer to the geometry
	G4String rndmFlag;     //flag for a rndm impact point
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// ====================================================================
// inline functions
// ====================================================================

inline void SeededGeneratorAction::SetGenerator(G4VPrimaryGenerator* gen) {
	currentGenerator = gen;
}

inline void SeededGeneratorAction::PipeData(){
	ofstream outputFile("temp.rndm");
	Int_t line1 = evt_-> status().x(), line2 = evt_-> status().y(),
			line3 = evt_-> seeds().x(), line4 = evt_-> seeds().y();
	outputFile << "Uvec\n" <<  line1 << "\n" << line2 << "\n"
			<< line3 << "\n" << line4;
}
inline void SeededGeneratorAction::SetGenerator(G4String genname) {
	std::map<G4String, G4VPrimaryGenerator*>::iterator pos = gentypeMap.find(
			genname);
	if (pos != gentypeMap.end()) {
		currentGenerator = pos->second;
		currentGeneratorName = genname;
	}
}

inline G4VPrimaryGenerator* SeededGeneratorAction::GetGenerator() const {
	return currentGenerator;
}

inline G4String SeededGeneratorAction::GetGeneratorName() const {
	return currentGeneratorName;
}
#endif

