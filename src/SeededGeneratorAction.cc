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

#include "SeededGeneratorAction.hh"

#include "DetectorConstruction.hh"
#include "SeededGeneratorMessenger.hh"
#include "StackingAction.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "globals.hh"

#include "HepMCG4AsciiReader.hh"
#include "HepMCG4PythiaInterface.hh"
#include "HGCSSGenParticle.hh"
#include "TTree.h"
#include "TFile.h"

#include <fstream>
#include <string>
#include <sstream>
#define PI 3.1415926535
#include <iostream>
using namespace std;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


SeededGeneratorAction::SeededGeneratorAction(G4int mod, std::string data) {
	model_ = mod;
	data_ = data;
	G4int n_particle = 1;
	evt_ = 0; inc_ = 0;
	data = "PFcal.root";
	file_ = TFile::Open(data.c_str());
	tree_  = (TTree*) file_->Get("HGCSSTree");
	tree_->SetBranchAddress("HGCSSIncAction",&inc_);
	tree_->SetBranchAddress("HGCSSEvent",&evt_);

	eventAction_ =
			(EventAction*) G4RunManager::GetRunManager()->GetUserEventAction();
	eventAction_->Add(
			((DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction())->getStructure());
	//stacker_ = (StackingAction*) G4RunManager::GetRunManager()->GetUserStackingAction();

	// default generator is particle gun.
	currentGenerator = particleGun = new G4ParticleGun(n_particle);
	currentGeneratorName = "particleGun";
	hepmcAscii = new HepMCG4AsciiReader();
#ifdef G4LIB_USE_PYTHIA
	pythiaGen= new HepMCG4PythiaInterface();
#else
	pythiaGen = 0;
#endif
	gentypeMap["particleGun"] = particleGun;
	gentypeMap["hepmcAscii"] = hepmcAscii;
	gentypeMap["pythia"] = pythiaGen;

	Detector =
			(DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction();

	//create a messenger for this class

	// default particle kinematic

	gunMessenger = new SeededGeneratorMessenger(this);

	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName;
	G4ParticleDefinition* particle = particleTable->FindParticle(particleName ="e-");

	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
	particleGun->SetParticleEnergy(4. * GeV);
	G4double position = -0.5 * (Detector->GetWorldSizeZ());
	particleGun->SetParticlePosition(G4ThreeVector(0. * cm, 0. * cm, position));

	G4cout << " -- Gun position set to: 0,0," << position << G4endl;

	rndmFlag = "off";

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SeededGeneratorAction::~SeededGeneratorAction() {
	delete particleGun;
	delete hepmcAscii;
	delete pythiaGen;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SeededGeneratorAction::GeneratePrimaries(G4Event* anEvent) {

	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName;

	G4ParticleDefinition* particle = particleTable->FindParticle(particleName =
			"e-");
	particleGun->SetParticleDefinition(particle);
	int currentEvt = anEvent->GetEventID();
	tree_->GetEntry(currentEvt);
	G4double et = 0.0;
	CLHEP::HepRandom::restoreEngineStatus ("temp.rndm");
	if (inc_->size() != 0){
		for(int i = 0; i < inc_->size(); i++){
			HGCSSGenParticle& incPart = (*inc_)[i];
			if (incPart.vertexKE() > 500){
				eventAction_->SetWait(true);
				et = 4.0;
			}
		}

	}
	else{
		PipeData();
		eventAction_->SetWait(false);
	}

	particleGun->SetParticleEnergy(et * GeV);
	particleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));

	G4double y0 = G4RandFlat::shoot(-10.,10);
	G4double x0 = G4RandFlat::shoot(-10.,10);
	G4double z0 = -0.5 * (Detector->GetWorldSizeZ());


	particleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
	HGCSSGenParticle genPart;
	genPart.vertexKE(et);
	TVector3 vec(x0,y0,z0);
	genPart.vertexPos(vec);
	int pdgid = particle->GetPDGEncoding();
	genPart.pdgid(pdgid);

	if (currentGenerator) {
		currentGenerator->GeneratePrimaryVertex(anEvent);
		eventAction_->genvec_.push_back(genPart);
	} else
		G4Exception("SeededGeneratorAction::GeneratePrimaries",
				"PrimaryGeneratorAction001", FatalException,
				"generator is not instanciated.");


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

