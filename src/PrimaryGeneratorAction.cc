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

#include "PrimaryGeneratorAction.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorMessenger.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "HepMCG4AsciiReader.hh"
#include "HepMCG4PythiaInterface.hh"
#include "HGCSSGenParticle.hh"

#include <fstream>
#include <string>
#include <sstream>
#define PI 3.1415926535

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<std::string> split(const std::string &s, char delim) {
	std::stringstream ss(s);
	std::string item;
	std::vector<std::string> tokens;
	while (std::getline(ss, item, delim)) {
		tokens.push_back(item);
	}
	return tokens;
}

PrimaryGeneratorAction::PrimaryGeneratorAction(G4int mod, bool signal,
		std::string data) {
	model_ = mod;
	signal_ = signal;
	data_ = data;
	G4int n_particle = 1;


	eventAction_ =
			(EventAction*) G4RunManager::GetRunManager()->GetUserEventAction();
	eventAction_->Add(
			((DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction())->getStructure());


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
	gunMessenger = new PrimaryGeneratorMessenger(this);

	// default particle kinematic

	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName;
	G4ParticleDefinition* particle = particleTable->FindParticle(particleName =
			"e-");
	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
	particleGun->SetParticleEnergy(4. * GeV);
	G4double position = -0.5 * (Detector->GetWorldSizeZ());
	particleGun->SetParticlePosition(G4ThreeVector(0. * cm, 0. * cm, position));

	G4cout << " -- Gun position set to: 0,0," << position << G4endl;

	rndmFlag = "off";

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
	delete particleGun;
	delete hepmcAscii;
	delete pythiaGen;
	delete gunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName;

	G4ParticleDefinition* particle = particleTable->FindParticle(particleName =
			"e-");
	particleGun->SetParticleDefinition(particle);
	G4double et = 4.0;
	//double engs3[16] = {.25,.5,.75,1.,1.25,1.5,1.75,2,2.25,2.5,2.75,3.,3.25,3.5,3.75,4.0};
	//unsigned iSecret = rand() % 16 ;
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
		G4Exception("PrimaryGeneratorAction::GeneratePrimaries",
				"PrimaryGeneratorAction001", FatalException,
				"generator is not instanciated.");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

