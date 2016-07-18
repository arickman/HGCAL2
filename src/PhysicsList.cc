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

#include "PhysicsList.hh"
#include "DetectorConstruction.hh"
#include "G4RunManager.hh"
#include "G4ProcessManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//PhysicsList::PhysicsList():  QGSP_FTFP_BERT() //G4VUserPhysicsList()
PhysicsList::PhysicsList() :
FTFP_BERT() //G4VUserPhysicsList()
{
	defaultCutValue = 0.03 * mm;
	SetVerboseLevel(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList() {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCuts() {
	if (verboseLevel > 0) {
		G4cout << "PhysicsList::SetCuts:";
		G4cout << "CutLength : " << G4BestUnit(defaultCutValue, "Length")
				<< G4endl;
			}

			// set cut values for gamma at first and for e- second and next for e+,
			// because some processes for e+/e- need cut values for gamma
			//

			SetCutValue(.7*mm, "gamma");
			SetCutValue(.7*mm, "e-");
			SetCutValue(.7*mm, "e+");
			SetCutValue(.7*mm, "proton");

			//set smaller cut for Si
			const std::vector<G4LogicalVolume*> & logSi = ((DetectorConstruction*)G4RunManager::GetRunManager()->GetUserDetectorConstruction())->getSiLogVol();

			for(size_t i=0; i<logSi.size(); i++) {
				G4Region* reg = logSi[i]->GetRegion();
				G4ProductionCuts* cuts = new G4ProductionCuts;
				cuts->SetProductionCut(defaultCutValue);
				reg->SetProductionCuts(cuts);
			}

			if (verboseLevel>0) DumpCutValuesTable();
		}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

