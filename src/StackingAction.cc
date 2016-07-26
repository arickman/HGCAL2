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
/// \file electromagnetic/TestEm1/src/StackingAction.cc
/// \brief Implementation of the StackingAction class
//
// $Id: StackingAction.cc 81776 2014-06-05 08:41:01Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "StackingAction.hh"
#include "G4Track.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::StackingAction(std::string data)
 : G4UserStackingAction()
{
	//wait_ = false;
	eventAction_ =
			(EventAction*) G4RunManager::GetRunManager()->GetUserEventAction();
	stepAction_ =
			(SteppingAction*) G4RunManager::GetRunManager()->GetUserSteppingAction();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::~StackingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
StackingAction::ClassifyNewTrack(const G4Track* lTrack)
{
	G4double kinEng = lTrack->GetKineticEnergy();
	G4int pdgID = lTrack->GetDefinition()->GetPDGEncoding();
	/*
	G4int trackID = lTrack->GetTrackID();
	unsigned int hadronTrackLoc = std::find(eventAction_->targetPartEngs.begin(),
			eventAction_->targetPartEngs.end(), trackID)
			- eventAction_->targetPartEngs.begin();
	G4double kineng = lTrack->GetKineticEnergy();
	G4int pdgId = lTrack->GetDefinition()->GetPDGEncoding();
	 */
	if (kinEng>10 && (abs(pdgID) != 11) && (abs(pdgID) != 22 )){
			HGCSSGenParticle genPart;
			genPart.vertexKE(lTrack->GetKineticEnergy());
			const G4ThreeVector &p = lTrack->GetMomentumDirection();
			const G4ThreeVector &pos = lTrack->GetPosition();
			TVector3 momVec(p[0], p[1], p[2]);
			genPart.vertexMom(momVec);
			TVector3 posVec(pos[0], pos[1], pos[2] );
			genPart.vertexPos(posVec);
			genPart.mass(lTrack->GetDefinition()->GetPDGMass());
			genPart.pdgid(lTrack->GetDefinition()->GetPDGEncoding());
			genPart.parentKE(stepAction_->stepKE);
			genPart.parentPdgId(stepAction_->stepPDGID);

			eventAction_->novelVec_.push_back(genPart);
			eventAction_->novelPartEngs.push_back(lTrack->GetTrackID());

		}

	if ( ((abs(pdgID) == 11) ||  (abs(pdgID) == 22)) && kinEng < 500) {
		if (!eventAction_->GetWait()){
			if (eventAction_->GetSpeed()){
				return fKill;
			}
			else{
				return fUrgent;
			}
		}
		else{
			return fWaiting;
		}
	}
	else{
		  return fUrgent;
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
