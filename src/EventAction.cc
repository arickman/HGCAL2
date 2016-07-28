#include "EventAction.hh"

#include "EventActionMessenger.hh"
#include "DetectorConstruction.hh"
#include "StackingAction.hh"
#include <iostream>

#include "HGCSSInfo.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>
#include "TVector3.h"
#include "TMath.h"
using namespace std;
//
EventAction::EventAction(G4bool doFast) {
	eventMessenger = new EventActionMessenger(this);
	printModulo = 100;
	doFast_ = doFast;
	outF_ = TFile::Open("PFcal.root", "RECREATE");
	//summedDep = 0;nSteps = 0;nMainSteps = 0;
	//depCut = 150;
	/*for (Int_t i = 0; i < 1000000;  i++){
		step[i] = i;
		stepMain[i] = i;
	}*/
	hadronicInts = 0;
	outF_->cd();
	double xysize =
			((DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction())->GetCalorSizeXY();
	initLayer = ((DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction())->initLayer();
	//stacker_  = (StackingAction*) G4RunManager::GetRunManager()->GetUserStackingAction();

	//save some info
	HGCSSInfo *info = new HGCSSInfo();
	info->calorSizeXY(xysize);
	info->cellSize(CELL_SIZE_X);
	info->model(
			((DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction())->getModel());
	info->version(
			((DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction())->getVersion());

	std::cout << " -- check Info: version = " << info->version() << " model = "
			<< info->model() << std::endl;
	outF_->WriteObjectAny(info, "HGCSSInfo", "Info");

	tree_ = new TTree("HGCSSTree", "HGC Standalone simulation tree");
	tree_->Branch("HGCSSEvent", "HGCSSEvent", &event_);
	tree_->Branch("HGCSSGenAction", "std::vector<HGCSSGenParticle>",
			&genvec_);
	tree_->Branch("HGCSSHadAction", "std::vector<HGCSSGenParticle>",
			&hadvec_);
	tree_->Branch("HGCSSIncAction", "std::vector<HGCSSGenParticle>",
			&incvec_);
	tree_->Branch("HGCSSEscapeAction", "std::vector<HGCSSGenParticle>",
			&escapevec_);
	tree_->Branch("HGCSSNovelAction", "std::vector<HGCSSGenParticle>",
			&novelVec_);
	//tree_->Branch("nSteps",&nSteps,"nSteps/I");
	//tree_->Branch("nMainSteps",&nMainSteps,"nMainSteps/I");

	//tree_->Branch("step",&step,"step[nSteps]/I");
	//tree_->Branch("stepDep",&stepDep,"stepDep[nSteps]/F");

	//tree_->Branch("stepMain",&stepMain,"stepMain[nMainSteps]/I");
	//tree_->Branch("mainKinEng",&mainKinEng,"mainKinEng[nMainSteps]/F");


	// }
}

//
EventAction::~EventAction() {
	outF_->cd();
	tree_->Write();
	outF_->Close();
	//fout_.close();
	delete eventMessenger;
}

//
void EventAction::BeginOfEventAction(const G4Event* evt) {
	evtNb_ = evt->GetEventID();
	if (evtNb_ % printModulo == 0) {
		G4cout << "\n---> Begin of event: " << evtNb_ << G4endl;
		//CLHEP::HepRandom::restoreEngineStatus("temp.rndm");
		//CLHEP::HepRandom::showEngineStatus();
	}

}

//
void EventAction::Detect(G4double eDepRaw, G4int trackID,G4double kinEng, G4VPhysicalVolume *volume) {
	std::pair<G4bool,G4bool> stopIter = std::make_pair(false,false);

	for (size_t i = initLayer; i < detector_->size(); i++){
		if (stopIter.first) break;
		stopIter = (*detector_)[i].add( eDepRaw, volume);
	}
	/*if (stopIter.second)
		summedDep += eDepRaw;
	stepDep[nSteps] = summedDep;
	if (trackID == 1 && kinEng > 500){
		stepMain [nMainSteps] = nSteps ;
		mainKinEng [nMainSteps] = kinEng ;
		nMainSteps = nMainSteps + 1;
	}
	if (summedDep > depCut) {
		nSteps = 0;
		nMainSteps = 0;
		G4RunManager::GetRunManager()->AbortEvent();
	}
	nSteps = nSteps + 1;
*/
}

//


void EventAction::EndOfEventAction(const G4Event* g4evt) {
	//return;

	event_.eventNumber(evtNb_);
	event_.steelThick(((DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction())->GetSteelThick());
	double totalSens = 0;
	double wgtTotalSens = 0;

	//if (summedDep < depCut){
	G4String fileN = "temp.rndm";
	CLHEP::HepRandom::saveEngineStatus(fileN);
	CLHEP::HepRandom::showEngineStatus();
	std::ifstream input(fileN);
	std::string currentLine;
	//Double_t stat_x = 0,stat_y = 0,seed_x = 0,seed_y = 0;
	Double_t stat_x;
	Double_t stat_y;
	Double_t seed_x;
	Double_t seed_y;
	for(int count = 0; count < 5; count++ ){
		getline( input, currentLine );
		if (count == 1){
			stat_x = std::atoi(currentLine.c_str());
			cout << "the current line is: "<<currentLine.c_str()<<endl;
			cout<< stat_x<<endl;
		}
		if (count == 2){
			stat_y = std::atoi(currentLine.c_str());
			cout<< stat_y<<endl;
		}

		if (count == 3){
			seed_x = std::atoi(currentLine.c_str());
			cout<< seed_x<<endl;
		}
		if (count == 4){
			seed_y = std::atoi(currentLine.c_str());
			cout<< "the seed I am putting in is:" <<seed_y<<endl;
		}
	//}
		TVector3 status(stat_x,stat_y,0);
		TVector3 seeds(seed_x,seed_y,0);

		event_.seeds(seeds);
		event_.status(status);
		//Changing initLayer because initial layers contain tracking sections.
		for (size_t i = ((DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction())->initLayer()
				; i < detector_->size(); i++) {
			Double_t weight = (i < 8) ? .8 : 1. ;
			totalSens += (*detector_)[i].getTotalSensE();
			wgtTotalSens += weight*(*detector_)[i].getTotalSensE();

			(*detector_)[i].resetCounters();
			} //loop on sensitive layers
		//G4cout << "This was a good event, the totalSens was " << totalSens << G4endl;

	}
	/*
	else{
		for (size_t i = ((DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction())->initLayer()
				; i < detector_->size(); i++) {
			Double_t weight = (i < 8) ? .8 : 1. ;
			totalSens += (*detector_)[i].getTotalSensE();
			wgtTotalSens += weight*(*detector_)[i].getTotalSensE();

			(*detector_)[i].resetCounters();
			} //loop on sensitive layers
		//G4cout << "This was a good event, the totalSens was " << totalSens << G4endl;
	}*/
	event_.dep(totalSens);
	event_.wgtDep(wgtTotalSens);
	SetWait(false);
	//G4cout << "The dep cut is " << depCut << " The totalSens is " << totalSens << " The summedDep is " << summedDep << G4endl;
	tree_->Fill();
	//summedDep = 0;
	//nSteps = 0;
	//nMainSteps = 0;
	hadronicInts = 0;
	//reset vectors
	genvec_.clear();
	hadvec_.clear();
	incvec_.clear();
	escapevec_.clear();
	targetPartEngs.clear();
	novelVec_.clear();
	novelPartEngs.clear();
}
