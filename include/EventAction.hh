#ifndef EventAction_h
#define EventAction_h 1

#include "SamplingSection.hh"

#include "G4ThreeVector.hh"
#include "G4UserEventAction.hh"
#include "globals.hh"

#include "TFile.h"
#include "TTree.h"
#include "SamplingSection.hh"
#include "G4SiHit.hh"
#include "HGCSSEvent.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSGenParticle.hh"
#include "HGCSSGeometryConversion.hh"

#include <vector>
#include <map>
#include "fstream"

class EventActionMessenger;
class StackingAction;

class EventAction: public G4UserEventAction {
public:
	EventAction(G4bool doFast);
	virtual ~EventAction();
	void BeginOfEventAction(const G4Event*);
	void EndOfEventAction(const G4Event*);

	void Detect(G4double eRawDep,G4int trackID,G4double kinEng,G4VPhysicalVolume *volume);

	void SetPrintModulo(G4int val) {
		printModulo = val;
	}
	;
	void Add(std::vector<SamplingSection> *newDetector) {
		detector_ = newDetector;
	}
	inline void SetWait(G4bool wait) {
		wait_ = wait;
	}
	;
	inline G4bool GetWait() {
		return wait_;
	}
	;
	inline G4bool GetSpeed() {
		return doFast_;
	}
	;
	HGCSSGenParticleVec genvec_,hadvec_,incvec_,escapevec_,novelVec_;
	G4int hadronicInts;
	std::vector<double> targetPartEngs;
	G4int evtNb_;
	std::vector<int> novelPartEngs;
private:
	std::vector<SamplingSection> *detector_;
	G4int printModulo,initLayer;
	G4bool wait_;
	G4bool doFast_;
	//int nSteps,nMainSteps;
	//G4double depCut;
	//G4double summedDep;
	HGCSSGeometryConversion* geomConv_;
	//Int_t step[1000000],stepMain[1000000];
	//Float_t stepDep[1000000],mainKinEng[1000000];
	//StackingAction* stacker_;
	TFile *outF_;
	TTree *tree_;
	HGCSSEvent event_;

	EventActionMessenger* eventMessenger;

	//std::ofstream fout_;

};

#endif

