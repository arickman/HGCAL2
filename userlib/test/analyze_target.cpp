#include<string>
#include<set>
#include<iostream>
#include<fstream>
#include<sstream>
#include <algorithm>
#include <boost/algorithm/string.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TNtuple.h"
#include "TH2D.h"
#include "TH2Poly.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "Math/Vector4D.h"

#include "fastjet/ClusterSequence.hh"

#include "HGCSSEvent.hh"
#include "HGCSSInfo.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSGenParticle.hh"
#include "HGCSSRecoJet.hh"
#include "HGCSSCalibration.hh"
#include "HGCSSDigitisation.hh"
#include "HGCSSDetector.hh"
#include "HGCSSGeometryConversion.hh"

#ifdef __MAKECINT__
#pragma link C++ class vector<float>+;
#endif
Double_t weights[26] = { 1.0, 1.7453056191650522, 1.0, 1.7453056191650522, 1.0,
		1.7453056191650522, 1.0, 1.7453056191650522, 1.3114320944160649,
		2.1620926545710764, 1.3114320944160649, 2.1620926545710764,
		1.3114320944160649, 2.1620926545710764, 1.3114320944160649,
		2.1620926545710764, 1.3114320944160649, 2.1620926545710764,
		1.8564382596441786, 2.8567377135811172, 1.8564382596441786,
		2.8567377135811172, 1.8564382596441786, 2.8567377135811172,
		1.8564382596441786, 2.8567377135811172 };

int main() {
	freopen("test_log_2.txt", "w", stdout);

	TFile *infile = TFile::Open(
			"/afs/cern.ch/work/m/mullin/signal/git_hexaV02-01-01/version_1/model_2/signal_1.000/run_1/HGcal__version1_model2_run1.root");
//			"/afs/cern.ch/work/o/ocolegro/signal_events.root");
	TTree *tree = (TTree*) infile->Get("RecoTree");

	std::vector<HGCSSRecoHit> * rechits = 0;
	std::vector<HGCSSGenParticle> * gen_partons = 0;

	tree->SetBranchAddress("HGCSSRecoHitVec", &rechits);
	tree->SetBranchAddress("HGCSSGenParticleVec", &gen_partons);

	unsigned nEvts = tree->GetEntries();
	//std::cout << "There are Nevents = " << nEvts  << std::endl;

	TFile hfile("Event_signal.root", "RECREATE");
	TTree t1("t1", "a simple Tree with simple variables");

	Int_t N_e, N_p, N_nf, N_ni, N_hits, Evt;
	t1.Branch("Event", &Evt, "Evt/I");
	t1.Branch("N_e", &N_e, "N_e/I");
	t1.Branch("N_p", &N_p, "N_p/I");
	t1.Branch("N_nf", &N_nf, "N_nf/I");
	t1.Branch("N_ni", &N_ni, "N_ni/I");
	t1.Branch("N_hits", &N_hits, "N_hits/I");

	Float_t Eng_dep, N_pvec_init, N_pvec_final, N_pscalar_final, N_pscalar_init;
	t1.Branch("Eng_dep", &Eng_dep, "Eng_dep/F");
	t1.Branch("N_pvec_init", &N_pvec_init, "N_pvec_init/F");
	t1.Branch("N_pvec_final", &N_pvec_final, "N_pvec_final/F");
	t1.Branch("N_pscalar_final", &N_pscalar_final, "N_pscalar_final/F");
	t1.Branch("N_pscalar_init", &N_pscalar_init, "N_pscalar_init/F");

	Float_t E_theta[50], E_pt[50], E_pz[50], P_theta[50], P_pt[50], P_pz[50],
			N_theta_final[500], N_pt_final[500], N_pz_final[500],
			N_theta_initial[500], N_pt_initial[500], N_pz_initial[500];

	t1.Branch("E_theta", E_theta, "E_theta[N_e]/F");
	t1.Branch("E_pt", E_pt, "E_pt[N_e]/F");
	t1.Branch("E_pz", E_pz, "E_pz[N_e]/F");
	t1.Branch("P_theta", P_theta, "P_theta[N_p]/F");
	t1.Branch("P_pt", P_pt, "P_pt[N_p]/F");
	t1.Branch("P_pz", P_pz, "P_pz[N_p]/F");
	t1.Branch("N_theta_final", N_theta_final, "N_theta_final[N_nf]/F");
	t1.Branch("N_pt_final", N_pt_final, "N_pt_final[N_nf]/F");
	t1.Branch("N_pz_final", N_pz_final, "N_pz_final[N_nf]/F");
	t1.Branch("N_theta_initial", N_theta_initial, "N_theta_initial[N_ni]/F");
	t1.Branch("N_pt_initial", N_pt_initial, "N_pt_initial[N_ni]/F");
	t1.Branch("N_pz_initial", N_pz_initial, "N_pz_initial[N_ni]/F");

	for (unsigned ievt(0); ievt < nEvts; ++ievt) { //loop on entries
		tree->GetEntry(ievt);
		Double_t eng = 0;
		Int_t len_hits = rechits->size();

		std::vector<Int_t> trackids;
		std::vector<HGCSSGenParticle> final_state;
		std::vector<HGCSSGenParticle> init_state;

		//Fill up the good partons :)
		for (Int_t j = 0; j < gen_partons->size(); j++) {
			HGCSSGenParticle& parton = (*gen_partons)[j];
			Int_t track_ID = parton.trackID();
			Int_t id_ = std::find(trackids.begin(), trackids.end(), track_ID)
					- trackids.begin();

			if ((id_ == trackids.size()) || (track_ID == 1)) {
				trackids.push_back(track_ID);
				final_state.push_back(parton);
				init_state.push_back(parton);

			} else {
				Double_t current_z = final_state.at(id_).z();
				if (current_z < parton.z())
					//Overwrite that parton
					final_state[id_] = parton;
				if (current_z > parton.z())
					//Overwrite that parton
					init_state[id_] = parton;
			}
		}

		N_e = 0;
		N_p = 0;
		N_nf = 0;
		N_ni = 0, N_pscalar_final = 0;
		N_pscalar_init = 0;
		TVector3 N_final = TVector3(0, 0, 0);
		//Calculate Parton Quantities
		for (Int_t j = 0; j < final_state.size(); j++) {
			if ((j == 0) && final_state.at(j).pz() < 4000)
				break;
			else if (j == 0)
				continue;
			HGCSSGenParticle parton = final_state.at(j);
			if (parton.pdgid() == 11) {
				E_theta[N_e] = parton.theta();
				E_pt[N_e] = parton.pt();
				E_pz[N_e] = parton.pz();
				N_e += 1;
			}
			if (parton.pdgid() == 22) {
				P_theta[N_p] = parton.theta();
				P_pt[N_p] = parton.pt();
				P_pz[N_p] = parton.pz();
				N_p += 1;
			}
			if (parton.pdgid() == 2112) {
				N_theta_final[N_nf] = parton.theta();
				N_pt_final[N_nf] = parton.pt();
				N_pz_final[N_nf] = parton.pz();
				TVector3 temp_neut = TVector3(parton.px(), parton.py(),
						parton.pz());
				N_final += temp_neut;
				N_pscalar_final += temp_neut.Mag();
				N_nf += 1;
			}
		}
		N_pvec_final = N_final.Mag();

		TVector3 N_init = TVector3(0, 0, 0);
		for (Int_t j = 0; j < init_state.size(); j++) {
			if (init_state.at(j).pdgid() != 2112)
				continue;
			HGCSSGenParticle parton = init_state.at(j);
			N_theta_initial[N_ni] = parton.theta();
			N_pt_initial[N_ni] = parton.pt();
			N_pz_initial[N_ni] = parton.pz();
			TVector3 temp_neut = TVector3(parton.px(), parton.py(),
					parton.pz());

			N_init += temp_neut;
			N_pscalar_init += temp_neut.Mag();

			N_ni += 1;
		}
		if (N_ni > 75)
			tree->Show(ievt, 1000000);
		N_pvec_init = N_init.Mag();

		//Calculate the energy deposited
		for (Int_t j = 0; j < len_hits; j++) {
			HGCSSRecoHit& hit = (*rechits)[j];
			eng += weights[hit.layer()] * hit.energy();
		}
		N_hits = len_hits;
		Eng_dep = eng;
		Evt = ievt;
		t1.Fill();
	}
	t1.Write();

	return 1;
}
