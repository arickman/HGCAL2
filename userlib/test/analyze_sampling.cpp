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
#include "HGCSSSamplingSection.hh"
#ifdef __MAKECINT__
#pragma link C++ class vector<float>+;
#endif

int main(int argc, char** argv) {
	std::cout << "Opening the file " << argv[1] << std::endl;
	TFile *infile = TFile::Open(argv[1]);
	TTree *tree = (TTree*) infile->Get("HGCSSTree");
	freopen("log.txt", "w", stdout);

	std::vector<HGCSSSamplingSection> * samplingVec = 0;
	tree->SetBranchAddress("HGCSSSamplingSectionVec", &samplingVec);

	std::vector<HGCSSSimHit> * hitVec = 0;
	tree->SetBranchAddress("HGCSSSimHitVec", &hitVec);

	std::vector<HGCSSGenParticle> * hadronVec = 0;
	tree->SetBranchAddress("HGCSSHadronVec", &hadronVec);

	std::vector<HGCSSGenParticle> * targetVec = 0;
	tree->SetBranchAddress("HGCSSTargetVec", &targetVec);
	Int_t firstLayer = 0;
	unsigned nEvts = tree->GetEntries();

	TFile hfile("analyzed_tuple.root", "RECREATE");
	TTree t1("sampling", "Sampling Study");

	Float_t summedTotalNonIon,summedTotal, summedSen, maxHadronKe,genKin,
			layerAvgHFlux,layerAvgNFlux,layerAvgEFlux,layerAvgGFlux,layerAvgMFlux,layerAvgEGFlux,
			layerAvgNCount,layerAvgHCount,layerAvgECount,layerAvgGCount,layerAvgMCount,
			layerAvgHFrac,layerAvgNFrac,layerAvgMFrac,layerAvgEFrac,layerAvgGFrac,layerAvgEGCount,layerAvgEGFrac,
			layerEShowerSizeAvg,layerHShowerSizeAvg,hardestEarlyHadron,
			layerHFlux[500], layerNFlux[500], layerMFlux[500],
			layerHShowerSize[500], layerEShowerSize[500], layerTotal[500],layerEGFrac[500],
			layerHFrac[500], layerNFrac[500], layerMFrac[500],layerGFrac[500],layerEFrac[500],
			layerSen[500],layerHCount[500],layerNCount[500],layerMCount[500],layerNonIon[500],
			layerEFlux[500],layerGFlux[500],layerECount[500],layerGCount[500],layerEGFlux[500],layerEGCount[500];

	Int_t layer[500], caloLen, maxHadronProdLayer,genCounter;
	t1.Branch("caloLen", &caloLen, "caloLen/I");

	t1.Branch("summedTotal", &summedTotal, "summedTotal/F");
	t1.Branch("summedSen", &summedSen, "summedSen/F");
	t1.Branch("summedTotalNonIon", &summedTotalNonIon, "summedTotalNonIon/F");

	t1.Branch("layerAvgHFlux", &layerAvgHFlux, "layerAvgHFlux/F");
	t1.Branch("layerAvgNFlux", &layerAvgNFlux, "layerAvgNFlux/F");
	t1.Branch("layerAvgMFlux", &layerAvgMFlux, "layerAvgMFlux/F");
	t1.Branch("layerAvgEFlux", &layerAvgEFlux, "layerAvgEFlux/F");
	t1.Branch("layerAvgGFlux", &layerAvgGFlux, "layerAvgGFlux/F");
	t1.Branch("layerAvgEGFlux", &layerAvgEGFlux, "layerAvgEGFlux/F");

	t1.Branch("layerAvgHCount", &layerAvgHCount, "layerAvgHCount/F");
	t1.Branch("layerAvgNCount", &layerAvgNCount, "layerAvgNCount/F");
	t1.Branch("layerAvgMCount", &layerAvgMCount, "layerAvgMCount/F");
	t1.Branch("layerAvgECount", &layerAvgECount, "layerAvgECount/F");
	t1.Branch("layerAvgGCount", &layerAvgGCount, "layerAvgGCount/F");
	t1.Branch("layerAvgEGCount", &layerAvgEGCount, "layerAvgEGCount/F");

	t1.Branch("layerAvgHFrac", &layerAvgHFrac, "layerAvgHFrac/F");
	t1.Branch("layerAvgNFrac", &layerAvgNFrac, "layerAvgNFrac/F");
	t1.Branch("layerAvgMFrac", &layerAvgMFrac, "layerAvgMFrac/F");
	t1.Branch("layerAvgEFrac", &layerAvgEFrac, "layerAvgEFrac/F");
	t1.Branch("layerAvgGFrac", &layerAvgGFrac, "layerAvgGFrac/F");

	t1.Branch("layerHShowerSizeAvg", &layerHShowerSizeAvg, "layerHShowerSizeAvg/F");
	t1.Branch("layerEShowerSizeAvg", &layerEShowerSizeAvg, "layerEShowerSizeAvg/F");


	t1.Branch("layerHFlux", &layerHFlux, "layerHFlux[caloLen]/F");
	t1.Branch("layerNFlux", &layerNFlux, "layerNFlux[caloLen]/F");
	t1.Branch("layerMFlux", &layerMFlux, "layerMFlux[caloLen]/F");
	t1.Branch("layerEFlux", &layerEFlux, "layerEFlux[caloLen]/F");
	t1.Branch("layerGFlux", &layerGFlux, "layerGFlux[caloLen]/F");
	t1.Branch("layerEGFlux", &layerEGFlux, "layerEGFlux[caloLen]/F");

	t1.Branch("layerHCount", &layerHCount, "layerHCount[caloLen]/F");
	t1.Branch("layerNCount", &layerNCount, "layerNCount[caloLen]/F");
	t1.Branch("layerMCount", &layerMCount, "layerMCount[caloLen]/F");
	t1.Branch("layerECount", &layerECount, "layerECount[caloLen]/F");
	t1.Branch("layerGCount", &layerGCount, "layerGCount[caloLen]/F");
	t1.Branch("layerEGCount", &layerEGCount, "layerEGCount[caloLen]/F");

	t1.Branch("layerHFrac", &layerHFrac, "layerHFrac[caloLen]/F");
	t1.Branch("layerNFrac", &layerNFrac, "layerNFrac[caloLen]/F");
	t1.Branch("layerMFrac", &layerMFrac, "layerMFrac[caloLen]/F");
	t1.Branch("layerEFrac", &layerEFrac, "layerEFrac[caloLen]/F");
	t1.Branch("layerGFrac", &layerGFrac, "layerGFrac[caloLen]/F");
	t1.Branch("layerEGFrac", &layerEGFrac, "layerEGFrac[caloLen]/F");

	t1.Branch("layerHShowerSize", &layerHShowerSize, "layerHShowerSize[caloLen]/F");
	t1.Branch("layerEShowerSize", &layerEShowerSize, "layerEShowerSize[caloLen]/F");

	t1.Branch("layerSen", &layerSen, "layerSen[caloLen]/F");
	t1.Branch("layerTotal", &layerTotal, "layerTotal[caloLen]/F");
	t1.Branch("layerNonIon", &layerNonIon, "layerNonIon[caloLen]/F");
	t1.Branch("layer", &layer, "layer[caloLen]/I");


	t1.Branch("maxHadronKe", &maxHadronKe, "maxHadronKe/F");
	t1.Branch("genKin", &genKin, "genKin/F");
	t1.Branch("genCounter", &genCounter, "genCounter/I");
	t1.Branch("maxHadronProdLayer", &maxHadronProdLayer, "maxHadronProdLayer/I");
	t1.Branch("hardestEarlyHadron", &hardestEarlyHadron, "hardestEarlyHadron/F");

	for (unsigned ievt(0); ievt < nEvts; ++ievt) { //loop on entries
		tree->GetEntry(ievt);

		summedSen = 0, summedTotal = 0, summedTotalNonIon = 0,caloLen = 0,
				layerAvgHFlux = 0, layerAvgNFlux = 0, layerAvgMFlux =0,layerAvgEFlux=0,layerAvgGFlux=0,layerAvgEGCount = 0,
				layerAvgHCount=0,layerAvgNCount=0,layerAvgMCount=0,layerAvgECount=0,layerAvgGCount=0,
				layerAvgHFrac=0,layerAvgNFrac=0,layerAvgMFrac=0,layerAvgEFrac=0,layerAvgGFrac=0,layerAvgEGFrac=0,
				maxHadronKe=0,genCounter = 0,genKin=0,
				layerHShowerSizeAvg=0,layerEShowerSizeAvg=0,hardestEarlyHadron = 0,layerAvgEGFlux=0;

		Float_t nSens = 3.0;
		Float_t nLayers = samplingVec->size();
		std::cout << "The number of layers is " << nLayers << std::endl;
		for (Int_t j = firstLayer; j < samplingVec->size(); j++) {
			HGCSSSamplingSection& sec = (*samplingVec)[j];
			summedSen += sec.sensDep();
			summedTotal += sec.totalDep();
			summedTotalNonIon += sec.totalNonIonDep();

			layerAvgHFlux += sec.hadKinFlux()/(nSens * nLayers);
			layerAvgNFlux += sec.neutronKinFlux()/(nSens * nLayers);
			layerAvgMFlux += sec.muKinFlux()/(nSens * nLayers);
			layerAvgEFlux += sec.eleKinFlux()/(nSens * nLayers);
			layerAvgGFlux += sec.gamKinFlux()/(nSens * nLayers);
			layerAvgEGFlux += (sec.eleKinFlux()+sec.gamKinFlux())/(nSens * nLayers);

			layerAvgHCount += sec.hadCount()/(nSens * nLayers);
			layerAvgNCount += sec.neutronCount()/(nSens * nLayers);
			layerAvgMCount += sec.muCount()/(nSens * nLayers);
			layerAvgECount += sec.eleCount()/(nSens * nLayers);
			layerAvgGCount += sec.gamCount()/(nSens * nLayers);
			layerAvgEGCount += ( sec.eleCount() +  sec.gamCount())/(nSens * nLayers);

			layerAvgHFrac += sec.hadDepFrac()/(nLayers);
			layerAvgNFrac += sec.neutronDepFrac()/(nLayers);
			layerAvgMFrac += sec.muDepFrac()/(nLayers);
			layerAvgEFrac += sec.eleDepFrac()/(nLayers);
			layerAvgGFrac += sec.gamDepFrac()/(nLayers);
			layerAvgEGFrac += (sec.gamDepFrac()+sec.gamDepFrac())/(nLayers);

			layerSen[j - firstLayer] = sec.sensDep()/nSens;
			layerTotal[j - firstLayer] = sec.totalDep()/nSens;
			layerNonIon[j - firstLayer] = sec.totalNonIonDep()/nSens;

			layerHFlux[j - firstLayer] = sec.hadKinFlux()/nSens;
			layerNFlux[j - firstLayer] = sec.neutronKinFlux()/nSens;
			layerMFlux[j - firstLayer] = sec.muKinFlux()/nSens;
			layerEFlux[j - firstLayer] = sec.eleKinFlux()/nSens;
			layerGFlux[j - firstLayer] = sec.gamKinFlux()/nSens;
			layerEGFlux[j - firstLayer] = (sec.eleKinFlux()+sec.gamKinFlux())/nSens;

			layerHCount[j - firstLayer] = sec.hadCount()/nSens;
			layerNCount[j - firstLayer] = sec.neutronCount()/nSens;
			layerMCount[j - firstLayer] = sec.muCount()/nSens;
			layerECount[j - firstLayer] = sec.eleCount()/nSens;
			layerGCount[j - firstLayer] = sec.gamCount()/nSens;
			layerEGCount[j - firstLayer] = (sec.eleCount()+sec.gamCount())/nSens;

			layerHFrac[j - firstLayer] = sec.hadDepFrac();
			layerNFrac[j - firstLayer] = sec.neutronDepFrac();
			layerMFrac[j - firstLayer] = sec.muDepFrac();
			layerEFrac[j - firstLayer] = sec.eleDepFrac();
			layerGFrac[j - firstLayer] = sec.gamDepFrac();
			layerEGFrac[j - firstLayer] = sec.gamDepFrac() + sec.eleDepFrac();

			layerHShowerSize[j - firstLayer] = sec.hadronShowerSize();
			layerEShowerSize[j - firstLayer] = sec.eleShowerSize();

			layerHShowerSizeAvg += sec.hadronShowerSize()/nLayers;
			layerEShowerSizeAvg += sec.eleShowerSize()/nLayers;

			layer[j - firstLayer] = j - firstLayer;
			caloLen = caloLen + 1;
		}
			for (Int_t j = 0; j < hadronVec->size(); j++) {
				HGCSSGenParticle& parton = (*hadronVec)[j];
				genCounter += 1;
				Float_t engK = parton.E() -parton.mass();
				genKin += engK;
				if (parton.layer() < 12){
					if (engK > hardestEarlyHadron)
						hardestEarlyHadron = engK;
				}
				if (engK> maxHadronKe){
					maxHadronKe = engK;
					maxHadronProdLayer  = parton.layer();
				}

		}

		t1.Fill();
	}
	t1.Write();
	return 1;
}
