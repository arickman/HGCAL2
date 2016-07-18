#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include <boost/algorithm/string.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TRandom3.h"

#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSParameters.hh"

int main(int argc, char** argv){//main  

  if (argc < 3) {
    std::cout << " Usage: " 
	      << argv[0] << " <nEvts to process (0=all)>"
	      << " <full path to input file>"
	      << " <file name (PFcal.root, or DigiPFcal.root)>" 
	      << std::endl;
    return 1;
  }

  //////////////////////////////////////////////////////////
  //// Hardcoded config ////////////////////////////////////
  //////////////////////////////////////////////////////////

  const unsigned nLayers = 64;//64;//Calice 54;//Scint 9;//HCAL 33;//All 64
  const unsigned nEcalLayers = 31;//31;
  const unsigned nHcalSiLayers = 24;//concept 24;//calice 47

  //double minX=-250,maxX=250;
  double minX=-1700,maxX=1700;
  //double minY=150,maxY=800;
  //double minY=1100,maxY=1500;
  double minY=-1700,maxY=1700;
  double minZ=3170,maxZ=5070;
  
  unsigned nEvtsOut = 10;

  //////////////////////////////////////////////////////////
  //// End Hardcoded config ////////////////////////////////////
  //////////////////////////////////////////////////////////

  const unsigned pNevts = atoi(argv[1]);
  std::string filePath = argv[2];
  std::string fileName = argv[3];//"PFcal.root";
  std::cout << " -- Input parameters: " << std::endl
	    << " -- Input file path: " << filePath << std::endl
	    << " -- Input file name: " << fileName << std::endl
	    << " -- Processing ";
  if (pNevts == 0) std::cout << "all events." << std::endl;
  else std::cout << pNevts << " events." << std::endl;
  
  TRandom3 lRndm(0);

  std::cout << " ----------------------------------- " << std::endl
	    << " ----------------------------------- " << std::endl
	    << " ----- Hardcoded configuration is : " << std::endl
	    << " ----------------------------------- " << std::endl
	    << " ----------------------------------- " << std::endl
	    << " -- N layers = " << nLayers << std::endl
	    << " -- N ECAL layers = " << nEcalLayers << std::endl
	    << " ----------------------------------- " << std::endl
	    << " ----------------------------------- " << std::endl
	    << " -----------------------------------" << std::endl;

  bool isG4Tree = true;
  
  std::ostringstream input;
  input << filePath ;
  input << "/" << fileName;
  TFile *inputFile = TFile::Open(input.str().c_str());

  if (!inputFile) {
    std::cout << " -- Error, input file " << input.str() << " cannot be opened. Exiting..." << std::endl;
    return 1;
  }

  TTree *lTree = (TTree*)inputFile->Get("HGCSSTree");
  if (!lTree){
    std::cout << " -- Error, tree HGCSSTree cannot be opened. Trying RecoTree instead." << std::endl;
    isG4Tree = false;
    lTree = (TTree*)inputFile->Get("RecoTree");
    if (!lTree){
      std::cout << " -- Error, tree RecoTree cannot be opened either. Exiting..." << std::endl;
      return 1;
    }
  }

  unsigned event = 0;
  float volNb = 0;
  std::vector<HGCSSSimHit> * simhitvec = 0;
  std::vector<HGCSSRecoHit> * rechitvec = 0;

  if (isG4Tree){
    //lTree->SetBranchAddress("event",&event);
    lTree->SetBranchAddress("volNb",&volNb);
    lTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
  }
  else {
    lTree->SetBranchAddress("event",&event);
    lTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
    lTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
  }
  
  const unsigned nEvts = 
    isG4Tree ? 
    ((pNevts > lTree->GetEntries()/nLayers || pNevts==0) ? static_cast<unsigned>(lTree->GetEntries()/nLayers) : pNevts) :
    ((pNevts > lTree->GetEntries() || pNevts==0) ? static_cast<unsigned>(lTree->GetEntries()) : pNevts) ;
  
  std::cout << "- Processing = " << nEvts  << " events out of " ;
  if (isG4Tree) std::cout << lTree->GetEntries()/nLayers << std::endl;
  else std::cout << lTree->GetEntries() << std::endl;
  

  TFile *outputFile = TFile::Open("140PU/PFcal_140PU_EEHE_full.root","RECREATE");

  if (!outputFile) {
    std::cout << " -- Error, cannot open output file. Exiting..." << std::endl;
    return 1;
  }

  TTree *outputTree = new TTree("PUTree","140 PU tree");
  HGCSSSimHitVec  pusimhitvec;
  outputTree->Branch("HGCSSSimHitVec","std::vector<HGCSSSimHit>",&pusimhitvec);
  //outputTree->Branch("HGCSSRecoHitVec","std::vector<HGCSSRecoHit>",&lRecoHits);

  for (unsigned evtOut(0); evtOut<nEvtsOut;++evtOut){//loop on output events
    
    std::vector<unsigned> ipuevt;
    unsigned lVtx = 0;
    
    //get poisson <140>
    lVtx = lRndm.Poisson(140);
    ipuevt.resize(lVtx,1);
    for (unsigned iV(0); iV<lVtx; ++iV){//loop on interactions
      ipuevt[iV] = lRndm.Integer(nEvts);
      //get random PU events among available;
      std::cout << " -- PU Random number #" << iV << " for event " << evtOut << " is " << ipuevt[iV] << std::endl;
    }
    
    //loop to have all
    HGCSSSimHitVec tmpvec[nLayers];

    for (unsigned iV(0); iV<lVtx; ++iV){//loop on interactions
      std::cout << " Processing interaction " << iV << std::endl;
      unsigned ievt = isG4Tree? nLayers*ipuevt[iV] : ipuevt[iV];	
      for (unsigned iL(0); iL<nLayers; ++iL){//loop on layers
	if (!isG4Tree && iL>0) continue;

	lTree->GetEntry(ievt+iL);
	//std::cout << iL << " " << (*simhitvec).size() << std::endl;
	for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){//loop on hits
	  HGCSSSimHit lHit = (*simhitvec)[iH];
	  double posx = lHit.get_x();
	  double posy = lHit.get_y();
	  double posz = lHit.get_z();
	  

	  bool inFid = posx > minX && posx < maxX &&
	    posy > minY && posy < maxY &&
	    posz > minZ && posz < maxZ;
	  
	  if (inFid) tmpvec[iL].push_back(lHit);
	  
	}//loop on hits
      }//loop on layers
    }//loop on interactions
    for (unsigned iL(0); iL<nLayers; ++iL){
      pusimhitvec.clear();
      //std::cout << " -- Layer " << iL << " has " << tmpvec[iL].size() << " hits." << std::endl;
      pusimhitvec.reserve(tmpvec[iL].size());
      for (unsigned iH(0); iH<tmpvec[iL].size(); ++iH){
	pusimhitvec.push_back(tmpvec[iL][iH]);
      }
      outputTree->Fill();
      tmpvec[iL].clear();
    }
  }//loop on out events

  outputFile->cd();
  outputTree->Write();
  outputFile->Close();

  return 0;
  
  
}//main
