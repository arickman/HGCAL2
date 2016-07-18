#include<string>
#include<iostream>
#include<fstream>
#include<sstream>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include "HGCSSRecoHit.hh"
#include "utilities.h"

int main(int argc, char** argv){//main  

  if (argc < 5) {
    std::cout << " Usage: " 
	      << argv[0] << " <nEntries to process (0=all)>"
	      << " <number of runs>"
	      << " <full path to input file: root://eoscms//eos/cms/store/cmst3/group/hgcal/HGCalEEGeant4/gitV00-03-07/MinBias/DigiPu200_IC2_version12_model2_BOFF_MinBias_>"
	      << " <path to output file>"
	      << " <optional: debug (default=0)>"
	      << std::endl;
    return 1;
  }

  const unsigned nLayers = 30;
  const unsigned nSectors = 12;

  const unsigned pNevts = atoi(argv[1]);
  const unsigned nRuns = atoi(argv[2]);
  std::string filePath = argv[3];
  std::string outPath = argv[4];
  unsigned debug = 0;
  if (argc >5) debug = atoi(argv[5]);

  std::cout << " -- Input parameters: " << std::endl
	    << " -- Input file: " << filePath << std::endl
	    << " -- Output file: " << outPath << std::endl
	    << " -- Processing ";
  if (pNevts == 0) std::cout << "all events." << std::endl;
  else std::cout << pNevts << " events." << std::endl;

  //TFile *inputFile = TFile::Open(filePath.c_str());
  
  TChain *lTree = new TChain("RecoTree");
  TFile * inputFile = 0;
  if (nRuns == 0){
    if (!testInputFile(filePath,inputFile)) return 1;
    lTree->AddFile(filePath.c_str());
  }
  else {
    for (unsigned i(0);i<nRuns;++i){
      std::ostringstream lstrrec;
      lstrrec << filePath << "_" << i+1 << ".root";
      if (!testInputFile(lstrrec.str(),inputFile)) continue;
      lTree->AddFile(lstrrec.str().c_str());
    }
  }

  if (!lTree){
    std::cout << " -- Error, tree RecoTree cannot be opened. Exiting..." << std::endl;
    return 1;
  }
  
  const unsigned nT = 20;
  double threshold[nT];
  std::ofstream hitsOut[nT];

  for (unsigned it(0); it<nT;++it){
    threshold[it] = it+1;
    std::ostringstream lstrrec;
    lstrrec << outPath << "_thresh" << threshold[it] << "mips.dat";
    hitsOut[it].open(lstrrec.str());
    
    if (!hitsOut[it].is_open()){
      std::cout << " -- Cannot open output file for writting ! Please create directory: " << lstrrec.str() << std::endl;
      return 1;
    }
  }

  std::vector<HGCSSRecoHit> * rechitvec = 0;
  unsigned nPuVtx = 0;
  lTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
  if (lTree->GetBranch("nPuVtx")) lTree->SetBranchAddress("nPuVtx",&nPuVtx);

  const unsigned nEvts = 
    (pNevts > lTree->GetEntries() || pNevts==0) ? 
    lTree->GetEntries() : 
    pNevts;
  
  std::cout << "- Processing = " << nEvts  << " entries out of " ;
  std::cout << lTree->GetEntries() << std::endl;
  
  
  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
    if (debug) std::cout << "... Processing event: " << ievt << std::endl;
    else if (ievt%100 == 0) std::cout << "... Processing entry: " << ievt << std::endl;
      
    lTree->GetEntry(ievt);

    if (debug){
      std::cout << "...Number of rechits  " << (*rechitvec).size() << "." << std::endl;
    }

    unsigned nAbove[nT][nLayers][nSectors];
    for (unsigned it(0); it<nT;++it){
      for (unsigned iL(0); iL<nLayers;++iL){
	for (unsigned iS(0); iS<nSectors;++iS){
	  nAbove[it][iL][iS] = 0;
	}
      }
    }
    for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on hits
      HGCSSRecoHit lHit = (*rechitvec)[iH];
      unsigned layer = lHit.layer();

      double posx = lHit.get_x();
      double posy = lHit.get_y();

      double phi = lHit.phi();
      unsigned sector = phi>=0? 6+static_cast<unsigned>(phi/(2*TMath::Pi())*nSectors) : 5-static_cast<unsigned>(fabs(phi)/(2*TMath::Pi())*nSectors);

      if (sector==12) sector=0;
      if (sector>=nSectors) {
	std::cout << " Wrong sector value! phi=" << phi << " sector=" << sector << std::endl;
	return 1;
      }

      //if (debug>1) std::cout << " - check phi: " << phi << " sector " << sector << std::endl;

      double energy = lHit.energy();
      if (debug>1) {
	std::cout << " --  RecHit " << iH << "/" << (*rechitvec).size() << " --" << std::endl;
	lHit.Print(std::cout);
      }
      
      for (unsigned it(0); it<nT;++it){
	if (energy>threshold[it]) {
	  nAbove[it][layer][sector]++;
	}
      }
      
    }//loop on hits

    for (unsigned it(0); it<nT;++it){
      hitsOut[it] <<  nPuVtx << " ";
      for (unsigned iL(0); iL<nLayers;++iL){
	for (unsigned iS(0); iS<nSectors;++iS){
	  hitsOut[it] << nAbove[it][iL][iS] << " ";
	}
      }
      hitsOut[it]<< std::endl;
    }
    
  }//loop on entries
  
  for (unsigned it(0); it<nT;++it){
    hitsOut[it].close();
  }
  
  return 0;
  
}//main
