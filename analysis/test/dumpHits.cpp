#include<string>
#include<iostream>
#include<fstream>
#include<sstream>

#include "TFile.h"
#include "TTree.h"

#include "HGCSSSimHit.hh"

int main(int argc, char** argv){//main  

  if (argc < 4) {
    std::cout << " Usage: " 
	      << argv[0] << " <nEntries to process (0=all) [1 entry = 1 layer]>"
	      << " <full path to input file: root://eoscms//eos/cms/store/user/amagnan/HGCalEEGeant4/Signal/blah.root>"
	      << " <full path to output file>"
	      << " <optional: debug (default=0)>"
	      << std::endl;
    return 1;
  }

  const unsigned pNevts = atoi(argv[1]);
  std::string filePath = argv[2];
  std::string outPath = argv[3];
  unsigned debug = 0;
  if (argc >4) debug = atoi(argv[4]);

  std::cout << " -- Input parameters: " << std::endl
	    << " -- Input file: " << filePath << std::endl
	    << " -- Output file: " << outPath << std::endl
	    << " -- Processing ";
  if (pNevts == 0) std::cout << "all events." << std::endl;
  else std::cout << pNevts << " events." << std::endl;

  TFile *inputFile = TFile::Open(filePath.c_str());
  
  if (!inputFile) {
    std::cout << " -- Error, input file " << filePath << " cannot be opened. Exiting..." << std::endl;
    return 1;
  }

  TTree *lTree = (TTree*)inputFile->Get("HGCSSTree");
  if (!lTree){
    std::cout << " -- Error, tree HGCSSTree cannot be opened. Exiting ..." << std::endl;
    return 1;
  }
  
  std::ofstream hitsOut;
  hitsOut.open(outPath);
  
  if (!hitsOut.is_open()){
    std::cout << " -- Cannot open output file for writting ! Please create directory: " << outPath << std::endl;
    return 1;
  }

  float event = 0;
  float volNb = 0;
  float volX0 = 0;
  std::vector<HGCSSSimHit> * simhitvec = 0;
  
  lTree->SetBranchAddress("event",&event);
  lTree->SetBranchAddress("volNb",&volNb);
  lTree->SetBranchAddress("volX0",&volX0);
  lTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
  
  const unsigned nEvts = 
    (pNevts > lTree->GetEntries() || pNevts==0) ? 
    lTree->GetEntries() : 
    pNevts;
  
  std::cout << "- Processing = " << nEvts  << " entries out of " ;
  std::cout << lTree->GetEntries() << std::endl;
  
  
  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
    if (debug) std::cout << "... Processing entry: " << ievt << std::endl;
    else if (ievt%1000 == 0) std::cout << "... Processing entry: " << ievt << std::endl;
      
    lTree->GetEntry(ievt);

    if (debug){
      std::cout << "... Processing layer " << volNb << " with " << (*simhitvec).size() << " simhits " << std::endl;
    }
    
    for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){//loop on hits
      HGCSSSimHit lHit = (*simhitvec)[iH];
      unsigned layer = lHit.layer();
      if (layer==volNb+1) {
	if (event==0) std::cout << " -- Warning, applying patch to layer number..." << std::endl; 
	layer = volNb;
      }

      double posx = lHit.get_x();
      double posy = lHit.get_y();
      double posz = lHit.get_z();
      double energy = lHit.energy();
      if (debug>1) {
	std::cout << " --  SimHit " << iH << "/" << (*simhitvec).size() << " --" << std::endl
		  << " --  position x,y " << posx << "," << posy << std::endl;
	lHit.Print(std::cout);
      }
      
      if (energy>0) {
	hitsOut << event << " " << layer << " " << posx << " " << posy << " " << posz << " " << energy << std::endl;
      }
	
    }//loop on hits

  }//loop on entries
  
  hitsOut.close();

  return 0;
  
}//main
