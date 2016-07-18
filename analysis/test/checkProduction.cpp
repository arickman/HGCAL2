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
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"

#include "HGCSSEvent.hh"
#include "HGCSSInfo.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSGenParticle.hh"
#include "HGCSSParameters.hh"
#include "HGCSSCalibration.hh"
#include "HGCSSDigitisation.hh"
#include "HGCSSDetector.hh"
#include "HGCSSGeometryConversion.hh"
#include "PositionFit.hh"

#include "Math/Vector3D.h"
#include "Math/Vector3Dfwd.h"
#include "Math/Point2D.h"
#include "Math/Point2Dfwd.h"



int main(int argc, char** argv){//main  

  if (argc < 5) {
    std::cout << " Usage: " 
	      << " <path to input files>"
	      << " <name of input sim file>"
	      << " <name of input reco file>"
	      << " <name of input reco PUmix file>"
	      << std::endl;
    return 1;
  }

  //////////////////////////////////////////////////////////
  //// End Hardcoded config ////////////////////////////////////
  //////////////////////////////////////////////////////////

  std::string filePath = argv[1];
  std::string simFileName = argv[2];
  std::string recoFileName = argv[3];
  std::string pumixFileName = argv[4];

  std::string inFilePath = filePath+simFileName;

  /////////////////////////////////////////////////////////////
  //input
  /////////////////////////////////////////////////////////////

  std::ostringstream input;
  input << filePath << "/" << simFileName;

  TFile *simFile = TFile::Open(input.str().c_str());

  if (!simFile) {
    std::cout << " -- Error, input file " << input.str() << " cannot be opened. Exiting..." << std::endl;
    return 1;
  }
  else std::cout << " -- input file " << simFile->GetName() << " successfully opened." << std::endl;
  
  TTree *lSimTree = (TTree*)simFile->Get("HGCSSTree");
  if (!lSimTree){
    std::cout << " -- Error, tree HGCSSTree cannot be opened. Exiting..." << std::endl;
    return 1;
  }

  input.str("");
  input << filePath << "/" << recoFileName;
  
  TFile *recFile = TFile::Open(input.str().c_str());

  if (!recFile) {
    std::cout << " -- Error, input file " << input.str() << " cannot be opened. Exiting..." << std::endl;
    return 1;
  }
  else std::cout << " -- input file " << recFile->GetName() << " successfully opened." << std::endl;

  TTree *lRecTree = (TTree*)recFile->Get("RecoTree");
  if (!lRecTree){
    std::cout << " -- Error, tree RecoTree cannot be opened. Exiting..." << std::endl;
    return 1;
  }

  input.str("");
  input << filePath << "/" << pumixFileName;
  
  TFile *pumixFile = TFile::Open(input.str().c_str());

  if (!pumixFile) {
    std::cout << " -- Error, input file " << input.str() << " cannot be opened. Exiting..." << std::endl;
    return 1;
  }
  else std::cout << " -- input file " << pumixFile->GetName() << " successfully opened." << std::endl;

  TTree *lPumixTree = (TTree*)pumixFile->Get("PUTree");
  if (!lPumixTree){
    std::cout << " -- Error, tree PumixoTree cannot be opened. Exiting..." << std::endl;
    return 1;
  }

  std::cout << "--- Sim evts: " << lSimTree->GetEntries()  
	    << " reco evts: " << lRecTree->GetEntries()
	    << " pumix evts: " << lPumixTree->GetEntries() << std::endl;

  
  return 0;
  

}//main
