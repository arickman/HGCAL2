#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include <boost/algorithm/string.hpp>
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "boost/format.hpp"
#include "boost/function.hpp"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
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
#include "TVector3.h"

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
#include "HGCSSPUenergy.hh"

using boost::lexical_cast;
namespace po=boost::program_options;

bool testInputFile(std::string input, TFile* & file){
  file = TFile::Open(input.c_str());
  
  if (!file) {
    std::cout << " -- Error, input file " << input.c_str() << " cannot be opened. Skipping..." << std::endl;
    return false;
  }
  else std::cout << " -- input file " << file->GetName() << " successfully opened." << std::endl;
  return true;
};
int main(int argc, char** argv){//main  

  //Input output and config options
  std::string cfg;
   std::string filePath;
  std::string simFileName;
  po::options_description preconfig("Configuration"); 
  preconfig.add_options()("cfg,c",po::value<std::string>(&cfg)->required());
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(preconfig).allow_unregistered().run(), vm);
  po::notify(vm);
  po::options_description config("Configuration");
  config.add_options()
    //Input output and config options //->required()
    ("filePath,i",     po::value<std::string>(&filePath)->required())
    ("simFileName,s",  po::value<std::string>(&simFileName)->required())
    ;
  po::store(po::command_line_parser(argc, argv).options(config).allow_unregistered().run(), vm);
  po::store(po::parse_config_file<char>(cfg.c_str(), config), vm);
  po::notify(vm);
  std::string inFilePath = filePath+simFileName;
  std::cout << " -- Input parameters: " << std::endl
	    << " -- Input file path: " << filePath << std::endl;
  std::ostringstream inputsim;
  inputsim << filePath << "/" << simFileName;
  TChain *lSimTree = new TChain("HGCSSTree");
  TFile * simFile = 0;
     if (!testInputFile(inputsim.str(),simFile)) return 1;
    lSimTree->AddFile(inputsim.str().c_str());
  if (!lSimTree){
    std::cout << " -- Error, tree HGCSSTree cannot be opened. Exiting..." << std::endl;
    return 1;
  }
  //loop on events
  HGCSSEvent * event = 0;
  std::vector<HGCSSSamplingSection> * ssvec = 0;
  std::vector<HGCSSSimHit> * simhitvec = 0;
   std::vector<HGCSSGenParticle> * genvec = 0;
  lSimTree->SetBranchAddress("HGCSSEvent",&event);
  lSimTree->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
  lSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
  lSimTree->SetBranchAddress("HGCSSGenParticleVec",&genvec);

  lSimTree->GetEntry(0);
  std::cout << " -- Absorber weights used for total energy:" << std::endl;
  for(unsigned iL(0); iL<(*ssvec).size(); iL++){
    double w = (*ssvec)[iL].volX0trans()/(*ssvec)[0].volX0trans();
    std::cout << "if (layer==" << iL << ") return " << w << ";" << std::endl;
  }

  return 0;

}//main
