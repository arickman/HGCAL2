#define PI 3.14159265

#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include <cmath>
#include <boost/algorithm/string.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TProfile2D.h"
#include "TChain.h"
#include "TF1.h"

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

int main(int argc, char** argv){//main  

  if (argc < 6) {
    std::cout << " Usage: " 
	      << argv[0] << " <nEvts to process (0=all)>"
	      << " <path to input files>"
	      << " <name of input reco file>"
	      << " <full path to output file>"
	      << " <number of si layers to consider: 1,2 or 3>" 
	      << " <optional: debug (default=0)>"
	      << std::endl;
    return 1;
  }

  //////////////////////////////////////////////////////////
  //// Hardcoded config ////////////////////////////////////
  //////////////////////////////////////////////////////////
  //for HGCAL, true means only 12 FHCAL layers considered (24 are simulated)
  bool concept = true;

  bool selectEarlyDecays = false;

  double minX=-1700,maxX=1700;
  double minY=-1700,maxY=1700;
  double minZ=3170,maxZ=5070;

  unsigned nX=(maxX-minX)/10,nY=(maxY-minY)/10;
  unsigned nZ=maxZ-minZ;

  //// End Hardcoded config ////////////////////////////////////
  //////////////////////////////////////////////////////////

  const unsigned pNevts = atoi(argv[1]);
  std::string filePath = argv[2];
  std::string recoFileName = argv[3];

  std::string outPath = argv[4];
  unsigned nSiLayers = 2;
  nSiLayers = atoi(argv[5]);

  unsigned debug = 0;
  if (argc >7) debug = atoi(argv[7]);

  bool isEM = false;

  if (selectEarlyDecays && isEM) {
    selectEarlyDecays = false;
  }

  std::cout << " -- Input parameters: " << std::endl
	    << " -- Input file path: " << filePath << std::endl
	    << " -- Output file path: " << outPath << std::endl
	    << " -- Requiring " << nSiLayers << " si layers." << std::endl
	    << " -- Processing ";
  if (pNevts == 0) std::cout << "all events." << std::endl;
  else std::cout << pNevts << " events." << std::endl;

  TRandom3 lRndm(1);
  std::cout << " -- Random number seed: " << lRndm.GetSeed() << std::endl;

  /////////////////////////////////////////////////////////////
  //input
  /////////////////////////////////////////////////////////////

  TFile *recFile = new TFile(); 
  TChain *lRecTree = new TChain("RecoTree");
  if(recoFileName.find("*")!=recoFileName.npos){ 
     ofstream myscript;
     myscript.open("eosls.sh");
     myscript<<"#!/bin/bash" << std::endl;
     myscript<<"/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select ls " << filePath << std::endl; 
     myscript.close();
     FILE *script = popen("bash eosls.sh", "r");
     char eoslsName[100];
     int nF =0;
     while(fgets(eoslsName, 100, script)) {
        std::ostringstream input;
        std::string temp = std::string(eoslsName).substr(0,strlen(eoslsName)-1);
        if(temp.find("HG")!=temp.npos)continue;
        input << filePath << temp;
        lRecTree->AddFile(input.str().c_str());
        if(nF==0)recFile=TFile::Open(input.str().c_str());
        nF+=1;
     }
     pclose(script);  
     system("rm ./eosls.sh");
   }
   else {
     std::ostringstream input;
     input << filePath << recoFileName;
     lRecTree->AddFile(input.str().c_str());
     recFile=TFile::Open(input.str().c_str());
   }
   if (!recFile) {
     std::cout << " -- Error, input file cannot be opened. Exiting..." << std::endl;
     return 1;
   }

  /////////////////////////////////////////////////////////////
  //Info
  /////////////////////////////////////////////////////////////

  HGCSSInfo * info=(HGCSSInfo*)recFile->Get("Info");
  const double cellSize = info->cellSize();
  const unsigned versionNumber = info->version();
  const unsigned model = info->model();
  
  //models 0,1 or 3.
  bool isCaliceHcal = versionNumber==23;//inFilePath.find("version23")!=inFilePath.npos || inFilePath.find("version_23")!=inFilePath.npos;


  std::cout << " -- Version number is : " << versionNumber 
	    << ", model = " << model
	    << ", cellSize = " << cellSize
	    << std::endl;


  //initialise detector
  HGCSSDetector & myDetector = theDetector();

  myDetector.buildDetector(versionNumber,concept,isCaliceHcal);

  //initialise calibration class
  HGCSSDigitisation myDigitiser;
  myDigitiser.setRandomSeed(lRndm.GetSeed());

  const unsigned nLayers = myDetector.nLayers();
  const unsigned nSections = myDetector.nSections();

  std::cout << " -- N layers = " << nLayers << std::endl
	    << " -- N sections = " << nSections << std::endl;

  ofstream EnDensity;
  EnDensity.open("../test/EnergyDensity.dat");

  TFile *outputFile = TFile::Open(outPath.c_str(),"RECREATE");
  
  if (!outputFile) {
    std::cout << " -- Error, output file " << outPath << " cannot be opened. Please create output directory. Exiting..." << std::endl;
    return 1;
  }
  else {
    std::cout << " -- output file " << outputFile->GetName() << " successfully opened." << std::endl;
  }
  
  
  std::cout << " -- 2-D histograms: " << std::endl
	    << " -- X: " << nX << " " << minX << " " << maxX << std::endl
	    << " -- Y: " << nY << " " << minY << " " << maxY << std::endl
	    << " -- Z: " << nZ << " " << minZ << " " << maxZ << std::endl
    ;
  outputFile->cd();

  TH1F *p_nRecHits = new TH1F("p_nRecHits","n(RecHits)",
			      1000,0,500000);
  p_nRecHits->StatOverflows();

  TProfile *p_EvsLayer = new TProfile("p_EvsLayer","Average E vs Layer",30,1,30);

  TProfile *p_AveE[nLayers];
  std::ostringstream p_AveE_Name;
  for(unsigned iL(0); iL<nLayers; iL++){
     p_AveE_Name.str("");
     p_AveE_Name << "p_AveE_layer" <<iL+1;
     p_AveE[iL] = new TProfile(p_AveE_Name.str().c_str(),"Average Energy vs eta",40,1.5,3.5);
  }

  TH2F *p_Occupancy[nLayers];
  for(unsigned iL(0); iL<nLayers; iL++){
     p_AveE_Name.str("");
     p_AveE_Name << "p_Occupancy_Layer" <<iL+1;
     p_Occupancy[iL] = new TH2F(p_AveE_Name.str().c_str(),"RecoHits Occupancy",340,-170,170,340,-170,170);
 }

  std::vector<HGCSSRecoHit> * rechitvec = 0;
  
  lRecTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);

  const unsigned nEvts = ((pNevts > lRecTree->GetEntries() || pNevts==0) ? static_cast<unsigned>(lRecTree->GetEntries()) : pNevts) ;
  
  std::cout << "- Processing = " << nEvts  << " events out of " << lRecTree->GetEntries() << std::endl;
  
  //Initialise histos
  //necessary to have overflows ?
  gStyle->SetOptStat(0);

  double Z_layer[nLayers];
  for(unsigned iL(0);iL<nLayers;iL++){
    Z_layer[iL]=0;
  }
  
  double CellNumber[nLayers][40];
  for(unsigned iL(0);iL<nLayers;iL++){
     for(unsigned iEta(0);iEta<40;iEta++){
        CellNumber[iL][iEta] =0;
      }
  }

  std::map<unsigned,bool> channelAlive;

  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
    if (debug) std::cout << "... Processing entry: " << ievt << std::endl;
    else if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;
    
    lRecTree->GetEntry(ievt);

    double EvsLayer[nLayers][40];
    double EtotRec[nLayers];
    for(unsigned iL(0);iL<nLayers;iL++){
       EtotRec[iL]=0;
       for(unsigned iEta(0);iEta<40;iEta++){
           EvsLayer[iL][iEta] =0;
       }
    } 
  
    for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on rechits
      HGCSSRecoHit lHit = (*rechitvec)[iH];
      
      double energy = lHit.energy();//in MIP already...
      unsigned layer = lHit.layer();
      int x = lHit.get_x()/10;
      int y = lHit.get_y()/10;
      double eta = lHit.eta();
      if (layer >= nLayers) {
	//std::cout << " WARNING! RecoHits with layer " << layer << " outside of detector's definition range ! Please fix the digitiser or the detector definition used here. Ignoring..." << std::endl;
	continue;
      }

      p_Occupancy[layer]->Fill(x,y);

      if(ievt==0){
         if(Z_layer[layer] ==0)Z_layer[layer]=lHit.get_z()/10.0;} 
  
      for(unsigned iEta(0);iEta<40;iEta++){
         if(eta>1.5+iEta*0.05 && eta<1.5+(iEta+1)*0.05)EvsLayer[layer][iEta]+=energy;
      }

      if(eta>1.5 && eta<3.5)EtotRec[layer] += energy;
      if (debug>1) std::cout << "-hit" << iH << "-" << layer << " " << energy << " " << EtotRec[layer];
    }//loop on rechits
    
    p_nRecHits->Fill((*rechitvec).size());
    for(unsigned iL(0); iL<nLayers; iL++){
       p_EvsLayer->Fill(iL+1, EtotRec[iL]);
    }

    if(ievt==0){
      for(unsigned iL(0);iL<nLayers;iL++){
         for(int iX(-170); iX < 170; iX++){
            for(int iY(-170);iY < 170; iY++){
              double radius = sqrt(iX*iX + iY*iY);
              if(radius > 170 || radius < 15)continue;
              double R = sqrt(radius*radius + Z_layer[iL]*Z_layer[iL]);
              double theta = acos(Z_layer[iL]/R);
              double eta = -log(tan(theta/2.));
              for(unsigned iEta(0); iEta < 40; iEta++){
                if(eta>1.5+iEta*0.05 && eta<1.5+(iEta+1)*0.05)CellNumber[iL][iEta]+=1;
              }
             }
          }
       }
     }

    for(unsigned iEta(0);iEta<40;iEta++){
       for(unsigned iL(0);iL<nLayers;iL++){
          p_AveE[iL]->Fill(1.5+iEta*0.05,EvsLayer[iL][iEta]/CellNumber[iL][iEta]);
       }
    }
    
  }//loop on entries

  outputFile->Write();
  //outputFile->Close();

  // Fit the Energy Density
  EnDensity << "f(x) = exp(p0+p1*x)" << std::endl;
  EnDensity << "layer  p0  p1" << std::endl;
  for(unsigned iL(0);iL<nLayers;iL++){
       TF1 *expofit = new TF1("density","expo",1.5,3.0); 
       p_AveE[iL]->Fit("density","R");
       double p0 = expofit->GetParameter(0);
       double p1 = expofit->GetParameter(1);
       EnDensity << iL << " " << p0 << " " << p1 << std::endl;  
  }
  EnDensity.close();
 
  return 0;


}//main
