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

#include "PositionFit.hh"
#include "SignalRegion.hh"

#include "Math/Vector3D.h"
#include "Math/Vector3Dfwd.h"
#include "Math/Point2D.h"
#include "Math/Point2Dfwd.h"

#include "utilities.h"

using boost::lexical_cast;
namespace po=boost::program_options;

int main(int argc, char** argv){//main  

  //Input output and config options
  std::string cfg;
  bool concept;
  unsigned pNevts;
  std::string inFilePath;
  std::string outFilePath;
  unsigned nSiLayers;
  unsigned debug;
  unsigned start;
  unsigned end;
  double etamean;
  double deta;


  po::options_description preconfig("Configuration"); 
  preconfig.add_options()("cfg,c",po::value<std::string>(&cfg)->required());
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(preconfig).allow_unregistered().run(), vm);
  po::notify(vm);
  po::options_description config("Configuration");
  config.add_options()
    //Input output and config options //->required()
    ("concept",        po::value<bool>(&concept)->default_value(true))
    ("pNevts,n",       po::value<unsigned>(&pNevts)->default_value(0))
    ("inFilePath,i",   po::value<std::string>(&inFilePath)->required())
    ("outFilePath,o",  po::value<std::string>(&outFilePath)->required())
    ("nSiLayers",      po::value<unsigned>(&nSiLayers)->default_value(2))
    ("debug,d",        po::value<unsigned>(&debug)->default_value(0))
    ("start",        po::value<unsigned>(&start)->default_value(1))
    ("end",        po::value<unsigned>(&end)->default_value(50))
    ("etamean",po::value<double>(&etamean)->default_value(2.85))
    ("deta",po::value<double>(&deta)->default_value(0.15))
    ;

  po::store(po::command_line_parser(argc, argv).options(config).allow_unregistered().run(), vm);
  po::store(po::parse_config_file<char>(cfg.c_str(), config), vm);
  po::notify(vm);


  std::cout << " -- Input parameters: " << std::endl
	    << " -- Input file path: " << inFilePath << std::endl
	    << " -- Output file path: " << outFilePath << std::endl
	    << " -- Requiring " << nSiLayers << " si layers." << std::endl
	    << " -- Processing ";
  if (pNevts == 0) std::cout << "all events." << std::endl;
  else std::cout << pNevts << " events." << std::endl;



  /////////////////////////////////////////////////////////////
  //input
  /////////////////////////////////////////////////////////////

  TChain *puTree = new TChain("HGCSSTree");
  unsigned nPuEvts = 0;
  HGCSSInfo * info;

  /*ofstream myscript;
  myscript.open("eosls.sh");
  myscript<<"#!/bin/bash" << std::endl;
  myscript<<"/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select ls " << inFilePath << std::endl; 
  myscript.close();
  FILE *script = popen("bash eosls.sh", "r");
  char eoslsName[100];
  bool first = true;
  while(fgets(eoslsName, 100, script)) {
    std::ostringstream puInput;
    std::string temp = std::string(eoslsName).substr(0,strlen(eoslsName)-1);
    //if(temp.find("Digi")!=temp.npos)continue;
    //if(temp.find(".root")==temp.npos)continue;
    puInput << inFilePath << temp << "/gitV00-03-00/e-/HGcal_version12_model2_BOFF.root";
    puTree->AddFile(puInput.str().c_str());
    std::cout << "Adding MinBias file:" << puInput.str().c_str() << std::endl;
    if (first){
      TFile *simFile = TFile::Open(puInput.str().c_str());
      if (simFile) {
	info =(HGCSSInfo*)simFile->Get("Info");
	first=false;
      }
    }
  }
  pclose(script);  
  system("rm ./eosls.sh");*/
  bool first = true;
  std::ostringstream puInput;
  if (inFilePath.find(".root") != inFilePath.npos){
    puInput << inFilePath;
    puTree->AddFile(puInput.str().c_str());
    std::cout << "Adding file:" << puInput.str().c_str() << std::endl;
    TFile *simFile = TFile::Open(puInput.str().c_str());
    if (simFile) {
      info =(HGCSSInfo*)simFile->Get("Info");
    }
  }
  else {
    for (unsigned idx(start); idx<end+1;++idx){
      puInput.str("");
      puInput << inFilePath << "HGcal_version12_model2_BOFF_MinBias_" << idx << ".root";
      //"/pile_" << idx << "/gitV00-03-00/e-/HGcal_version12_model2_BOFF.root";
      puTree->AddFile(puInput.str().c_str());
      std::cout << "Adding MinBias file:" << puInput.str().c_str() << std::endl;
      if (first){
	TFile *simFile = TFile::Open(puInput.str().c_str());
	if (simFile) {
	  info =(HGCSSInfo*)simFile->Get("Info");
	  if (info) first=false;
	}
      }
    }
  }

  nPuEvts = puTree->GetEntries();
  std::cout << "- Number of PU events available: " << nPuEvts  << std::endl;


  /////////////////////////////////////////////////////////////
  //Info
  /////////////////////////////////////////////////////////////

  const double cellSize = info->cellSize();
  const unsigned versionNumber = info->version();
  const unsigned model = info->model();
  
  //models 0,1 or 3.
  bool isTBsetup = (model != 2);
  bool isCaliceHcal = versionNumber==23;//inFilePath.find("version23")!=inFilePath.npos || inFilePath.find("version_23")!=inFilePath.npos;

  //extract input energy

  std::cout << " -- Version number is : " << versionNumber 
	    << ", model = " << model
	    << ", cellSize = " << cellSize
	    << std::endl;


  //initialise detector
  HGCSSDetector & myDetector = theDetector();
 
  myDetector.buildDetector(versionNumber,concept,isCaliceHcal);

  HGCSSCalibration mycalib(inFilePath);

  const unsigned nLayers = myDetector.nLayers();
  const unsigned nSections = myDetector.nSections();

  std::cout << " -- N layers = " << nLayers << std::endl
	    << " -- N sections = " << nSections << std::endl;


  HGCSSGeometryConversion geomConv(inFilePath,model,cellSize);
  //set granularity to get cellsize for PU subtraction
  std::vector<unsigned> granularity;
  granularity.resize(nLayers,2);
  geomConv.setGranularity(granularity);
  geomConv.initialiseHistos();

  HGCSSDigitisation myDigitiser;
  std::vector<double> pNoiseInMips;
  pNoiseInMips.resize(nLayers,0.);
  std::vector<unsigned> pThreshInADC;
  pThreshInADC.resize(nLayers,0);

  TRandom3 lRndm(0);
  std::cout << " -- Random number seed: " << lRndm.GetSeed() << std::endl;
  myDigitiser.setRandomSeed(lRndm.GetSeed());

  //////////////////////////////////////////////////
  //////////////////////////////////////////////////
  ///////// Output File // /////////////////////////
  //////////////////////////////////////////////////
  //////////////////////////////////////////////////


  std::ostringstream label;
  label << outFilePath 
	<< "/MipStudy_eta" << (etamean-deta) << "_" << (etamean+deta) 
	<< "_MinBias" << start << "_" << end 
	<< ".root";
  
    TFile *outputFile = TFile::Open(label.str().c_str(),"RECREATE");
  
  if (!outputFile) {
    std::cout << " -- Error, output file " << label.str() << " cannot be opened. Please create output directory. Exiting..." << std::endl;
    return 1;
  }
  else {
    std::cout << " -- output file " << outputFile->GetName() << " successfully opened." << std::endl;
  }
  outputFile->cd();


  HGCSSInfo *lInfo = new HGCSSInfo();
  lInfo->cellSize(cellSize);
  lInfo->version(versionNumber);
  lInfo->version(model);

  TTree *outputTree = new TTree("RecoTree","HGC Standalone simulation reco tree");
  HGCSSRecoHitVec lRecoHits;
  HGCSSEvent lEvent;
  outputTree->Branch("HGCSSEvent",&lEvent);
  outputTree->Branch("HGCSSRecoHitVec","std::vector<HGCSSRecoHit>",&lRecoHits);
  TH1F * p_noise = new TH1F("noiseCheck",";noise (MIPs)",100,-5,5);

  unsigned maxRecHits = 0;

  const unsigned nEvts = ((pNevts > puTree->GetEntries() || pNevts==0) ? static_cast<unsigned>(puTree->GetEntries()) : pNevts) ;
  

  //loop on events
  HGCSSEvent * event = 0;
  std::vector<HGCSSSamplingSection> * ssvec = 0;
  std::vector<HGCSSSimHit> * hitvec = 0;
  std::vector<HGCSSGenParticle> * genvec = 0;

  puTree->SetBranchAddress("HGCSSEvent",&event);
  puTree->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
  puTree->SetBranchAddress("HGCSSSimHitVec",&hitvec);
  puTree->SetBranchAddress("HGCSSGenParticleVec",&genvec);


  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
    if (debug) std::cout << "... Processing entry: " << ievt << std::endl;
    else if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;

    puTree->GetEntry(ievt);
    lEvent.eventNumber(event->eventNumber());
 
    unsigned prevLayer = 10000;
    DetectorEnum type = DetectorEnum::FECAL;
    unsigned subdetLayer=0;
    for (unsigned iH(0); iH<(*hitvec).size(); ++iH){//loop on hits
      HGCSSSimHit lHit = (*hitvec)[iH];
      unsigned layer = lHit.layer();
      double leta = lHit.eta();
      double energy = lHit.energy()*mycalib.MeVToMip(layer);//,leta);

      //select specific eta bands:
      bool passeta = fabs(leta-etamean)<deta;
	//fabs(leta-1.7)<0.1 ||
	//fabs(leta-2.0)<0.1 ||
	//fabs(leta-2.5)<0.1;

      if (!passeta || energy==0) continue;

      if (layer != prevLayer){
	const HGCSSSubDetector & subdet = myDetector.subDetectorByLayer(layer);
	type = subdet.type;
	subdetLayer = layer-subdet.layerIdMin;
	prevLayer = layer;
	if (debug > 1) std::cout << " - layer " << layer << " " << subdet.name << " " << subdetLayer << std::endl;
      }
      double posx = lHit.get_x(cellSize);
      double posy = lHit.get_y(cellSize);
      double posz = zPos(layer);//lHit.get_z();
      double realtime = mycalib.correctTime(lHit.time(),posx,posy,posz);
      bool passTime = myDigitiser.passTimeCut(type,realtime);
      if (!passTime) continue;
      double lR = 0;//sqrt(pow(posx,2)+pow(posy,2));
      if (energy>0 && passeta &&
	  lHit.silayer() < geomConv.getNumberOfSiLayers(type,lR) 
	  ){
	if (debug > 1) std::cout << " hit " << iH 
				 << " lay " << layer  
				 << " x " << posx 
				 << " y " << posy
				 << " z " << posz
				 << " t " << lHit.time() << " " << realtime
				 << std::endl;
	geomConv.fill(type,subdetLayer,energy,realtime,posx,posy,posz);
      }

    }//loop on input simhits

    unsigned nTotBins = 0;
    for (unsigned iL(0); iL<nLayers; ++iL){//loop on layers
      TH2D *histE = geomConv.get2DHist(iL,"E");
      TH2D *histTime = geomConv.get2DHist(iL,"Time");
      //TH2D *histZ = geomConv.get2DHist(iL,"Z");
      const HGCSSSubDetector & subdet = myDetector.subDetectorByLayer(iL);
      DetectorEnum adet = subdet.type;
      bool isScint = subdet.isScint;
      bool isSi = subdet.isSi;
      nTotBins += histE->GetNbinsX()*histE->GetNbinsY();

      double meanZpos = zPos(iL);//geomConv.getAverageZ(iL);

      if (debug>1){
	std::cout << " -- Layer " << iL << " " << subdet.name << " z=" << meanZpos
		  << " totbins = " << nTotBins << " histE entries = " << histE->GetEntries() << std::endl;
      }

      //cell-to-cell cross-talk for scintillator
      if (isScint){
	//2.5% per 30-mm edge
	myDigitiser.setIPCrossTalk(0.025*histE->GetXaxis()->GetBinWidth(1)/30.);
      }
      else {
	myDigitiser.setIPCrossTalk(0);
      }

      for (int iX(1); iX<histE->GetNbinsX()+1;++iX){
	for (int iY(1); iY<histE->GetNbinsY()+1;++iY){
	  double digiE = 0;
	  double simE = histE->GetBinContent(iX,iY);

	  //noise=0, want to save only hits in previous eta bands
	  if (simE==0) continue;

	  //double time = 0;
	  //if (simE>0) time = histTime->GetBinContent(iX,iY)/simE;

	  //fill vector with neighbours and calculate cross-talk
	  double xtalkE = simE;
	  if (isScint){
	    std::vector<double> simEvec;
	    simEvec.push_back(simE);
	    if (iX>1) simEvec.push_back(histE->GetBinContent(iX-1,iY));
	    if (iX<histE->GetNbinsX()) simEvec.push_back(histE->GetBinContent(iX+1,iY));
	    if (iY>1) simEvec.push_back(histE->GetBinContent(iX,iY-1));
	    if (iY<histE->GetNbinsY()) simEvec.push_back(histE->GetBinContent(iX,iY+1));
	    xtalkE = myDigitiser.ipXtalk(simEvec);
	  }

	  //bool passTime = myDigitiser.passTimeCut(adet,time);
	  //if (!passTime) continue;

	  //double posz = 0;
	  //for noise only hits
	  //if (simE>0) posz = histZ->GetBinContent(iX,iY)/simE;
	  //else posz = meanZpos;

	  double x = histE->GetXaxis()->GetBinCenter(iX);
	  double y = histE->GetYaxis()->GetBinCenter(iY);

	  //if (fabs(x) > 500 || fabs(y)>500) std::cout << " x=" << x << ", y=" << y << std::endl;

	  //correct for particle angle in conversion to MIP
	  double simEcor = isTBsetup ? xtalkE : myDigitiser.mipCor(xtalkE,x,y,meanZpos);
	  digiE = simEcor;

	  if (isScint && simEcor>0) {
	    digiE = myDigitiser.digiE(simEcor);
	  }
	  myDigitiser.addNoise(digiE,iL,p_noise);

	  double noiseFrac = 1.0;
	  if (simEcor>0) noiseFrac = (digiE-simEcor)/simEcor;

	  //for silicon-based Calo
	  unsigned adc = 0;
	  //if (isSi){
	  //adc = myDigitiser.adcConverter(digiE,adet);
	  //digiE = myDigitiser.adcToMIP(adc,adet);
	  //}
	  bool aboveThresh = true;//digiE > 0.5;
	  //(isSi && adc > pThreshInADC[iL]) ||
	  //(isScint && digiE > pThreshInADC[iL]*myDigitiser.adcToMIP(1,adet,false));
	  //histE->SetBinContent(iX,iY,digiE);
	  if (aboveThresh)
	    {//save hits
	      //double calibE = myDigitiser.MIPtoGeV(subdet,digiE);
	      HGCSSRecoHit lRecHit;
	      lRecHit.layer(iL);
	      lRecHit.energy(digiE);
	      lRecHit.adcCounts(adc);
	      lRecHit.x(x);
	      lRecHit.y(y);
	      lRecHit.z(meanZpos);
	      lRecHit.noiseFraction(noiseFrac);
	      //unsigned x_cell = static_cast<unsigned>(fabs(x)/(cellSize*granularity[iL]));
	      //unsigned y_cell = static_cast<unsigned>(fabs(y)/(cellSize*granularity[iL]));
	      //lRecHit.encodeCellId(x>0,y>0,x_cell,y_cell,granularity[iL]);

	      lRecoHits.push_back(lRecHit);

	    }//save hits
	}//loop on y
      }//loop on x
    }//loop on layers

    if (debug) {
      std::cout << " **DEBUG** sim-digi-reco hits = " << (*hitvec).size() 
		<< std::endl;
    }
    
    outputTree->Fill();
    //reserve necessary space and clear vectors.
    if (lRecoHits.size() > maxRecHits) {
      maxRecHits = 2*lRecoHits.size();
      std::cout << " -- INFO: event " << ievt << " maxRecHits updated to " << maxRecHits << std::endl;
    }
    lRecoHits.clear();
    geomConv.initialiseHistos();
    lRecoHits.reserve(maxRecHits);

  }//loop on entries

  outputFile->cd();
  outputFile->WriteObjectAny(lInfo,"HGCSSInfo","Info");
  outputFile->Write();
  p_noise->Write();
  outputFile->Close();

  return 0;
  
}//main
