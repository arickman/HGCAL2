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

unsigned encodeChannelId(const unsigned iL,const unsigned rebin,
			 const unsigned binX,const unsigned binY)
{

  return 
    (rebin & 0x000F) | 
    ((binX & 0x00FF)<<4) |
    ((binY & 0x00FF)<<12) |
    ((iL & 0x0FFF)<<20)
    ;

}

int main(int argc, char** argv){//main  

  if (argc < 5) {
    std::cout << " Usage: " 
	      << argv[0] << " <nEvts to process (0=all)>"
	      << " <path to input files>"
	      << " <name of input sim file>"
	      << " <full path to output file>"
	      << " <optional: debug (default=0)>"
	      << std::endl;
    return 1;
  }

  const unsigned pNevts = atoi(argv[1]);
  std::string filePath = argv[2];
  std::string simFileName = argv[3];

  std::string inFilePath = filePath+simFileName;

  std::string outPath = argv[4];
  unsigned debug = 0;
  if (argc >5) debug = atoi(argv[5]);


  bool isEM = false;

  if (inFilePath.find("e-")!=inFilePath.npos || 
      inFilePath.find("e+")!=inFilePath.npos) isEM = true;



  std::cout << " -- Input parameters: " << std::endl
	    << " -- Input file path: " << filePath << std::endl
	    << " -- Output file path: " << outPath << std::endl
	    << " -- Processing ";
  if (pNevts == 0) std::cout << "all events." << std::endl;
  else std::cout << pNevts << " events." << std::endl;

  TRandom3 lRndm(1);
  std::cout << " -- Random number seed: " << lRndm.GetSeed() << std::endl;

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

  /////////////////////////////////////////////////////////////
  //Info
  /////////////////////////////////////////////////////////////

  HGCSSInfo * info=(HGCSSInfo*)simFile->Get("Info");
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

  myDetector.buildDetector(versionNumber,false,isCaliceHcal);

  //initialise calibration class
  HGCSSCalibration mycalib(inFilePath);
  HGCSSDigitisation myDigitiser;
  myDigitiser.setRandomSeed(lRndm.GetSeed());

  HGCSSGeometryConversion geomConv(inFilePath,model,cellSize);
  std::vector<unsigned> granularity;
  std::vector<double> pNoiseInMips;
 
  const unsigned nLayers = myDetector.nLayers();
  const unsigned nSections = myDetector.nSections();

  std::cout << " -- N layers = " << nLayers << std::endl
	    << " -- N sections = " << nSections << std::endl;

  granularity.resize(nLayers,12);
  pNoiseInMips.resize(nLayers,0.15);
  
  std::cout << " -- Granularities are setup like this:" << std::endl;
  for (unsigned iL(0); iL<nLayers; ++iL){
    std::cout << "Layer " ;
    if (iL<10) std::cout << " ";
    std::cout << iL << " : " << granularity[iL] << ", ";
    if (iL%5==4) std::cout << std::endl;
    myDigitiser.setNoise(iL,pNoiseInMips[iL]);
  }
  std::cout << std::endl;
  geomConv.setGranularity(granularity);
  geomConv.initialiseHistos();
  
  
  TFile *outputFile = TFile::Open(outPath.c_str(),"RECREATE");
  
  if (!outputFile) {
    std::cout << " -- Error, output file " << outPath << " cannot be opened. Please create output directory. Exiting..." << std::endl;
    return 1;
  }
  else {
    std::cout << " -- output file " << outputFile->GetName() << " successfully opened." << std::endl;
  }
  
  outputFile->cd();
    //tree
  TTree *outtree = new TTree("Estudy","Tree to study energy resolution");
  std::vector<double> energies;
  energies.resize(6,0);
  outtree->Branch("G4",&energies[0]);
  outtree->Branch("G4mipcut",&energies[1]);
  outtree->Branch("G4XT",&energies[2]);
  outtree->Branch("G4XTRandN3",&energies[3]);
  outtree->Branch("G4XTRandN3Noise",&energies[4]);
  outtree->Branch("G4XTRandN6Noise",&energies[5]);

  TH1F *p_simhitEnergy = new TH1F("p_simhitEnergy",";E_{hit} (MIPs)",200,0,10);
  TH1F *p_xtalkhitEnergy = new TH1F("p_xtalkhitEnergy",";E_{hit} (MIPs)",200,0,10);

  TH2F *p_pixvspe = new TH2F("p_pixvspe",";npe;npixels",500,0,5000,1160,0,1160);
  TH1F *p_npixels = new TH1F("p_npixels",";pixels",1160,0,1160);
  TH1F *p_npixelssmeared3 = new TH1F("p_npixelssmeared3",";pixels smeared #sigma=3",1160,0,1160);
  TH1F *p_npixelssmeared6 = new TH1F("p_npixelssmeared6",";pixels smeared #sigma=6",1160,0,1160);
  TH2F *p_outvsnpix3 =  new TH2F("p_outvsnpix3",";pixels smeared #sigma=3;E_{hit} (MIPs)",1160,0,1160,1000,0,500);
  TH2F *p_outvsnpix6 =  new TH2F("p_outvsnpix6",";pixels smeared #sigma=6;E_{hit} (MIPs)",1160,0,1160,1000,0,500);

  TH1F *p_rand3hitEnergy = new TH1F("p_rand3hitEnergy",";E_{hit} (MIPs)",200,0,10);
  TH1F *p_rand6hitEnergy = new TH1F("p_rand6hitEnergy",";E_{hit} (MIPs)",200,0,10);
  TH1F *p_rechitEnergy = new TH1F("p_rechitEnergy",";E_{hit} (MIPs)",200,0,10);

  TH2F *p_outvsinEnergy = new TH2F("p_outvsinEnergy",";E_{hit}^{sim} (MIPs);E_{hit}^{rec} (MIPs)",500,0,500,500,0,500);

  TH1F * p_noise = new TH1F("noiseCheck",";noise (MIPs)",100,-1,1);

  TH1F *p_simhitEnergy_tail = new TH1F("p_simhitEnergy_tail",";E_{hit} (MIPs)",200,0,1000);
  TH1F *p_xtalkhitEnergy_tail = new TH1F("p_xtalkhitEnergy_tail",";E_{hit} (MIPs)",200,0,1000);
  TH1F *p_rand3hitEnergy_tail = new TH1F("p_rand3hitEnergy_tail",";E_{hit} (MIPs)",200,0,1000);
  TH1F *p_rand6hitEnergy_tail = new TH1F("p_rand6hitEnergy_tail",";E_{hit} (MIPs)",200,0,1000);
  TH1F *p_rechitEnergy_tail = new TH1F("p_rechitEnergy_tail",";E_{hit} (MIPs)",200,0,1000);



  HGCSSEvent * event = 0;
  std::vector<HGCSSSamplingSection> * ssvec = 0;
  std::vector<HGCSSSimHit> * simhitvec = 0;
  
  lSimTree->SetBranchAddress("HGCSSEvent",&event);
  lSimTree->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
  lSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
    const unsigned nEvts = ((pNevts > lSimTree->GetEntries() || pNevts==0) ? static_cast<unsigned>(lSimTree->GetEntries()) : pNevts) ;
  
  std::cout << "- Processing = " << nEvts  << " events out of " << lSimTree->GetEntries() << std::endl;
  
  //Initialise histos
  //necessary to have overflows ?
  gStyle->SetOptStat(1111111);

  bool firstEvent = true;
  std::map<unsigned,bool> channelAlive;

  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
    if (debug) std::cout << "... Processing entry: " << ievt << std::endl;
    else if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;
    
    lSimTree->GetEntry(ievt);

    //to get simhit energy in final granularity
    unsigned prevLayer = 10000;
    DetectorEnum type = DetectorEnum::FECAL;
    unsigned subdetLayer=0;
    
    for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){//loop on hits
      HGCSSSimHit lHit = (*simhitvec)[iH];

      unsigned layer = lHit.layer();

      if (layer >= nLayers) {
	std::cout << " WARNING! SimHits with layer " << layer << " outside of detector's definition range ! Please fix the digitiser or the detector definition used here. Ignoring..." << std::endl;
	continue;
      }
      if (layer != prevLayer){
	const HGCSSSubDetector & subdet = myDetector.subDetectorByLayer(layer);
	type = subdet.type;
	subdetLayer = layer-subdet.layerIdMin;
	prevLayer = layer;
	if (debug > 1) std::cout << " - layer " << layer << " " << subdet.name << " " << subdetLayer << std::endl;
      }     

      double posx = lHit.get_x(cellSize);
      double posy = lHit.get_y(cellSize);
      double posz = lHit.get_z();

      double lRealTime = mycalib.correctTime(lHit.time(),posx,posy,posz);

      bool passTime = myDigitiser.passTimeCut(type,lRealTime);
      if (!passTime) continue;
      
      //fill map to have simhits in final granularity
      geomConv.fill(type,subdetLayer,lHit.energy(),lRealTime,posx,posy,posz);
      
      if (debug>1) {
	std::cout << " --  SimHit " << iH << "/" << (*simhitvec).size() << " --" << std::endl
		  << " --  position x,y " << posx << "," << posy << std::endl;
	lHit.Print(std::cout);
      }
      
    }//loop on hits

    for (unsigned iE(0);iE<energies.size();++iE) energies[iE] = 0;

    for (unsigned iL(0); iL<nLayers; ++iL){//loop on layers
      TH2D *histE = (TH2D*)geomConv.get2DHist(iL,"E");//->Clone();
      
      double MeVtoMip = mycalib.MeVToMip(iL);
      double absweight = ((*ssvec)[iL].volX0trans())/((*ssvec)[0].volX0trans());
      //3*3 cells
      
      for (int iX(1); iX<histE->GetNbinsX()+1;++iX){
	for (int iY(1); iY<histE->GetNbinsY()+1;++iY){
	  double simE = histE->GetBinContent(iX,iY);
	  //0 energy hits could still have some from cross-talk and noise...	      		
	  //if (simE==0) continue;
	  std::vector<double> simEvec;
	  simEvec.push_back(simE);
	  if (iX>1) simEvec.push_back(histE->GetBinContent(iX-1,iY));
	  if (iX<histE->GetNbinsX()) simEvec.push_back(histE->GetBinContent(iX+1,iY));
	  if (iY>1) simEvec.push_back(histE->GetBinContent(iX,iY-1));
	  if (iY<histE->GetNbinsY()) simEvec.push_back(histE->GetBinContent(iX,iY+1));
	  
	  //assemble in real granularity
	  double posx = histE->GetXaxis()->GetBinCenter(iX);
	  double posy = histE->GetYaxis()->GetBinCenter(iY);

	  unsigned channelId = encodeChannelId(iL,1,iX,iY);
	  //discard 2% randomly
	  if (firstEvent){
	    double keep = lRndm.Rndm();
	    if (keep<0.02) {
	      channelAlive[channelId] = false;
	      //std::cout << " Channel " << iL << " " << posx << " " << posy << " set to dead." << std::endl;
	    }
	    else channelAlive[channelId] = true;
	  }
	  if (!channelAlive[channelId]) continue;
	  
	  //sim
	  double sim = simEvec[0]*MeVtoMip;
	  p_simhitEnergy->Fill(sim);
	  p_simhitEnergy_tail->Fill(sim);
	  energies[0] += sim;
	  if (sim>0.5) energies[1] += sim;

	  //inter-pixel crosstalk+noise
	  myDigitiser.setIPCrossTalk(0.025);
	  double xtalkE = myDigitiser.ipXtalk(simEvec)*MeVtoMip;
	  double digiE = xtalkE;
	  p_xtalkhitEnergy->Fill(digiE);
	  p_xtalkhitEnergy_tail->Fill(digiE);
	  if (digiE>0.5) energies[2] += digiE;

	  //randomisation
	  myDigitiser.setCrossTalk(0.25);
	  myDigitiser.setNTotalPixels(1156);
	  myDigitiser.setSigmaPix(3);
	  digiE = myDigitiser.digiE(xtalkE,p_pixvspe,p_npixels,p_npixelssmeared3,p_outvsnpix3);
	  p_rand3hitEnergy->Fill(digiE);
	  p_rand3hitEnergy_tail->Fill(digiE);
	  if (digiE>0.5) energies[3] += digiE;

	  //noise
	  myDigitiser.setNoise(iL,0.15);
	  myDigitiser.addNoise(digiE,iL,p_noise);	      
	  p_rechitEnergy->Fill(digiE);
	  p_rechitEnergy_tail->Fill(digiE);
	  if (digiE>0.5) energies[4] += digiE;

	  p_outvsinEnergy->Fill(sim,digiE);

	  //sigma=6
	  myDigitiser.setSigmaPix(6);
	  TH2F *ht1 = 0;
	  TH1F *ht2 = 0;
	  digiE = myDigitiser.digiE(xtalkE,ht1,ht2,p_npixelssmeared6,p_outvsnpix6);
	  p_rand6hitEnergy->Fill(digiE);
	  p_rand6hitEnergy_tail->Fill(digiE);
	  if (digiE>0.5) energies[5] += digiE;


	}//binY
      }//binX
    }//loop on layers

    geomConv.initialiseHistos(true,false);
    outtree->Fill();

    if (firstEvent){
      std::cout << " Check of channelAlive size: " << channelAlive.size() << "/" << static_cast<unsigned>(216*30+141*24) << std::endl;
    }
    firstEvent = false;
  }//loop on entries
  

  outputFile->Write();
  //outputFile->Close();
  
  return 0;
  
  
}//main
