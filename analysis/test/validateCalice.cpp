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

double getResolution(unsigned genEn,unsigned versionNumber){
  double s = 0.73;
  double c = 0.044;
  double n = 0.02;
  if (versionNumber==23){
    s = 0.48;
    c = 0.045;
    n = 0.18;
  }
  return genEn*sqrt(pow(s/sqrt(genEn),2)+pow(n/genEn,2)+pow(c,2));
};

double showerEnergy(const double Cglobal,
		    const double Etotcal,
		    const bool correctLinearity){
  double a1 = 1.06;
  double a2 = 0.001;
  double a3 = 0;//1.2e-5;
  double Esh = Cglobal*Etotcal;
  if (!correctLinearity) return Esh;
  return Esh*(a1+a2*Esh+a3*Esh*Esh);
};

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

void fillHitHistos(TH1F *p_simhitEnergy,
		   TH1F *p_simnoisehitEnergy,
		   TH1F *p_xtalknoisehitEnergy,
		   TH1F *p_rechitEnergy,
		   TH2F *p_outvsinEnergy,
		   const std::vector<double> & simE, 
		   const double & MeVToMip,
		   const double & absweight,
		   const unsigned layer,
		   HGCSSDigitisation & myDigitiser,
		   TH1F *p_noise,
		   const unsigned rebin,
		   std::vector<double> & energies,
		   std::vector<double> & hitEnergies
		   ){

  double sim = simE[0]*MeVToMip;
  energies[0] += sim*absweight;
  p_simhitEnergy->Fill(sim);//log10(simE));
  if (sim>0.5){
    energies[13] += sim*absweight;
  }

  //just noise
  double digiE = sim;
  myDigitiser.setNoise(layer,0.12);
  myDigitiser.addNoise(digiE,layer,p_noise);
  if (digiE > 0.5){
    energies[1] += digiE*absweight;
    p_simnoisehitEnergy->Fill(digiE);
  }
  digiE = sim;
  myDigitiser.setNoise(layer,0.15);
  myDigitiser.addNoise(digiE,layer,p_noise);
  if (digiE > 0.5){
    energies[2] += digiE*absweight;
  }
  digiE = sim;
  myDigitiser.setNoise(layer,0.2);
  myDigitiser.addNoise(digiE,layer,p_noise);
  if (digiE > 0.5){
    energies[3] += digiE*absweight;
  }

  //inter-pixel crosstalk+noise
  myDigitiser.setNoise(layer,0.12);
  myDigitiser.setIPCrossTalk(0.025*rebin);
  double xtalkE = myDigitiser.ipXtalk(simE)*MeVToMip;
  digiE = xtalkE;
  myDigitiser.addNoise(digiE,layer,p_noise);
  if (digiE > 0.5){
    energies[4] += digiE*absweight;
    p_xtalknoisehitEnergy->Fill(digiE);
  }
  myDigitiser.setIPCrossTalk(0.035*rebin);
  xtalkE = myDigitiser.ipXtalk(simE)*MeVToMip;
  digiE = xtalkE;
  myDigitiser.addNoise(digiE,layer,p_noise);
  if (digiE > 0.5){
    energies[5] += digiE*absweight;
  }
  myDigitiser.setIPCrossTalk(0.05*rebin);
  xtalkE = myDigitiser.ipXtalk(simE)*MeVToMip;
  digiE = xtalkE;
  myDigitiser.addNoise(digiE,layer,p_noise);
  if (digiE > 0.5){
    energies[6] += digiE*absweight;
  }

  //randomisation
  myDigitiser.setNoise(layer,0.12);
  myDigitiser.setCrossTalk(0.25);
  myDigitiser.setNTotalPixels(1156);
  myDigitiser.setSigmaPix(3);
  digiE = myDigitiser.digiE(sim);
  myDigitiser.addNoise(digiE,layer,p_noise);	      
  if (digiE > 0.5){
    energies[7] += digiE*absweight;
  }

  myDigitiser.setSigmaPix(6);
  digiE = myDigitiser.digiE(sim);
  myDigitiser.addNoise(digiE,layer,p_noise);	      
  if (digiE > 0.5){
    energies[8] += digiE*absweight;
  }

  myDigitiser.setNTotalPixels(925);
  myDigitiser.setSigmaPix(3);
  digiE = myDigitiser.digiE(sim);
  myDigitiser.addNoise(digiE,layer,p_noise);	      
  if (digiE > 0.5){
    energies[9] += digiE*absweight;
  }

  myDigitiser.setSigmaPix(6);
  digiE = myDigitiser.digiE(sim);
  myDigitiser.addNoise(digiE,layer,p_noise);	      
  if (digiE > 0.5){
    energies[10] += digiE*absweight;
  }

  myDigitiser.setIPCrossTalk(0.025*rebin);
  xtalkE = myDigitiser.ipXtalk(simE)*MeVToMip;
  
  myDigitiser.setNTotalPixels(1156);
  myDigitiser.setSigmaPix(3);
  digiE = myDigitiser.digiE(xtalkE);
  myDigitiser.addNoise(digiE,layer,p_noise);	      
  if (digiE > 0.5){
    if (layer<38) hitEnergies.push_back(digiE);
    p_rechitEnergy->Fill(digiE);
    energies[11] += digiE*absweight;
    p_outvsinEnergy->Fill(sim,digiE);
  }

  myDigitiser.setIPCrossTalk(0.035*rebin);
  xtalkE = myDigitiser.ipXtalk(simE)*MeVToMip;
  myDigitiser.setNTotalPixels(925);
  myDigitiser.setSigmaPix(6);
  myDigitiser.setNoise(layer,0.15);
  digiE = myDigitiser.digiE(xtalkE);
  myDigitiser.addNoise(digiE,layer,p_noise);	      
  if (digiE > 0.5){
    energies[12] += digiE*absweight;
  }

}

void fillBHHistos(TH2D *histE,
		  const unsigned iL,
		  bool firstEvent,
		  TRandom3 & lRndm,
		  std::map<unsigned,bool> & channelAlive,
		  TH1F *p_simhitEnergy,
		  TH1F *p_simnoisehitEnergy,
		  TH1F *p_xtalknoisehitEnergy,
		  TH1F *p_rechitEnergy,
		  TH2F *p_outvsinEnergy,
		  const double & MeVtoMip,
		  const double & absweight,
		  HGCSSDigitisation & myDigitiser,
		  TH1F *p_noise,
		  std::vector<double> & Etotcal,
		  std::vector<double> & hitEnergies,
		  double minx, double maxx,
		  double miny, double maxy
		  ){
  //std::cout << " rebin 4: " << histE->GetNbinsX() << " " << histE->GetNbinsY() << std::endl;
  for (int iX(1); iX<histE->GetNbinsX()+1;++iX){
    for (int iY(1); iY<histE->GetNbinsY()+1;++iY){
      //discard corners...
      if ((iX==1 && (iY==1 || iY==histE->GetNbinsY())) ||
	  (iX==histE->GetNbinsX() && (iY==1 || iY==histE->GetNbinsY()))
	  ) continue;
     
      double simE = histE->GetBinContent(iX,iY);
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
      
      if (posx > minx &&
	  posx < maxx && 
	  posy > miny && 
	  posy < maxy){

	//if (iL==0) std::cout << iL << " " << counter << " " << histE->GetName() << " " << iX << " " << iY << " " << posx << " " << posy << std::endl;
	
	unsigned channelId = encodeChannelId(iL,4,iX,iY);
	//discard 2% randomly
	if (firstEvent){
	  double keep = lRndm.Rndm();
	  if (keep<0.02) {
	    channelAlive[channelId] = false;
	    std::cout << " Channel " << iL << " " << posx << " " << posy << " set to dead." << std::endl;
	  }
	  else channelAlive[channelId] = true;
	}
	if (!channelAlive[channelId]) continue;
	
	fillHitHistos(p_simhitEnergy,
		      p_simnoisehitEnergy,
		      p_xtalknoisehitEnergy,
		      p_rechitEnergy,
		      p_outvsinEnergy,
		      simEvec,MeVtoMip,absweight,iL,
		      myDigitiser,p_noise,4,Etotcal,
		      hitEnergies);
      }
      
    }//iY   
  }//iX
}//method

int main(int argc, char** argv){//main  

  if (argc < 7) {
    std::cout << " Usage: " 
	      << argv[0] << " <nEvts to process (0=all)>"
	      << " <path to input files>"
	      << " <name of input sim file>"
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

  bool selectEarlyDecays = true;

  double minX=-1700,maxX=1700;
  double minY=-1700,maxY=1700;
  double minZ=3170,maxZ=5070;
  //double minX=-510,maxX=510;
  //double minY=-510,maxY=510;
  //double minZ=-1000,maxZ=1000;

  unsigned nX=(maxX-minX)/10,nY=(maxY-minY)/10;
  unsigned nZ=maxZ-minZ;

  double HcalEMCalib = 1;//39.81;//38;
  double HcalEMOffset = 0;//1.9;//-15;
  double HcalPionCalib = 1;//0.901;//1./0.9;//1/0.846;
  double HcalPionOffset = 0;//-0.81;
  // choose a jet definition
  //double R = 0.5;
  //JetDefinition jet_def(antikt_algorithm, R);

  bool plotHitSpectra = true;

  const unsigned nLimits = 10;//5;
  const double pElim[nLimits] = {3,4,5,6,7,8,9,10,11,12};
  const unsigned idxRef = 2;

  //////////////////////////////////////////////////////////
  //// End Hardcoded config ////////////////////////////////////
  //////////////////////////////////////////////////////////

  const unsigned pNevts = atoi(argv[1]);
  std::string filePath = argv[2];
  std::string simFileName = argv[3];
  std::string recoFileName = argv[4];

  std::string inFilePath = filePath+simFileName;

  std::string outPath = argv[5];
  unsigned nSiLayers = 2;
  nSiLayers = atoi(argv[6]);

  unsigned debug = 0;
  if (argc >7) debug = atoi(argv[7]);


  bool isEM = false;

  if (inFilePath.find("e-")!=inFilePath.npos || 
      inFilePath.find("e+")!=inFilePath.npos) isEM = true;

  if (selectEarlyDecays && isEM) {
    selectEarlyDecays = false;
    HcalPionCalib = 1;
    HcalPionOffset = 0;
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


  /////////////////////////////////////////////////////////////
  //Info
  /////////////////////////////////////////////////////////////

  HGCSSInfo * info=(HGCSSInfo*)simFile->Get("Info");
  const double cellSize = info->cellSize();
  const unsigned versionNumber = info->version();
  const unsigned model = info->model();
  
  //models 0,1 or 3.
  bool isTBsetup = (model != 2);
  bool isCaliceHcal = versionNumber==23;//inFilePath.find("version23")!=inFilePath.npos || inFilePath.find("version_23")!=inFilePath.npos;


  std::cout << " -- Version number is : " << versionNumber 
	    << ", model = " << model
	    << ", cellSize = " << cellSize
	    << std::endl;


  //initialise detector
  HGCSSDetector & myDetector = theDetector();

  myDetector.buildDetector(versionNumber,concept,isCaliceHcal);

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


  if (isCaliceHcal && plotHitSpectra) {
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
  }

  //fill histos for 12*12 cells
  TH2D *bottom = new TH2D("bottom",";x;y",
			  7,-450,390,
			  7,-450,390);
  TH2D *right = new TH2D("right",";x;y",
			 7,-390,450,
			 7,-450,390);
  TH2D *top = new TH2D("top",";x;y",
		       7,-390,450,
		       7,-390,450);
  TH2D *left = new TH2D("left",";x;y",
			7,-450,390,
			7,-390,450);

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

  //tree
  TTree *outtree = new TTree("Estudy","Tree to study energy resolution");
  std::vector<double> energies;
  energies.resize(14,0);
  outtree->Branch("G4",&energies[0]);
  outtree->Branch("G4mipcut",&energies[13]);
  outtree->Branch("G4Noise12",&energies[1]);
  outtree->Branch("G4Noise15",&energies[2]);
  outtree->Branch("G4Noise20",&energies[3]);
  outtree->Branch("G4XT2d5Noise",&energies[4]);
  outtree->Branch("G4XT3d5Noise",&energies[5]);
  outtree->Branch("G4XT5Noise",&energies[6]);
  outtree->Branch("G4Rand1156N3Noise",&energies[7]);
  outtree->Branch("G4Rand1156N6Noise",&energies[8]);
  outtree->Branch("G4Rand925N3Noise",&energies[9]);
  outtree->Branch("G4Rand925N6Noise",&energies[10]);
  outtree->Branch("G4XT2d5Rand1156N3Noise12",&energies[11]);
  outtree->Branch("G4XT3d5Rand925N6Noise15",&energies[12]);
  double gcCorEnergy = 0;
  std::vector<double> gcCglobal;
  gcCglobal.resize(nLimits,0);
  outtree->Branch("GlobalCorE",&gcCorEnergy);
  for (unsigned ilim(0);ilim<nLimits;++ilim){
    std::ostringstream bname;
    bname << "Cglobal_" << static_cast<unsigned>(pElim[ilim]) << "mip";
    outtree->Branch(bname.str().c_str(),&gcCglobal[ilim]);
  }

  TH2F *p_EsimvsLayer = new TH2F("p_EsimvsLayer",";layer ; Esim (MIPs)",
				 nLayers,0,nLayers,
				 1000,0,5000);
  TH2F *p_ErecovsLayer = new TH2F("p_ErecovsLayer",";layer ; Ereco (MIPs)",
				  nLayers,0,nLayers,
				  1000,0,5000);
  TH1F *p_timeSim = new TH1F("p_timeSim",";G4 time (ns)",1000,0,1000);
  TH2F *p_HCALvsECAL = new TH2F("p_HCALvsECAL",";ECAL (GeV);HCAL (GeV)",
				500,0,500,
				500,0,500);
  TH2F *p_BHCALvsFHCAL = new TH2F("p_BHCALvsFHCAL",";FHCAL (GeV);BHCAL (GeV)",
				  500,0,500,
				  500,0,500);

  TH1F *p_nGenPart = new TH1F("p_nGenPart",";n(genParticles)",200,0,200);
  TH1F *p_genPartId = new TH1F("p_genPartId",";pdgid",12000,-6000,6000);

  TH1F *p_firstInteraction = new TH1F("p_firstInteraction",";layer with 1st nucl. int.",nLayers,0,nLayers);

  TH1F *p_nSimHits = new TH1F("p_nSimHits","n(SimHits)",
			      1000,0,500000);
  p_nSimHits->StatOverflows();
  
  TH1F *p_nRecHits = new TH1F("p_nRecHits","n(RecHits)",
			      1000,0,5000);
  p_nRecHits->StatOverflows();

  TH1F *p_EsimTotal = new TH1F("p_EsimTotal",";Esim (MIPs)",15000,0,300000);
  TH1F *p_ErecoTotal = new TH1F("p_ErecoTotal",";Ereco (GeV)",20000,0,5000);
  TH1F *p_EcorTotal = new TH1F("p_EcorTotal",";Ereco (GeV)",20000,0,5000);
  p_EsimTotal->StatOverflows();
  p_ErecoTotal->StatOverflows();
  p_EcorTotal->StatOverflows();
  TH1F *p_nHitsFHCAL = new TH1F("p_nHitsFHCAL","n_{Hits} (FHCAL)",
			   1000,0,1000);
  p_nHitsFHCAL->StatOverflows();
  TH2F *p_HCALvsCglobal = new TH2F("p_HCALvsCglobal",
				   ";Cglobal (e_{lim}=5 MIP);E_{HCAL} (GeV)",
				   200,0.,2.0,
				   100,0,100);

  TH1F *p_simhitEnergy = new TH1F("p_simhitEnergy",";E_{hit} (MIPs)",200,0,10);
  TH1F *p_simnoisehitEnergy = new TH1F("p_simnoisehitEnergy",";E_{hit} (MIPs)",200,0,10);
  TH1F *p_xtalknoisehitEnergy = new TH1F("p_xtalknoisehitEnergy",";E_{hit} (MIPs)",200,0,10);
  TH1F *p_rechitEnergy = new TH1F("p_rechitEnergy",";E_{hit} (MIPs)",500,0,20);
  p_rechitEnergy->StatOverflows();
  TH2F *p_outvsinEnergy = new TH2F("p_outvsinEnergy",";E_{hit}^{sim} (MIPs);E_{hit}^{rec} (MIPs)",500,0,500,500,0,500);

  TH1F * p_noise = new TH1F("noiseCheck",";noise (MIPs)",100,-5,5);

  //TH2F *p_xy[nLayers];
  //TH2F *p_recoxy[nLayers];
  // TH1F *p_EfracSim[nLayers];
  // TH1F *p_EfracReco[nLayers];
  TH1F *p_Esim[nSections];
  TH1F *p_Ereco[nSections];

  std::ostringstream lName;
  /*for (unsigned iL(0); iL<nLayers; ++iL){
    lName.str("");
    lName << "p_xy_" << iL;
    p_xy[iL] = new TH2F(lName.str().c_str(),";x(mm);y(mm)",
   			nX,minX,maxX,
   			nY,minY,maxY);
    lName.str("");
    lName << "p_recoxy_" << iL;
    p_recoxy[iL] = new TH2F(lName.str().c_str(),";x(mm);y(mm)",
   			    nX,minX,maxX,
   			    nY,minY,maxY);
  //   lName.str("");
  //   lName << "p_EfracSim_" << iL;
  //   p_EfracSim[iL] = new TH1F(lName.str().c_str(),";integrated sim E_{layer}/E_{total}",101,0,1.01);
  //   lName.str("");
  //   lName << "p_EfracReco_" << iL;
  //   p_EfracReco[iL] = new TH1F(lName.str().c_str(),";integrated reco E_{layer}/E_{total}",101,0,1.01); 
  }
  */

  for (unsigned iD(0); iD<nSections; ++iD){
    lName.str("");
    lName << "p_Esim_" << myDetector.detName(iD);
    if (myDetector.detType(iD)==DetectorEnum::BHCAL1 || myDetector.detType(iD)==DetectorEnum::BHCAL2) p_Esim[iD] = new TH1F(lName.str().c_str(),";Esim (MIPs)",2000,0,20000);
    else p_Esim[iD] = new TH1F(lName.str().c_str(),";Esim (MIPs)",20000,0,200000);
    p_Esim[iD]->StatOverflows();
    lName.str("");
    lName << "p_Ereco_" << myDetector.detName(iD);
    if (myDetector.detType(iD)==DetectorEnum::BHCAL1 || myDetector.detType(iD)==DetectorEnum::BHCAL2)  p_Ereco[iD] = new TH1F(lName.str().c_str(),";Ereco (MIPs)",200,0,2000);
    else p_Ereco[iD] = new TH1F(lName.str().c_str(),";Ereco (MIPs)",1000,0,10000);
    p_Ereco[iD]->StatOverflows();
  }



  HGCSSEvent * event = 0;
  std::vector<HGCSSSamplingSection> * ssvec = 0;
  std::vector<HGCSSSimHit> * simhitvec = 0;
  std::vector<HGCSSRecoHit> * rechitvec = 0;
  std::vector<HGCSSGenParticle> * genvec = 0;
  
  lSimTree->SetBranchAddress("HGCSSEvent",&event);
  lSimTree->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
  lSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
  lSimTree->SetBranchAddress("HGCSSGenParticleVec",&genvec);
  
  lRecTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);

  const unsigned nEvts = ((pNevts > lSimTree->GetEntries() || pNevts==0) ? static_cast<unsigned>(lSimTree->GetEntries()) : pNevts) ;
  
  std::cout << "- Processing = " << nEvts  << " events out of " << lSimTree->GetEntries() << std::endl;
  
  //Initialise histos
  //necessary to have overflows ?
  gStyle->SetOptStat(1111111);
  double EtotSim[nLayers];
  double EtotRec[nLayers];
  
  for (unsigned iL(0);iL<nLayers;++iL){
    EtotSim[iL] = 0;
    EtotRec[iL] = 0;
  }
  double Esim[nSections];
  double Ereco[nSections];
  for (unsigned iD(0); iD<nSections; ++iD){
    Esim[iD] = 0;
    Ereco[iD] = 0;
  }

  bool firstEvent = true;
  std::map<unsigned,bool> channelAlive;

  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
    if (debug) std::cout << "... Processing entry: " << ievt << std::endl;
    else if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;
    
    lSimTree->GetEntry(ievt);
    lRecTree->GetEntry(ievt);

    if (debug){
      std::cout << "... Size of hit vectors: sim = " <<  (*simhitvec).size() << ", reco = " << (*rechitvec).size()<< std::endl;
    }

    for (unsigned iP(0); iP<(*genvec).size(); ++iP){//loop on gen particles
      p_genPartId->Fill((*genvec)[iP].pdgid());
    }//loop on gen particles

    p_nGenPart->Fill((*genvec).size());

    unsigned firstInteraction = 0;

    // double refThicknessOdd = (*ssvec)[1].volX0trans();
    // if (refThicknessOdd == 0) {
    //   std::cerr << " ERROR, ref thickness odd is " << refThicknessOdd << ", setting to 1..." << std::endl;
    //   refThicknessOdd = 1;
    // }
    // double refThicknessEven = (*ssvec)[2].volX0trans();
    // if (refThicknessEven == 0) {
    //   std::cerr << " ERROR, ref thickness odd is " << refThicknessEven << ", setting to 1..." << std::endl;
    //   refThicknessEven = 1;
    // }

    //to get simhit energy in final granularity
    unsigned prevLayer = 10000;
    DetectorEnum type = DetectorEnum::FECAL;
    unsigned subdetLayer=0;

    for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){//loop on hits
      HGCSSSimHit lHit = (*simhitvec)[iH];

      //discard some si layers...
      if (lHit.silayer() >= nSiLayers) continue; 

      unsigned layer = lHit.layer();

      if (layer >= nLayers) {
	//std::cout << " WARNING! SimHits with layer " << layer << " outside of detector's definition range ! Please fix the digitiser or the detector definition used here. Ignoring..." << std::endl;
	continue;
      }
      if (layer != prevLayer){
	const HGCSSSubDetector & subdet = myDetector.subDetectorByLayer(layer);
	type = subdet.type;
	subdetLayer = layer-subdet.layerIdMin;
	prevLayer = layer;
	if (debug > 1) std::cout << " - layer " << layer << " " << subdet.name << " " << subdetLayer << std::endl;
      }     

      unsigned sec =  myDetector.getSection(layer);

      if ( firstInteraction == 0 &&
	   (lHit.nNeutrons()>0 || 
	    lHit.nProtons()>0 ||
	    lHit.nHadrons()>0 ) && 
	   lHit.mainParentTrackID() > 0
	   ) firstInteraction = layer;

      double posx = lHit.get_x(cellSize);
      double posy = lHit.get_y(cellSize);
      double posz = lHit.get_z();
      //double radius = sqrt(posx*posx+posy*posy);
      double lRealTime = mycalib.correctTime(lHit.time(),posx,posy,posz);
      double energy = lHit.energy()*mycalib.MeVToMip(layer);

      bool passTime = myDigitiser.passTimeCut(type,lRealTime);
      if (!passTime) continue;

      //fill map to have simhits in final granularity
      if (energy>0 && isCaliceHcal && plotHitSpectra){
	geomConv.fill(type,subdetLayer,lHit.energy(),lRealTime,posx,posy,posz);
      }
      if (debug>1) {
	std::cout << " --  SimHit " << iH << "/" << (*simhitvec).size() << " --" << std::endl
		  << " --  position x,y " << posx << "," << posy << std::endl;
	lHit.Print(std::cout);
      }

      //p_xy[layer]->Fill(posx,posy,energy);
      //correct for time of flight
      p_timeSim->Fill(lRealTime);
      
      EtotSim[layer] += energy;
      if (debug>1) std::cout << "-hit" << iH << "-" << layer << " " << energy << " " << EtotSim[layer];

      double absweight = myDetector.subDetectorByLayer(layer).absWeight;
      //double absweight = (*ssvec)[layer].volX0trans()/(*ssvec)[0].volX0trans();

      //if (versionNumber==12){
	//absweight = layer%2==0 ?
	//(*ssvec)[layer].volX0trans()/refThicknessEven : 
	//(*ssvec)[layer].volX0trans()/refThicknessOdd;
      //std::cout << layer << " " << absweight << std::endl;
	//}
      Esim[sec] += energy*absweight;
      
    }//loop on hits

    for (unsigned iE(0);iE<energies.size();++iE) energies[iE] = 0;

    std::vector<double> hitEnergies;
    hitEnergies.reserve(1000);

    if (isCaliceHcal && plotHitSpectra){
      //fill hit energy in final granularity

      bool doCoarse = false;
      //for (unsigned iL(0); iL<myDetector.nLayers(DetectorEnum::FHCAL); ++iL){//loop on layers

      for (unsigned iL(0); iL<nLayers; ++iL){//loop on layers
	TH2D *histE = (TH2D*)geomConv.get2DHist(iL,"E");//->Clone();
	if (iL>=30) doCoarse = true;
	
	double MeVtoMip = mycalib.MeVToMip(iL);
	//double absweight = ((*ssvec)[iL].volX0trans())/((*ssvec)[0].volX0trans());
	double absweight = myDetector.subDetectorByLayer(iL).absWeight;
	if (debug>1) std::cout << iL << " " 
			       << (*ssvec)[iL].volX0trans() << " / " 
			       << (*ssvec)[0].volX0trans() << " = " 
			       << absweight
			       << std::endl;
	//3*3 cells
	if (!doCoarse){
	  //std::cout << iL << " rebin 1: " << histE->GetNbinsX() << " " << histE->GetNbinsY() << std::endl;
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
	      if (fabs(posx)<150 && fabs(posy)<150){
		//if (iL==0) std::cout << iL << " " << counter << " " << iX << " " << iY << " " << posx << " " << posy << std::endl;

		unsigned channelId = encodeChannelId(iL,1,iX,iY);
		//discard 2% randomly
		if (firstEvent){
		  double keep = lRndm.Rndm();
		  if (keep<0.02) {
		    channelAlive[channelId] = false;
		    std::cout << " Channel " << iL << " " << posx << " " << posy << " set to dead." << std::endl;
		  }
		  else channelAlive[channelId] = true;
		}
		if (!channelAlive[channelId]) continue;


		fillHitHistos(p_simhitEnergy,
			      p_simnoisehitEnergy,
			      p_xtalknoisehitEnergy,
			      p_rechitEnergy,
			      p_outvsinEnergy,
			      simEvec,MeVtoMip,absweight,iL,
			      myDigitiser,p_noise,1,
			      energies,
			      hitEnergies);
	      }
	    }
	  }
	}
	//6*6 cells
	histE->Rebin2D(2,2);
	//std::cout << " rebin 2: " << histE->GetNbinsX() << " " << histE->GetNbinsY() << std::endl;

	for (int iX(1); iX<histE->GetNbinsX()+1;++iX){
	  for (int iY(1); iY<histE->GetNbinsY()+1;++iY){

	    double simE = histE->GetBinContent(iX,iY);
	    std::vector<double> simEvec;
	    simEvec.push_back(simE);
	    if (iX>1) simEvec.push_back(histE->GetBinContent(iX-1,iY));
	    if (iX<histE->GetNbinsX()) simEvec.push_back(histE->GetBinContent(iX+1,iY));
	    if (iY>1) simEvec.push_back(histE->GetBinContent(iX,iY-1));
	    if (iY<histE->GetNbinsY()) simEvec.push_back(histE->GetBinContent(iX,iY+1));

	    //assemble in real granularity
	    double posx = histE->GetXaxis()->GetBinCenter(iX);
	    double posy = histE->GetYaxis()->GetBinCenter(iY);

	    bottom->Fill(posx,posy,simE);
	    right->Fill(posx,posy,simE);
	    top->Fill(posx,posy,simE);
	    left->Fill(posx,posy,simE);

	    if ( (doCoarse && fabs(posx)<330 && fabs(posy)<330) ||
		 (!doCoarse && 
		  ((fabs(posx)>150 && fabs(posx)<330 && fabs(posy)<330) ||
		   (fabs(posy)>150 && fabs(posy)<330 && fabs(posx)<330))
		  )){
	      //if (iL==0) std::cout << iL << " " << counter << " " << iX << " " << iY << " " << posx << " " << posy << std::endl;
	      unsigned channelId = encodeChannelId(iL,2,iX,iY);
	      //discard 2% randomly
	      if (firstEvent){
		double keep = lRndm.Rndm();
		if (keep<0.02) {
		  channelAlive[channelId] = false;
		  std::cout << " Channel " << iL << " " << posx << " " << posy << " set to dead." << std::endl;
		}
		else channelAlive[channelId] = true;
	      }
	      if (!channelAlive[channelId]) continue;
	      fillHitHistos(p_simhitEnergy,
			    p_simnoisehitEnergy,
			    p_xtalknoisehitEnergy,
			    p_rechitEnergy,
			    p_outvsinEnergy,
			    simEvec,MeVtoMip,absweight,iL,
			    myDigitiser,p_noise,2,energies,
			    hitEnergies);
	    }
	  }
	}

	
	//for bottom row
	fillBHHistos(bottom,iL,firstEvent,lRndm,channelAlive,
		     p_simhitEnergy,
		     p_simnoisehitEnergy,
		     p_xtalknoisehitEnergy,
		     p_rechitEnergy,
		     p_outvsinEnergy,
		     MeVtoMip,absweight,
		     myDigitiser,p_noise,energies,hitEnergies,
		     -330,270,-450,-330);
	//for right column
	fillBHHistos(right,iL,firstEvent,lRndm,channelAlive,
		     p_simhitEnergy,
		     p_simnoisehitEnergy,
		     p_xtalknoisehitEnergy,
		     p_rechitEnergy,
		     p_outvsinEnergy,
		     MeVtoMip,absweight,
		     myDigitiser,p_noise,energies,hitEnergies,
		     330,450,-330,270);
	//for top row
	fillBHHistos(top,iL,firstEvent,lRndm,channelAlive,
		     p_simhitEnergy,
		     p_simnoisehitEnergy,
		     p_xtalknoisehitEnergy,
		     p_rechitEnergy,
		     p_outvsinEnergy,
		     MeVtoMip,absweight,
		     myDigitiser,p_noise,energies,hitEnergies,
		     -270,330,330,450);
	//for left column
	fillBHHistos(left,iL,firstEvent,lRndm,channelAlive,
		     p_simhitEnergy,
		     p_simnoisehitEnergy,
		     p_xtalknoisehitEnergy,
		     p_rechitEnergy,
		     p_outvsinEnergy,
		     MeVtoMip,absweight,
		     myDigitiser,p_noise,energies,hitEnergies,
		     -450,-330,-270,330);
	

	bottom->Reset();
	right->Reset();
	top->Reset();
	left->Reset();


      }//loop on layers


      geomConv.initialiseHistos(true,false);

    }//isCaliceHcal

    p_nSimHits->Fill((*simhitvec).size());
    p_firstInteraction->Fill(firstInteraction);

    if (debug)  std::cout << std::endl;

    for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on rechits
      HGCSSRecoHit lHit = (*rechitvec)[iH];
      if (debug>1) {
	std::cout << " --  RecoHit " << iH << "/" << (*rechitvec).size() << " --" << std::endl
		  << " --  position x,y " << lHit.get_x() << "," << lHit.get_y() << std::endl;
	lHit.Print(std::cout);
      }
      
      //double posx = lHit.get_x();
      //double posy = lHit.get_y();
      //double posz = lHit.get_z();

      double energy = lHit.energy();//in MIP already...
      unsigned layer = lHit.layer();

      //p_rechitEnergy->Fill(log10(energy/mycalib.MeVToMip(layer)));

      if (layer >= nLayers) {
	//std::cout << " WARNING! RecoHits with layer " << layer << " outside of detector's definition range ! Please fix the digitiser or the detector definition used here. Ignoring..." << std::endl;
	continue;
      }
      unsigned sec =  myDetector.getSection(layer);
      
      //p_recoxy[layer]->Fill(posx,posy,energy);
      EtotRec[layer] += energy;
      if (debug>1) std::cout << "-hit" << iH << "-" << layer << " " << energy << " " << EtotRec[layer];

      Ereco[sec] += energy;
    }//loop on rechits
    
    p_nRecHits->Fill((*rechitvec).size());

    double Eecal = 0;
    if (myDetector.section(DetectorEnum::FECAL)<nSections) Eecal += myDigitiser.MIPtoGeV(myDetector.subDetectorByEnum(DetectorEnum::FECAL),Ereco[myDetector.section(DetectorEnum::FECAL)]);
    if (myDetector.section(DetectorEnum::MECAL)<nSections) Eecal += myDigitiser.MIPtoGeV(myDetector.subDetectorByEnum(DetectorEnum::MECAL),Ereco[myDetector.section(DetectorEnum::MECAL)]);
    if (myDetector.section(DetectorEnum::BECAL)<nSections) Eecal += myDigitiser.MIPtoGeV(myDetector.subDetectorByEnum(DetectorEnum::BECAL),Ereco[myDetector.section(DetectorEnum::BECAL)]);
    double Efhcal = 0;
    if (myDetector.section(DetectorEnum::FHCAL)<nSections) Efhcal += myDigitiser.MIPtoGeV(myDetector.subDetectorByEnum(DetectorEnum::FHCAL),Ereco[myDetector.section(DetectorEnum::FHCAL)]);
    double Ebhcal = 0;
    if (myDetector.section(DetectorEnum::BHCAL1)<nSections) Ebhcal += myDigitiser.MIPtoGeV(myDetector.subDetectorByEnum(DetectorEnum::BHCAL1),Ereco[myDetector.section(DetectorEnum::BHCAL1)]);
    if (myDetector.section(DetectorEnum::BHCAL2)<nSections) Ebhcal += myDigitiser.MIPtoGeV(myDetector.subDetectorByEnum(DetectorEnum::BHCAL2),Ereco[myDetector.section(DetectorEnum::BHCAL2)]);

    //double Etotcal = Eecal+(Efhcal+Ebhcal-HcalPionOffset)/HcalPionCalib;

    bool doFill = true;
    if (selectEarlyDecays && firstInteraction>5) doFill = false;


    //fill histos
    //double EtmpSim = 0;//[nSections];
    //double EtmpRec = 0;//[nSections];  
    //for (unsigned iD(0); iD<nSections; ++iD){
    //EtmpSim[iD] = 0;
    //EtmpRec[iD] = 0;
    //}
    
    double etotmips = 0;
    for (unsigned iD(0); iD<nSections; ++iD){
      etotmips += Esim[iD];//*(versionNumber==12?1:myDetector.subDetectorBySection(iD).absWeight);
    }
    
    for (unsigned iL(0);iL<nLayers;++iL){//loop on layers
      if (debug) std::cout << " -- Layer " << iL 
			   << " total sim E = " << EtotSim[iL] 
			   << " total rec E = " << EtotRec[iL] 
			   << std::endl;
      //unsigned sec =  myDetector.getSection(iL);
      if (doFill) p_EsimvsLayer->Fill(iL,EtotSim[iL]);
      if (doFill) p_ErecovsLayer->Fill(iL,EtotRec[iL]);
      //EtmpSim += EtotSim[iL];
      //EtmpRec += EtotRec[iL];
      // if (doFill) {
      // 	if (etotmips>0) p_EfracSim[iL]->Fill(EtmpSim/etotmips);
      // 	else p_EfracSim[iL]->Fill(0);
      // 	if (Etotcal>0) p_EfracReco[iL]->Fill(EtmpRec/Etotcal);
      // 	else p_EfracReco[iL]->Fill(0);
      // }
      EtotSim[iL] = 0;
      EtotRec[iL] = 0;
	
    }//loop on layers
    if (doFill) {
      p_HCALvsECAL->Fill(Eecal,Efhcal+Ebhcal);
      p_BHCALvsFHCAL->Fill(Efhcal,Ebhcal);
    }

    for (unsigned iD(0); iD<nSections; ++iD){
      if (doFill)p_Esim[iD]->Fill(Esim[iD]);
      if (doFill)p_Ereco[iD]->Fill(Ereco[iD]);
      Esim[iD]=0;
      Ereco[iD]=0;
    }
    
    if (doFill){
      double Etotcal = energies[11];
      double EtotcalEM = (Etotcal-HcalEMOffset)/HcalEMCalib;
      double EtotcalHAD = (EtotcalEM-HcalPionOffset)/HcalPionCalib;

      unsigned nHits = hitEnergies.size();
      p_nHitsFHCAL->Fill(nHits);
      //std::cout << " Check: nHits = " << nHits << std::endl;
      double efhcal = 0;
      for (unsigned iH(0); iH<nHits;++iH){
	efhcal += hitEnergies[iH];
      }
      double EmipMean = efhcal/nHits;
      //double EmipMean = p_rechitEnergy->GetMean();
      unsigned nHitsCountNum[nLimits];
      //int binMin = 0;
      //int binMax = p_rechitEnergy->FindBin(EmipMean);
      unsigned nHitsCountDen = 0;//p_rechitEnergy->Integral(binMin,binMax);
      for (unsigned iLim(0); iLim<nLimits;++iLim){
	nHitsCountNum[iLim] = 0;
	for (unsigned iH(0); iH<nHits;++iH){
	  //binMax = p_rechitEnergy->FindBin(pElim[iLim]);
	  //nHitsCountNum[iLim] = p_rechitEnergy->Integral(binMin,binMax);
	  if (hitEnergies[iH]<=pElim[iLim]) nHitsCountNum[iLim]++;
	  if (iLim==0 && hitEnergies[iH]<=EmipMean) nHitsCountDen++;
	}
	double Cglobal = 0;
	if (nHitsCountDen>0) Cglobal = nHitsCountNum[iLim]*1.0/nHitsCountDen;
	gcCglobal[iLim] = Cglobal;
	if (iLim==idxRef) {
	  gcCorEnergy = showerEnergy(Cglobal,EtotcalHAD,true);
	  p_HCALvsCglobal->Fill(Cglobal,(efhcal-10)/(38));
	}
      }
      outtree->Fill();
      p_EsimTotal->Fill(etotmips);
      p_ErecoTotal->Fill(EtotcalHAD);
      p_EcorTotal->Fill(gcCorEnergy);
    }
    
    if (firstEvent){
      std::cout << " Check of channelAlive size: " << channelAlive.size() << "/" << static_cast<unsigned>(216*30+141*24) << std::endl;
    }
    firstEvent = false;
  }//loop on entries
  
  //write
  for (unsigned iD(0); iD<nSections; ++iD){
    std::cout << " -- Summary of sim energies " 
	      << myDetector.detName(iD) << std::endl
	      <<  p_Esim[iD]->GetEntries() 
	      << " mean " << p_Esim[iD]->GetMean() 
	      << " rms " << p_Esim[iD]->GetRMS() 
	      << " underflows " << p_Esim[iD]->GetBinContent(0)
	      << " overflows " << p_Esim[iD]->GetBinContent(p_Esim[iD]->GetNbinsX()+1)
	      << std::endl;
    std::cout << " -- Summary of reco energies " 
	      << myDetector.detName(iD) << std::endl
	      <<  p_Ereco[iD]->GetEntries() 
	      << " mean " << p_Ereco[iD]->GetMean() 
	      << " rms " << p_Ereco[iD]->GetRMS() 
	      << " underflows " << p_Ereco[iD]->GetBinContent(0)
	      << " overflows " << p_Ereco[iD]->GetBinContent(p_Ereco[iD]->GetNbinsX()+1)
	      << std::endl;
  }
    std::cout << " -- Total Esim in MIPS: "
	      <<  p_EsimTotal->GetEntries() 
	      << " mean " << p_EsimTotal->GetMean() 
	      << " rms " << p_EsimTotal->GetRMS() 
	      << " rms/mean " << p_EsimTotal->GetRMS()/p_EsimTotal->GetMean()
	      << " underflows " << p_EsimTotal->GetBinContent(0)
	      << " overflows " << p_EsimTotal->GetBinContent(p_EsimTotal->GetNbinsX()+1)
	      << std::endl;

    std::cout << " -- Total Ereco in GeV: "
	      <<  p_ErecoTotal->GetEntries() 
	      << " mean " << p_ErecoTotal->GetMean() 
	      << " rms " << p_ErecoTotal->GetRMS() 
	      << " rms/mean " << p_ErecoTotal->GetRMS()/p_ErecoTotal->GetMean()
	      << " underflows " << p_ErecoTotal->GetBinContent(0)
	      << " overflows " << p_ErecoTotal->GetBinContent(p_ErecoTotal->GetNbinsX()+1)
	      << std::endl;

  outputFile->Write();
  //outputFile->Close();
  
  return 0;


}//main
