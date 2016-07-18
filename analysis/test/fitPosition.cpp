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
#include "Math/Vector3D.h"
#include "Math/Vector3Dfwd.h"
#include "Math/Point2D.h"
#include "Math/Point2Dfwd.h"



int main(int argc, char** argv){//main  

  if (argc < 7) {
    std::cout << " Usage: " 
	      << argv[0] << " <nEvts to process (0=all)>"
	      << " <path to input files>"
	      << " <name of input sim file>"
	      << " <name of input reco file>"
	      << " <full path to output file>"
	      << " <number of si layers to consider: 1,2 or 3>" 
      //<< " <generated E>"
	      << " <optional: debug (default=0)>"
	      << std::endl;
    return 1;
  }

  //////////////////////////////////////////////////////////
  //// Hardcoded config ////////////////////////////////////
  //////////////////////////////////////////////////////////
  bool concept = true;

  //for xvsy plots
  double minX=-1700,maxX=1700;
  double minY=-1700,maxY=1700;
  double minZ=3170,maxZ=3370;
  //double minX=-510,maxX=510;
  //double minY=-510,maxY=510;
  //double minZ=-1000,maxZ=1000;

  unsigned nX=(maxX-minX)/10,nY=(maxY-minY)/10;
  unsigned nZ=maxZ-minZ;

  //size of signal region to perform Chi2 position fit.
  //in units of 2.5mm cells to accomodate different granularities
  unsigned nSR = 12;

  //maximum value of residuals to use in error matrix: discard positions that are too far away 
  double residualMax = 25;//mm

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

  //unsigned genEn;
  //genEn = atoi(argv[7]);

  unsigned debug = 0;
  if (argc >7) debug = atoi(argv[7]);

  size_t end=outPath.find_last_of(".");
  std::string outFolder = outPath.substr(0,end);

  std::cout << " -- Input parameters: " << std::endl
	    << " -- Input file path: " << filePath << std::endl
	    << " -- Output file path: " << outPath << std::endl
	    << " -- Output folder: " << outFolder << std::endl
    //<< " -- Generated energy: " << genEn << std::endl
	    << " -- Requiring " << nSiLayers << " si layers." << std::endl
	    << " -- Number cells in signal region for fit: " << nSR << " *2.5*2.5 mm^2 cells" << std::endl
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

  //extract input energy

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

  const unsigned nLayers = myDetector.nLayers();
  const unsigned nSections = myDetector.nSections();

  std::cout << " -- N layers = " << nLayers << std::endl
	    << " -- N sections = " << nSections << std::endl;

  HGCSSGeometryConversion geomConv(inFilePath,model,cellSize);
  //assemble in 7.5*7.5 to fill maxE
  std::vector<unsigned> granularity;
  granularity.resize(nLayers,4);
  geomConv.setGranularity(granularity);
  geomConv.initialiseHistos(false,"_10");

  //////////////////////////////////////////////////
  //////////////////////////////////////////////////
  ///////// Output histos  /////////////////////////
  //////////////////////////////////////////////////
  //////////////////////////////////////////////////

  TFile *outputFile = TFile::Open(outPath.c_str(),"RECREATE");
  
  if (!outputFile) {
    std::cout << " -- Error, output file " << outPath << " cannot be opened. Please create output directory. Exiting..." << std::endl;
    return 1;
  }
  else {
    std::cout << " -- output file " << outputFile->GetName() << " successfully opened." << std::endl;
  }
    outputFile->cd();
  TH1F *p_nSimHits = new TH1F("p_nSimHits","n(SimHits)",
			      1000,0,500000);
  p_nSimHits->StatOverflows();
  
  TH1F *p_nRecHits = new TH1F("p_nRecHits","n(RecHits)",
			      1000,0,5000);
  p_nRecHits->StatOverflows();

  TH1F *p_EsimTotal = new TH1F("p_EsimTotal",";Esim (MIPs)",5000,0,50000);
  TH1F *p_ErecoTotal = new TH1F("p_ErecoTotal",";Ereco (MIPs)",5000,0,50000);
  p_EsimTotal->StatOverflows();
  p_ErecoTotal->StatOverflows();

  TH2F *p_genxy[nLayers];
  TH2F *p_xy[nLayers];
  TH2F *p_recoxy[nLayers];

  TH1F *p_residuals_x = new TH1F("p_residuals_x",";xreco-xtruth (mm)",1000,-50,50);
  TH1F *p_residuals_y = new TH1F("p_residuals_y",";yreco-ytruth (mm)",1000,-50,50);
  p_residuals_x->StatOverflows();
  p_residuals_y->StatOverflows();

  std::ostringstream lName;
  for (unsigned iL(0); iL<nLayers; ++iL){
    lName.str("");
    lName << "p_genxy_" << iL;
    p_genxy[iL] = new TH2F(lName.str().c_str(),";x(mm);y(mm)",
			   nX*10,minX,maxX,
			   nY*10,minY,maxY);
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
  }

    //////////////////////////////////////////////////
    //////////////////////////////////////////////////
    ///////// Event loop /////////////////////////////
    //////////////////////////////////////////////////
    //////////////////////////////////////////////////

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
  //initialisation for error matrix
  double mean[2][nLayers];//sum residuals for x and y
  double sigma[2][nLayers][nLayers];//sum square
  unsigned nL_mean[nLayers];//number of valid layers
  unsigned nL_sigma[nLayers][nLayers];

  std::vector<double> avgZ;
  avgZ.resize(nLayers,0);


  for (unsigned iL(0);iL<nLayers;++iL){
    EtotSim[iL] = 0;
    EtotRec[iL] = 0;
    nL_mean[iL] = 0;
    mean[0][iL] = 0;
    mean[1][iL] = 0;
    for (unsigned jL(0);jL<nLayers;++jL){
      nL_sigma[iL][jL] = 0;
      sigma[0][iL][jL] = 0;
      sigma[1][iL][jL] = 0;
    }
  }



  bool firstEvent = true;

  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
    if (debug) std::cout << "... Processing entry: " << ievt << std::endl;
    else if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;
    
    lSimTree->GetEntry(ievt);
    lRecTree->GetEntry(ievt);

    if (debug){
      std::cout << "... Size of hit vectors: sim = " <<  (*simhitvec).size() << ", reco = " << (*rechitvec).size()<< std::endl;
    }

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////
    //////// output files to save position for chi2 fit //////////
    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////

    std::ofstream fout;
    std::ostringstream foutname;
    foutname << outFolder << "_initialPos_evt" << ievt << ".dat";
    fout.open(foutname.str());
    if (!fout.is_open()){
      std::cout << " Cannot open outfile " << foutname.str() << " for writing ! Exiting..." << std::endl;
      return 1;
    }


    //////////////////////////////////////////////////
    //////////////////////////////////////////////////
    ///////// SimHits ////////////////////////////////
    //////////////////////////////////////////////////
    //////////////////////////////////////////////////

    //loop on simhits
    double etotmips = 0;
    unsigned prevLayer = 10000;
    DetectorEnum type = DetectorEnum::FECAL;
    unsigned subdetLayer=0;

    TH2F *etavsphi = new TH2F("etavsphi",";#phi;#eta;hits",150,-3.1416,3.1416,160,1.4,3.0);

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

      //unsigned sec =  myDetector.getSection(layer);

      double posx = lHit.get_x(cellSize);
      double posy = lHit.get_y(cellSize);
      double posz = lHit.get_z();
      //double radius = sqrt(posx*posx+posy*posy);
      double lRealTime = mycalib.correctTime(lHit.time(),posx,posy,posz);
      double energy = lHit.energy()*mycalib.MeVToMip(layer);

      //if (energy>1) std::cout << "Hit " << layer << " " << posx << " " << posy << " " << posz << " " << energy << std::endl;

      if (type == DetectorEnum::FECAL ||
	  type == DetectorEnum::MECAL ||
	  type == DetectorEnum::BECAL){
	//correct for si thickness
	//default for 200um
	energy *= 2./nSiLayers;
      }
      
      geomConv.fill(type,subdetLayer,energy,lRealTime,posx,posy,posz);

      bool passTime = myDigitiser.passTimeCut(type,lRealTime);
      if (!passTime) continue;

      if (debug>1) {
	std::cout << " --  SimHit " << iH << "/" << (*simhitvec).size() << " --" << std::endl
		  << " --  position x,y " << posx << "," << posy << std::endl;
	lHit.Print(std::cout);
      }

      EtotSim[layer] += energy;
      p_xy[layer]->Fill(posx,posy,energy);

      ROOT::Math::XYZVector pos(posx,posy,posz);
      etavsphi->Fill(pos.phi(),pos.eta(),energy);

      //double absweight = myDetector.subDetectorByLayer(layer).absWeight;
      double absweight = (*ssvec)[layer].volX0trans()/(*ssvec)[1].volX0trans();

      //if (versionNumber==12){
	//absweight = layer%2==0 ?
	//(*ssvec)[layer].volX0trans()/refThicknessEven : 
	//(*ssvec)[layer].volX0trans()/refThicknessOdd;
	//}
      etotmips += energy*absweight;
      
    }//loop on hits

    p_nSimHits->Fill((*simhitvec).size());
 
    if (debug)  std::cout << std::endl;


    //get position of maximum E tower
    int maxbin = etavsphi->GetMaximumBin();
    int binx,biny,binz;
    etavsphi->GetBinXYZ(maxbin,binx,biny,binz);
    double phimax =etavsphi->GetXaxis()->GetBinCenter(binx); 
    double etamax =etavsphi->GetYaxis()->GetBinCenter(biny); 

    //std::cout << " MaxE cell eta,phi = " << etamax << " " << phimax << std::endl;
    //////////////////////////////////////////////////
    //////////////////////////////////////////////////
    ///////// GenParticles////////////////////////////
    //////////////////////////////////////////////////
    //////////////////////////////////////////////////

    //get truth position

    std::vector<ROOT::Math::XYPoint> truthPos;
    truthPos.resize(nLayers,ROOT::Math::XYPoint());

    for (unsigned iP(0); iP<(*genvec).size(); ++iP){//loop on gen particles    
      //if ((*genvec).size()!= 1) (*genvec)[iP].Print(std::cout);
      if ((*genvec)[iP].trackID()==1){
	double x0 = (*genvec)[iP].x();
	double y0 = (*genvec)[iP].y();
	double z0 = (*genvec)[iP].z();
	double p = sqrt(pow((*genvec)[iP].px(),2)+pow((*genvec)[iP].py(),2)+pow((*genvec)[iP].pz(),2));
	//double energy = sqrt(pow((*genvec)[iP].mass(),2)+pow(p,2));
	//p_genxy[0]->Fill(x0,y0,energy);
	//std::cout << "init : " << x0 << " " << y0 << " " << z0 << std::endl;
	//fill layers by propagating with momentum
	ROOT::Math::XYZVector unit((*genvec)[iP].px()/p,(*genvec)[iP].py()/p,(*genvec)[iP].pz()/p);

	//std::cout << " Gen particle eta,phi = " << unit.eta() << " " << unit.phi() << std::endl;

	for (unsigned iL(0); iL<nLayers; ++iL){
	  if (avgZ[iL]<z0) avgZ[iL] = geomConv.getAverageZ(iL);
	  if (avgZ[iL]>z0) {
	    double xy = (avgZ[iL]-z0)/sinh(unit.eta());
	    double x = xy*cos(unit.phi())+x0;
	    double y = xy*sin(unit.phi())+y0;
	    
	    //std::cout << "Lay " << iL << ": " << x << " " << y << " " << avgZ[iL] << std::endl;
	    p_genxy[iL]->Fill(x,y,1);
	    truthPos[iL] = ROOT::Math::XYPoint(x,y);
	  }
	}

      }

      //p_genPartId->Fill((*genvec)[iP].pdgid());
    }//loop on gen particles

    //////////////////////////////////////////////////
    //////////////////////////////////////////////////
    ///////// RecHits ////////////////////////////////
    //////////////////////////////////////////////////
    //////////////////////////////////////////////////

    double Etotcal = 0;
    std::vector<double> xmax;
    xmax.resize(nLayers,0);
    std::vector<double> ymax;
    ymax.resize(nLayers,0);
    std::vector<double> dRmin;
    dRmin.resize(nLayers,10);

    for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on rechits
      HGCSSRecoHit lHit = (*rechitvec)[iH];
      if (debug>1) {
	std::cout << " --  RecoHit " << iH << "/" << (*rechitvec).size() << " --" << std::endl
		  << " --  position x,y " << lHit.get_x() << "," << lHit.get_y() << std::endl;
	lHit.Print(std::cout);
      }
      
      double energy = lHit.energy();//in MIP already...
      unsigned layer = lHit.layer();

      if (layer >= nLayers) {
	std::cout << " WARNING! RecoHits with layer " << layer << " outside of detector's definition range ! Please fix the digitiser or the detector definition used here. Ignoring..." << std::endl;
	continue;
      }

      double posx = lHit.get_x();
      double posy = lHit.get_y();
      double posz = lHit.get_z();
      ROOT::Math::XYZVector pos(posx,posy,posz);
      double deta = fabs(pos.eta()-etamax);
      double dphi = fabs(pos.phi()-phimax);
      double dR = sqrt(pow(deta,2)+pow(dphi,2));
      if (dR<dRmin[layer]) {
	dRmin[layer] = dR;
	xmax[layer] = posx;
	ymax[layer] = posy;
      }


      //unsigned sec =  myDetector.getSection(layer);
      
      p_recoxy[layer]->Fill(posx,posy,energy);
      EtotRec[layer] += energy;

      if (debug>1) std::cout << "-hit" << iH << "-" << layer << " " << energy << " " << EtotRec[layer];

      double absweight = (*ssvec)[layer].volX0trans()/(*ssvec)[1].volX0trans();

      Etotcal += energy*absweight;
    }//loop on rechits
    
    p_nRecHits->Fill((*rechitvec).size());

    p_EsimTotal->Fill(etotmips);
    p_ErecoTotal->Fill(Etotcal);
    
    //for (unsigned iL(0);iL<nLayers;++iL){//loop on layers
      //p_EsimvsLayer->Fill(iL,EtotSim[iL]);
      // p_ErecovsLayer->Fill(iL,EtotRec[iL]);
      //}

    //get energy-weighted position around maximum
    std::vector<ROOT::Math::XYPoint> recoPos;
    recoPos.resize(nLayers,ROOT::Math::XYPoint(0,0));
    std::vector<double> eSum;
    eSum.resize(nLayers,0);
    std::vector<unsigned> nHits;
    nHits.resize(nLayers,0);
    for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on rechits
      HGCSSRecoHit lHit = (*rechitvec)[iH];
      double energy = lHit.energy();//in MIP already...
      unsigned layer = lHit.layer();
      double posx = lHit.get_x();
      double posy = lHit.get_y();
      double posz = lHit.get_z();
      double step = cellSize*nSR/2.+0.1;//+0.1 to accomodate double precision
      if (fabs(posx-xmax[layer]) < step && 
	  fabs(posy-ymax[layer]) < step){
	recoPos[layer].SetX(recoPos[layer].X() + posx*energy);
	recoPos[layer].SetY(recoPos[layer].Y() + posy*energy);
	eSum[layer] += energy;
	if (energy>0) nHits[layer]++;
      }

    }//loop on rechits


    //fill error matrix
    if (debug) std::cout << " Summary of reco and truth positions:" << std::endl;
    for (unsigned iL(0);iL<nLayers;++iL){//loop on layers
      if (nHits[iL]==0) continue;
      recoPos[iL].SetX(recoPos[iL].X()/eSum[iL]);
      recoPos[iL].SetY(recoPos[iL].Y()/eSum[iL]);
      if (debug) std::cout << iL << " nHits=" << nHits[iL] << " Max=(" << xmax[iL] << "," << ymax[iL] << ")\t Reco=(" << recoPos[iL].X() << "," << recoPos[iL].Y() << ")\t Truth=(" << truthPos[iL].X() << "," << truthPos[iL].Y() << ")" << std::endl;
      fout << iL << " " << recoPos[iL].X() << " " << recoPos[iL].Y() << " " << avgZ[iL] << " " << truthPos[iL].X() << " " << truthPos[iL].Y() << std::endl;
    }
    for (unsigned iL(0);iL<nLayers;++iL){//loop on layers
      if (nHits[iL]==0) continue;
      double residual_xi = recoPos[iL].X()-truthPos[iL].X();
      double residual_yi = recoPos[iL].Y()-truthPos[iL].Y();
      p_residuals_x->Fill(residual_xi);
      p_residuals_y->Fill(residual_yi);
      if (fabs(residual_xi)>residualMax || fabs(residual_yi)>residualMax) continue;
      mean[0][iL] += residual_xi;
      mean[1][iL] += residual_yi;
      ++nL_mean[iL];
      for (unsigned jL(0);jL<nLayers;++jL){//loop on layers
	if (nHits[jL]==0) continue;
	double residual_xj = recoPos[jL].X()-truthPos[jL].X();
	double residual_yj = recoPos[jL].Y()-truthPos[jL].Y();
	if (fabs(residual_xj)>residualMax || fabs(residual_yj)>residualMax) continue;
	double sigma_x = residual_xi*residual_xj;
	double sigma_y = residual_yi*residual_yj;
	sigma[0][iL][jL] += sigma_x;
	sigma[1][iL][jL] += sigma_y;
	++nL_sigma[iL][jL];
      }//loop on layers
    }//loop on layers

    geomConv.initialiseHistos();
    etavsphi->Delete();

    fout.close();

    firstEvent = false;
  }//loop on entries
  std::cout << " -- Total Esim in MIPS: "
	    <<  p_EsimTotal->GetEntries() 
	    << " mean " << p_EsimTotal->GetMean() 
	    << " rms " << p_EsimTotal->GetRMS() 
	    << " rms/mean " << p_EsimTotal->GetRMS()/p_EsimTotal->GetMean()
	    << " underflows " << p_EsimTotal->GetBinContent(0)
	    << " overflows " << p_EsimTotal->GetBinContent(p_EsimTotal->GetNbinsX()+1)
	    << std::endl;
  
  std::cout << " -- Total Ereco in MIPS: "
	    <<  p_ErecoTotal->GetEntries() 
	    << " mean " << p_ErecoTotal->GetMean() 
	    << " rms " << p_ErecoTotal->GetRMS() 
	    << " rms/mean " << p_ErecoTotal->GetRMS()/p_ErecoTotal->GetMean()
	    << " underflows " << p_ErecoTotal->GetBinContent(0)
	    << " overflows " << p_ErecoTotal->GetBinContent(p_ErecoTotal->GetNbinsX()+1)
	    << std::endl;
  

  //finalise error matrix
  std::ofstream fmatrix;
  std::ostringstream fmatrixname;
  fmatrixname << outFolder << "_errorMatrix_ref.dat";
  fmatrix.open(fmatrixname.str());
  if (!fmatrix.is_open()){
    std::cout << " Cannot open outfile " << fmatrixname.str() << " for writing ! Exiting..." << std::endl;
    exit(1);
  }



  TMatrixD matrix(nLayers,nLayers);
  outputFile->cd();
  TH2F *p_errorMatrix = new TH2F("p_errorMatrix",";i;j;M_{ij}",
			    nLayers,0,nLayers,
			    nLayers,0,nLayers);
  TH2F *p_corrMatrix = new TH2F("p_corrMatrix",";i;j;M_{ij}",
				nLayers,0,nLayers,
				nLayers,0,nLayers);

  for (unsigned iL(0);iL<nLayers;++iL){//loop on layers
    mean[0][iL] = mean[0][iL]/nL_mean[iL];
    mean[1][iL] = mean[1][iL]/nL_mean[iL];
  }
  for (unsigned iL(0);iL<nLayers;++iL){//loop on layers
    for (unsigned jL(0);jL<nLayers;++jL){//loop on layers
      sigma[0][iL][jL] = sigma[0][iL][jL]/nL_sigma[iL][jL];
      sigma[1][iL][jL] = sigma[1][iL][jL]/nL_sigma[iL][jL];
      //consider average of both x and y in one matrix
      matrix[iL][jL] = 0.5*(sigma[0][iL][jL]-mean[0][iL]*mean[0][jL]+
			    sigma[1][iL][jL]-mean[1][iL]*mean[1][jL]);
      //matrix[jL][iL] = matrix[iL][jL];
      p_errorMatrix->Fill(iL,jL,matrix[iL][jL]);
 
      fmatrix << iL << " " << jL << " " << std::setprecision(17) << matrix[iL][jL] << std::endl;

      //if (iL!=jL){
      //p_matrix->Fill(jL,iL,matrix[iL][jL]);
      //}
    }
  }
  fmatrix.close();

  //fill correlation matrix
  for (unsigned iL(0);iL<nLayers;++iL){//loop on layers
    for (unsigned jL(0);jL<nLayers;++jL){//loop on layers
      if (matrix[iL][iL]!=0 && matrix[jL][jL]!= 0) 
	p_corrMatrix->Fill(iL,jL,matrix[iL][jL]/sqrt(matrix[iL][iL]*matrix[jL][jL]));
    }
  }

  //get back data for each event and perform chi2 fit:
  std::cout << " -- Performing chi2 fit for each event" << std::endl;

  unsigned nInvalidFits=0;
  TH1F *p_chi2[2];
  TH1F *p_chi2overNDF[2];

  p_chi2[0] = new TH1F("p_chi2",";#chi^{2};n_{events}",100,0,500);
  p_chi2overNDF[0] = new TH1F("p_chi2overNDF",";#chi^{2}/NDF;n_{events}",100,0,10);
  p_chi2[1] = new TH1F("p_chi2_truth",";#chi^{2};n_{events}",100,0,500);
  p_chi2overNDF[1] = new TH1F("p_chi2overNDF_truth",";#chi^{2}/NDF;n_{events}",100,0,10);
  for (unsigned rt(0); rt<2;++rt){
  p_chi2[rt]->StatOverflows();
  p_chi2overNDF[rt]->StatOverflows();
  }

  TH1F *p_impactX[2];
  p_impactX[0] = new TH1F("p_impactX",";x front face impact (mm);n_{events}",200,-500,500);
  p_impactX[1] = new TH1F("p_impactX_truth",";x front face impact (mm);n_{events}",200,-500,500);
  TH1F *p_impactY[2];
  p_impactY[0] = new TH1F("p_impactY",";y front face impact (mm);n_{events}",240,300,1500);
  p_impactY[1] = new TH1F("p_impactY_truth",";y front face impact (mm);n_{events}",240,300,1500);
  TH1F *p_angleX[2];
  p_angleX[0] = new TH1F("p_angleX",";x direction angle (rad);n_{events}",150,-3.1416,3.1416);
  p_angleX[1] = new TH1F("p_angleX_truth",";x direction angle (rad);n_{events}",150,-3.1416,3.1416);
  TH1F *p_angleY[2];
  p_angleY[0] = new TH1F("p_angleY",";y direction angle (rad);n_{events}",150,-3.1416,3.1416);
  p_angleY[1] = new TH1F("p_angleY_truth",";y direction angle (rad);n_{events}",150,-3.1416,3.1416);

  TH1F *p_positionReso[2];
  TH1F *p_angularReso[2];
  p_positionReso[0] = new TH1F("p_positionResoX",";#sigma_{x,y} (mm);n_{events}",100,0,20);
  p_positionReso[1] = new TH1F("p_positionResoY",";#sigma_{x,y} (mm);n_{events}",100,0,20);
  p_angularReso[0] = new TH1F("p_angularResoX",";#sigma_{#theta} (rad);n_{events}",100,0,1);
  p_angularReso[1] = new TH1F("p_angularResoY",";#sigma_{#theta} (rad);n_{events}",100,0,1);

  //open new file to save accurate positions
  std::ofstream fout;
  std::ostringstream foutname;
  foutname << outFolder << "_accuratePos.dat";
  fout.open(foutname.str());
  if (!fout.is_open()){
    std::cout << " Cannot open outfile " << foutname.str() << " for writing ! Exiting..." << std::endl;
    return 1;
  }

  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
    if (debug) std::cout << "... Processing entry: " << ievt << std::endl;
    else if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;

    std::ifstream fin;
    std::ostringstream finname;
    finname << outFolder << "_initialPos_evt" << ievt << ".dat";
    fin.open(finname.str());
    if (!fin.is_open()){
      std::cout << " Cannot open input file " << finname.str() << "! Exiting..." << std::endl;
      return 1;
    }

    std::vector<unsigned> layerId;
    std::vector<double> posx;
    std::vector<double> posy;
    std::vector<double> posz;
    std::vector<double> posxtruth;
    std::vector<double> posytruth;
    layerId.reserve(nLayers);
    posx.reserve(nLayers);
    posy.reserve(nLayers);
    posz.reserve(nLayers);
    posxtruth.reserve(nLayers);
    posytruth.reserve(nLayers);

    while (!fin.eof()){
      unsigned l=nLayers;
      double xr=0,yr=0,z=0,xt=0,yt=0;
      fin>>l>>xr>>yr>>z>>xt>>yt;
      if (l<nLayers){
	layerId.push_back(l);
	posx.push_back(xr);
	posy.push_back(yr);
	posz.push_back(z);
	posxtruth.push_back(xt);
	posytruth.push_back(yt);	
      }
    }

    fin.close();

    const unsigned nL = layerId.size();

    //for (unsigned iL(0); iL<nL;++iL){
      //std::cout << layerId[iL] << " " << posx[iL] << " " << posy[iL] << " " << posz[iL] << std::endl;
    //}

    //if less than 3 valid layers: no point doing a fit !!
    if (nL<3){
      nInvalidFits++;
      continue;
    }

    //number of points: x and y per layer minus number of parameters: 2 for x + 2 for y.
    double ndf = 2*nL-4;

    //Get error matrix removing lines with zero hits
    TMatrixDSym e(nL);
    TVectorD u(nL),z(nL),x(nL),y(nL);
    
    for(unsigned i(0);i<nL;++i) {
      u(i)=1.0;
      z(i)=posz[i];
      //std::cout << "fit() z(" << i << ") = " << z(i) << std::endl;
      
      for(unsigned j(i);j<nL;++j) {
	e(i,j)=matrix(layerId[i],layerId[j]);
	e(j,i)=matrix(layerId[j],layerId[i]);
      }
    }

    e.Invert();

    
    //do fit for reco and truth

    for (unsigned rt(0); rt<2;++rt){
      if (debug) {
	std::cout << "... Processing ";
	if (rt==0) std::cout << " fit to reco position.";
	else std::cout << " fit to truth position.";
	std::cout << std::endl;
      }
      double chiSq(0.0);
      double position[2];
      double positionFF[2];
      double TanAngle[2];
      
      TMatrixD fitMatrix(4,4);
      
      //resolve equation for x and y separately
      for(unsigned xy(0);xy<2;xy++) {//loop on x or y
	if (debug) {
	  std::cout << "... Processing ";
	  if (xy==0) std::cout << " fit to x position.";
	  else std::cout << " fit to y position.";
	  std::cout << std::endl;
	}
	for(unsigned i(0);i<nL;i++) {
	  x(i)= rt==0 ? ((xy==0) ? posx[i] : posy[i]) : ((xy==0) ? posxtruth[i] : posytruth[i]);
	  //std::cout << "fit() x(" << i << ") = " << x(i) << std::endl;
	}
	
	TMatrixD w(2,2);
	TVectorD v(2),p(2);
	
	w(0,0)=u*(e*u);
	w(0,1)=u*(e*z);
	w(1,0)=z*(e*u);
	w(1,1)=z*(e*z);

	v(0)=u*(e*x);
	v(1)=z*(e*x);
	
	w.Invert();
	
	p=w*v;
	if (debug) {
	  std::cout << "fit() w(0,0) = " << w(0,0) << std::endl;
	  std::cout << "fit() w(0,1) = " << w(0,1) << std::endl;
	  std::cout << "fit() w(1,0) = " << w(1,0) << std::endl;
	  std::cout << "fit() w(1,1) = " << w(1,1) << std::endl;	
	  std::cout << "fit() p(0) = " << p(0) << std::endl;
	  std::cout << "fit() p(1) = " << p(1) << std::endl;
	}

	position[xy] = p(0);
	positionFF[xy] = p(0)+p(1)*posz[0];
	TanAngle[xy] = p(1);
	
	fitMatrix[2*xy][2*xy]=w(0,0);
	fitMatrix[2*xy][2*xy+1]=w(0,1);
	fitMatrix[2*xy+1][2*xy]=w(1,0);
	fitMatrix[2*xy+1][2*xy+1]=w(1,1);
	
	
	TVectorD dp(nL);
	for(unsigned i(0);i<nL;i++) {
	  dp(i)=x(i)-p(0)-p(1)*z(i);
	}
	
	chiSq+=dp*(e*dp);
      }//loop on x or y
      
      p_chi2[rt]->Fill(chiSq);
      p_chi2overNDF[rt]->Fill(chiSq/ndf);
      p_impactX[rt]->Fill(positionFF[0]);
      p_angleX[rt]->Fill(atan(TanAngle[0]));
      p_impactY[rt]->Fill(positionFF[1]);
      p_angleY[rt]->Fill(atan(TanAngle[1]));

      if (rt==0) {
	p_positionReso[0]->Fill(sqrt(fabs(fitMatrix[0][0])));
	p_positionReso[1]->Fill(sqrt(fabs(fitMatrix[2][2])));
	p_angularReso[0]->Fill(sqrt(fabs(fitMatrix[1][1])));
	p_angularReso[1]->Fill(sqrt(fabs(fitMatrix[3][3])));

	fout << ievt << " " 
	     << position[0] << " " 
	     << sqrt(fabs(fitMatrix[0][0])) << " " 
	     << TanAngle[0] << " " 
	     << sqrt(fabs(fitMatrix[1][1])) << " "
	     << position[1] << " " 
	     << sqrt(fabs(fitMatrix[2][2])) << " "
	     << TanAngle[1] << " "
	     << sqrt(fabs(fitMatrix[3][3]))
	     << std::endl;
      }

    }//reco or truth

  }//loop on entries
  
  fout.close();    


  std::cout << " -- Number of invalid fits: " << nInvalidFits << std::endl;


  outputFile->Write();
  //outputFile->Close();
  
  return 0;


}//main
