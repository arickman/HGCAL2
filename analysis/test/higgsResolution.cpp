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

#include "PositionFit.hh"
#include "SignalRegion.hh"
#include "HiggsMass.hh"

#include "Math/Vector3D.h"
#include "Math/Vector3Dfwd.h"
#include "Math/Point2D.h"
#include "Math/Point2Dfwd.h"
#include "Math/VectorUtil.h"
#include "TLorentzVector.h"
#include "utilities.h"

using boost::lexical_cast;
namespace po=boost::program_options;

double getCalibratedE(const std::vector<double> & Evec, const double eta){
  double Etot = 0;
  unsigned nL = Evec.size();
  for (unsigned iL(0); iL<nL;++iL){
    Etot += Evec[iL]*absWeight(iL);
  }
  //calibration for signal region 2: 3*3 cm^2
  return calibratedE(Etot,eta);
};

std::string getMatrixFolder(const double & Erec, const double & aEta){
  unsigned pt[17] = {3,5,7,10,20,30,40,50,60,70,80,90,100,125,150,175,200};
  const unsigned neta=7;
  unsigned eta[neta] = {17,19,21,23,25,27,29};
  //unsigned eta[neta] = {17,21,25,29};
  int leta = static_cast<int>(aEta*10+0.5);
  double min = 1000;
  double Egenmin = 0;
  std::ostringstream folder;
  int mineta = 1000;
  unsigned etapoint = 0;
  for (unsigned ieta(0);ieta<neta;++ieta){
    if (abs(leta-eta[ieta])<mineta){
      mineta = abs(leta-eta[ieta]);
      etapoint = eta[ieta];
    }
  }
  for (unsigned ipt(0);ipt<17;++ipt){
    if ((etapoint==17 && pt[ipt]==20) || 
	(etapoint==19 && pt[ipt]==60)) continue;
    double Egen = E(pt[ipt],leta);
    if (fabs(Erec-Egen)<min){
      min = fabs(Erec-Egen);
      folder.str("");
      folder << "eta" << etapoint << "_et" <<  pt[ipt];
      Egenmin = Egen;
    }
  }
  std::cout << " -- Found minDelta = " << min << " Egen=" << Egenmin << " Ereco=" << Erec << " etarec " << leta << " " << folder.str() << std::endl;
  return folder.str();
};

int main(int argc, char** argv){//main  

  //Input output and config options
  std::string cfg;
  bool concept;
  //size of signal region to perform Chi2 position fit.
  //in units of 2.5mm cells to accomodate different granularities
  unsigned nSR;
  //maximum value of residuals to use in error matrix: discard positions that are too far away 
  double residualMax;//mm
  unsigned pNevts;
  std::string filePath;
  unsigned nRuns;
  std::string simFileName;
  std::string recoFileName;
  std::string outPath;
  unsigned nSiLayers;
  //0:do fit+energies, 1: redo initial pos, 2: redo zpos
  unsigned redoStep;
  unsigned debug;
  bool applyPuMixFix;
  std::string singleGammaPath;
  unsigned nVtx;

  po::options_description preconfig("Configuration"); 
  preconfig.add_options()("cfg,c",po::value<std::string>(&cfg)->required());
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(preconfig).allow_unregistered().run(), vm);
  po::notify(vm);
  po::options_description config("Configuration");
  config.add_options()
    //Input output and config options //->required()
    ("concept",        po::value<bool>(&concept)->default_value(true))
    ("nSR",            po::value<unsigned>(&nSR)->default_value(3))
    ("residualMax",    po::value<double>(&residualMax)->default_value(25))
    ("pNevts,n",       po::value<unsigned>(&pNevts)->default_value(0))
    ("filePath,i",     po::value<std::string>(&filePath)->required())
    ("nRuns",        po::value<unsigned>(&nRuns)->default_value(0))
    ("simFileName,s",  po::value<std::string>(&simFileName)->required())
    ("recoFileName,r", po::value<std::string>(&recoFileName)->required())
    ("outPath,o",      po::value<std::string>(&outPath)->required())
    ("nSiLayers",      po::value<unsigned>(&nSiLayers)->default_value(2))
    ("redoStep",       po::value<unsigned>(&redoStep)->default_value(0))
    ("debug,d",        po::value<unsigned>(&debug)->default_value(0))
    ("applyPuMixFix",  po::value<bool>(&applyPuMixFix)->default_value(false))
    ("singleGammaPath",     po::value<std::string>(&singleGammaPath)->required())
    ("nVtx",        po::value<unsigned>(&nVtx)->default_value(0))
    ;

  // ("output_name,o",            po::value<std::string>(&outputname)->default_value("tmp.root"))

  po::store(po::command_line_parser(argc, argv).options(config).allow_unregistered().run(), vm);
  po::store(po::parse_config_file<char>(cfg.c_str(), config), vm);
  po::notify(vm);



  std::string inFilePath = filePath+simFileName;

  size_t end=outPath.find_last_of(".");
  std::string outFolder1 = outPath.substr(0,end)+"_gamma1";
  std::string outFolder2 = outPath.substr(0,end)+"_gamma2";


  std::cout << " -- Input parameters: " << std::endl
	    << " -- Input file path: " << filePath << std::endl
	    << " -- Output file path: " << outPath << std::endl
	    << " -- Output folders: " << outFolder1 << " " << outFolder2 << std::endl
	    << " -- Requiring " << nSiLayers << " si layers." << std::endl
	    << " -- Number cells in signal region for fit: " << nSR << " cells" << std::endl
	    << " -- Residual max considered for filling matrix and fitting: " << residualMax << " mm" << std::endl
	    << " -- Apply PUMix fix? " << applyPuMixFix << std::endl
	    << " -- Processing ";
  if (pNevts == 0) std::cout << "all events." << std::endl;
  else std::cout << pNevts << " events." << std::endl;

  TRandom3 lRndm(0);
  std::cout << " -- Random number seed: " << lRndm.GetSeed() << std::endl;

  /////////////////////////////////////////////////////////////
  //input
  /////////////////////////////////////////////////////////////

  std::ostringstream inputsim;
  inputsim << filePath << "/" << simFileName;
  std::ostringstream inputrec;
  inputrec << filePath << "/" << recoFileName;

  //std::cout << inputsim.str() << " " << inputrec.str() << std::endl;

  HGCSSInfo * info;

  TChain *lSimTree = new TChain("HGCSSTree");
  TChain *lRecTree = 0;
  
  TFile * simFile = 0;
  TFile * recFile = 0;

  if (recoFileName.find("Digi") != recoFileName.npos) 
    lRecTree = new TChain("RecoTree");
  else lRecTree = new TChain("PUTree");

  if (nRuns == 0){
    if (!testInputFile(inputsim.str(),simFile)) return 1;
    if (simFile) info =(HGCSSInfo*)simFile->Get("Info");
    else {
      std::cout << " -- Error in getting information from simfile!" << std::endl;
      return 1;
    }
    lSimTree->AddFile(inputsim.str().c_str());
    if (!testInputFile(inputrec.str(),recFile)) return 1;
    lRecTree->AddFile(inputrec.str().c_str());
  }
  else {
    for (unsigned i(0);i<nRuns;++i){
      std::ostringstream lstr;
      lstr << inputsim.str() << "_run" << i << ".root";
      if (testInputFile(lstr.str(),simFile)){  
	if (simFile) info =(HGCSSInfo*)simFile->Get("Info");
	else {
	  std::cout << " -- Error in getting information from simfile!" << std::endl;
	  return 1;
	}
      }
      else continue;
      lSimTree->AddFile(lstr.str().c_str());
      lstr.str("");
      lstr << inputrec.str() << "_run" << i << ".root";
      if (!testInputFile(lstr.str(),recFile)) continue;
      lRecTree->AddFile(lstr.str().c_str());
    }
  }

  if (!lSimTree){
    std::cout << " -- Error, tree HGCSSTree cannot be opened. Exiting..." << std::endl;
    return 1;
  }

  if (!lRecTree){
    std::cout << " -- Error, tree RecoTree cannot be opened. Exiting..." << std::endl;
    return 1;
  }


  /////////////////////////////////////////////////////////////
  //Info
  /////////////////////////////////////////////////////////////

  const double cellSize = info->cellSize();
  const unsigned versionNumber = info->version();
  const unsigned model = info->model();
  
  //models 0,1 or 3.
  //bool isTBsetup = (model != 2);
  bool isCaliceHcal = versionNumber==23;//inFilePath.find("version23")!=inFilePath.npos || inFilePath.find("version_23")!=inFilePath.npos;

  //extract input energy

  std::cout << " -- Version number is : " << versionNumber 
	    << ", model = " << model
	    << ", cellSize = " << cellSize
	    << std::endl;


  //initialise detector
  HGCSSDetector & myDetector = theDetector();
 
  //bypass radius
  myDetector.buildDetector(versionNumber,concept,isCaliceHcal,true);

  const unsigned nLayers = myDetector.nLayers();
  const unsigned nSections = myDetector.nSections();

  std::cout << " -- N layers = " << nLayers << std::endl
	    << " -- N sections = " << nSections << std::endl;


  HGCSSGeometryConversion geomConv(inFilePath,model,cellSize);
  //set granularity to get cellsize for PU subtraction
  std::vector<unsigned> granularity;
  granularity.resize(nLayers,4);
  geomConv.setGranularity(granularity);

  //////////////////////////////////////////////////
  //////////////////////////////////////////////////
  ///////// Output File // /////////////////////////
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


  ///initialise PU density object

  HGCSSPUenergy puDensity("data/EnergyDensity.dat");


    //////////////////////////////////////////////////
    //////////////////////////////////////////////////
    ///////// positionFit /////////////////////////////
    //////////////////////////////////////////////////
    //////////////////////////////////////////////////
  
  const unsigned nEvts = ((pNevts > lSimTree->GetEntries() || pNevts==0) ? static_cast<unsigned>(lSimTree->GetEntries()) : pNevts) ;

  //higgs mass
  
  HiggsMass hM;
  hM.initialiseHistograms(outputFile,"HiggsMass");
  HiggsMass hMnofit;
  hMnofit.initialiseHistograms(outputFile,"HiggsMassNoFit");
  //Photon 1
  PositionFit lGamma1(nSR,residualMax,nLayers,nSiLayers,applyPuMixFix,debug,false);
  lGamma1.initialise(outputFile,"Gamma1Fit",outFolder1,geomConv,puDensity);

  if (!lGamma1.getZpositions(versionNumber))
    lGamma1.getZpositions(versionNumber,lSimTree,nEvts);
  
  
  //Photon 2
  PositionFit lGamma2(nSR,residualMax,nLayers,nSiLayers,applyPuMixFix,debug,false);
  lGamma2.initialise(outputFile,"Gamma2Fit",outFolder2,geomConv,puDensity);

  if (!lGamma2.getZpositions(versionNumber))
    lGamma2.getZpositions(versionNumber,lSimTree,nEvts);

  std::vector<double> zPos = lGamma1.getZpositions();


  bool doFit1 = false;
  bool doFit2 = false;
  //initialise
  if (redoStep>0) {
    //lGamma1.getInitialPositions(lSimTree,lRecTree,nEvts,1);
    //lGamma2.getInitialPositions(lSimTree,lRecTree,nEvts,2);
    lGamma1.initialiseClusterHistograms();
    lGamma1.initialisePositionHistograms();
    lGamma2.initialiseClusterHistograms();
    lGamma2.initialisePositionHistograms();
    lGamma1.initialiseLeastSquareFit();
    lGamma2.initialiseLeastSquareFit();
    doFit1 = true;
    doFit2 = true;
  }
  
  //initialise signal regions
  bool doOld = false;
  if (doOld) {
    doFit1 = true;
    doFit2 = true;
  }
  SignalRegion Signal1(outFolder1, nLayers, nEvts, geomConv, puDensity,applyPuMixFix,versionNumber);
  Signal1.initialise(outputFile,"Gamma1E");
  SignalRegion Signal1nofit(outFolder1, nLayers, nEvts, geomConv, puDensity,applyPuMixFix,versionNumber);
  Signal1nofit.initialise(outputFile,"Gamma1Enofit");
  if (!doOld && !Signal1.initialiseFitPositions()) {
    std::cout << " -- Redo fit for photon 1" << std::endl;
    doFit1 = true;
    if (redoStep==0) {
      lGamma1.initialiseClusterHistograms();
      lGamma1.initialisePositionHistograms();
      lGamma1.initialiseLeastSquareFit();
    }
  }

  SignalRegion Signal2(outFolder2, nLayers, nEvts, geomConv, puDensity,applyPuMixFix,versionNumber);
  Signal2.initialise(outputFile,"Gamma2E");
  SignalRegion Signal2nofit(outFolder2, nLayers, nEvts, geomConv, puDensity,applyPuMixFix,versionNumber);
  Signal2nofit.initialise(outputFile,"Gamma2Enofit");
  if (!doOld && !Signal2.initialiseFitPositions()) {
    std::cout << " -- Redo fit for photon 2" << std::endl;
    doFit2 = true;
    if (redoStep==0){
      lGamma2.initialiseClusterHistograms();
      lGamma2.initialisePositionHistograms();
      lGamma2.initialiseLeastSquareFit();
    }
  }

  unsigned nTwoPhotons = 0;
  unsigned nMatrixNotFound1 = 0;
  unsigned nMatrixNotFound2 = 0;
  unsigned nTooFar1 = 0;
  unsigned nTooFar2 = 0;

  std::cout << " --- Number of events: " << nEvts << std::endl;

  //loop on events
  HGCSSEvent * event = 0;
  std::vector<HGCSSSamplingSection> * ssvec = 0;
  std::vector<HGCSSSimHit> * simhitvec = 0;
  std::vector<HGCSSRecoHit> * rechitvec = 0;
  std::vector<HGCSSGenParticle> * genvec = 0;
  unsigned nPuVtx = 0;

  lSimTree->SetBranchAddress("HGCSSEvent",&event);
  lSimTree->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
  lSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
  lSimTree->SetBranchAddress("HGCSSGenParticleVec",&genvec);

  lRecTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
  if (lRecTree->GetBranch("nPuVtx")) lRecTree->SetBranchAddress("nPuVtx",&nPuVtx);


  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
    if (debug) std::cout << "... Processing entry: " << ievt << std::endl;
    else if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;

    lSimTree->GetEntry(ievt);
    lRecTree->GetEntry(ievt);

    if (!lGamma1.setTruthInfo(genvec,1)) continue;
    const Direction & truthDir1 = lGamma1.truthDir();
    const ROOT::Math::XYZPoint & truthVtx1 = lGamma1.truthVtx();
    const double truthE1 = lGamma1.truthE();

    if (!lGamma2.setTruthInfo(genvec,2)) continue;
    const Direction & truthDir2 = lGamma2.truthDir();
    const ROOT::Math::XYZPoint & truthVtx2 = lGamma2.truthVtx();
    const double truthE2 = lGamma2.truthE();

    bool found1 = false;
    bool found2 = false;
    Direction recoDir1;
    Direction recoDir2;
    Direction recoDir1nofit;
    Direction recoDir2nofit;
    ROOT::Math::XYZPoint posFF1nofit;
    ROOT::Math::XYZPoint posFF2nofit;
    FitResult fit1;
    FitResult fit2;

    if (debug) std::cout << " dofit = " << doFit1 << " " << doFit2 << std::endl;

    if (doFit1 || doFit2){
      //get initial position
      bool good1 = lGamma1.getInitialPosition(ievt,nPuVtx,rechitvec,nTooFar1);
      bool good2 = lGamma2.getInitialPosition(ievt,nPuVtx,rechitvec,nTooFar2);

      if (good1) lGamma1.getOutTree()->Fill();
      if (good2) lGamma2.getOutTree()->Fill();
      
      if (!good1 || !good2) continue;

      //get first guess at energy
      std::vector<unsigned> layerId1;
      std::vector<double> posx1;
      std::vector<double> posy1;
      std::vector<double> posz1;
      std::vector<double> posxtruth1;
      std::vector<double> posytruth1;
      std::vector<unsigned> layerId2;
      std::vector<double> posx2;
      std::vector<double> posy2;
      std::vector<double> posz2;
      std::vector<double> posxtruth2;
      std::vector<double> posytruth2;
      std::vector<double> Ereco1;
      std::vector<double> Ereco2;
      layerId1.reserve(nLayers);
      posx1.reserve(nLayers);
      posy1.reserve(nLayers);
      posz1.reserve(nLayers);
      posxtruth1.reserve(nLayers);
      posytruth1.reserve(nLayers);
      layerId2.reserve(nLayers);
      posx2.reserve(nLayers);
      posy2.reserve(nLayers);
      posz2.reserve(nLayers);
      posxtruth2.reserve(nLayers);
      posytruth2.reserve(nLayers);
      Ereco1.reserve(nLayers);
      Ereco2.reserve(nLayers);
      if (!lGamma1.getPositionFromFile(ievt,
				      layerId1,posx1,posy1,posz1,
				      posxtruth1,posytruth1,
				      Ereco1,
				      true,false) ||
	  !lGamma2.getPositionFromFile(ievt,
				      layerId2,posx2,posy2,posz2,
				      posxtruth2,posytruth2,
				      Ereco2,
				       true,false)) continue;
      
      if (Ereco1.size() != nLayers || Ereco2.size() != nLayers){
	std::cout << " Error! Not all layers are found... Fix code..." << std::endl;
	return 1;
      }

      FitResult recposfit1;
      //take layers around shower max
      unsigned minid = 0;
      unsigned maxid = nLayers-1;
      for (unsigned iL(0); iL<posx1.size();++iL){
	if (layerId1[iL]==10) minid = iL;
	if (layerId1[iL]==20) maxid = iL;
      }
      recposfit1.tanangle_x = (posx1[maxid]-posx1[minid])/(posz1[maxid]-posz1[minid]);
      recposfit1.tanangle_y = (posy1[maxid]-posy1[minid])/(posz1[maxid]-posz1[minid]);
      recposfit1.pos_x = posx1[0]-recposfit1.tanangle_x*posz1[0];
      recposfit1.pos_y = posy1[0]-recposfit1.tanangle_y*posz1[0];
      recposfit1.found=true;
      std::vector<ROOT::Math::XYZPoint> eventPos;
      eventPos.resize(nLayers,ROOT::Math::XYZPoint(0,0,0));
      for (unsigned iL(0); iL<posx1.size();++iL){
	eventPos[layerId1[iL]] = ROOT::Math::XYZPoint(posx1[iL],posy1[iL],posz1[iL]);
      }

      Signal1nofit.fillEnergies(ievt,(*ssvec),(*simhitvec),(*rechitvec),nPuVtx,eventPos);
      posFF1nofit = Signal1nofit.getAccuratePos(recposfit1,0);
      recoDir1nofit = Direction(recposfit1.tanangle_x,recposfit1.tanangle_y);

      FitResult recposfit2;
      minid = 0;
      maxid = nLayers-1;
      for (unsigned iL(0); iL<posx2.size();++iL){
	if (layerId2[iL]==10) minid = iL;
	if (layerId2[iL]==20) maxid = iL;
      }
      recposfit2.tanangle_x = (posx2[maxid]-posx2[minid])/(posz2[maxid]-posz2[minid]);
      recposfit2.tanangle_y = (posy2[maxid]-posy2[minid])/(posz2[maxid]-posz2[minid]);
      recposfit2.pos_x = posx2[0]-recposfit2.tanangle_x*posz2[0];
      recposfit2.pos_y = posy2[0]-recposfit2.tanangle_y*posz2[0];
      recposfit2.found=true;
      eventPos.resize(nLayers,ROOT::Math::XYZPoint(0,0,0));
      for (unsigned iL(0); iL<posx2.size();++iL){
	eventPos[layerId2[iL]] = ROOT::Math::XYZPoint(posx2[iL],posy2[iL],posz2[iL]);
      }

      Signal2nofit.fillEnergies(ievt,(*ssvec),(*simhitvec),(*rechitvec),nPuVtx,eventPos);
      posFF2nofit = Signal2nofit.getAccuratePos(recposfit2,0);
      recoDir2nofit = Direction(recposfit2.tanangle_x,recposfit2.tanangle_y);

      if (doFit1){
	//CAMM - Fix with using recoinfo...
	double eta1 = truthDir1.dir().Eta();
	double e1 = getCalibratedE(Ereco1,eta1);

	std::ostringstream mFolder;
	mFolder << singleGammaPath << getMatrixFolder(e1,eta1) << "_pu0";// << nVtx;
	
	//std::cout << " -- Ereco = " << e1 << " - matrix folders: " << mFolder.str() << std::endl;
	
	lGamma1.setMatrixFolder(mFolder.str());
	if (!lGamma1.fillMatrixFromFile(doOld)) {
	  nMatrixNotFound1++;
	  continue;
	}
	if ( lGamma1.performLeastSquareFit(ievt,fit1)==0){
	  found1 = Signal1.fillEnergies(ievt,(*ssvec),(*simhitvec),(*rechitvec),nPuVtx,fit1);
	} else std::cout << " -- Fit failed for photon 1." << std::endl;
      }
      else {
	found1 = Signal1.fillEnergies(ievt,(*ssvec),(*simhitvec),(*rechitvec),nPuVtx);
	fit1 = Signal1.getAccurateFit(ievt);
      }

      if (doFit2){
	//fix with using reco info
	double eta2 = truthDir2.dir().Eta();
	double e2 = getCalibratedE(Ereco2,eta2);

	std::ostringstream mFolder;
	mFolder << singleGammaPath << getMatrixFolder(e2,eta2) << "_pu0";// << nVtx;
	
	//std::cout << " -- Ereco = " << e2 << " - matrix folders: " << mFolder.str() << std::endl;
	
	lGamma2.setMatrixFolder(mFolder.str());
	if (!lGamma2.fillMatrixFromFile(doOld)) {
	  nMatrixNotFound2++;
	  continue;
	}
	if ( lGamma2.performLeastSquareFit(ievt,fit2)==0){
	  found2 = Signal2.fillEnergies(ievt,(*ssvec),(*simhitvec),(*rechitvec),nPuVtx,fit2);
	} else std::cout << " -- Fit failed for photon 2." << std::endl;
      }
      else {
	found2 = Signal2.fillEnergies(ievt,(*ssvec),(*simhitvec),(*rechitvec),nPuVtx);
	fit2 = Signal2.getAccurateFit(ievt);
      }
    }//if do fits
    else {
	found1 = Signal1.fillEnergies(ievt,(*ssvec),(*simhitvec),(*rechitvec),nPuVtx);
	found2 = Signal2.fillEnergies(ievt,(*ssvec),(*simhitvec),(*rechitvec),nPuVtx);
	fit1 = Signal1.getAccurateFit(ievt);
	fit2 = Signal2.getAccurateFit(ievt);
    }

    if (!found1 || !found2) continue;

    //get Higgs mass
    ROOT::Math::XYZPoint posFF1 = Signal1.getAccuratePos(fit1,0);
    recoDir1 = Direction(fit1.tanangle_x,fit1.tanangle_y);

    ROOT::Math::XYZPoint posFF2 = Signal2.getAccuratePos(fit2,0);
    recoDir2 = Direction(fit2.tanangle_x,fit2.tanangle_y);

    if (debug) {
      std::cout << " - Photon 1 direction:" << std::endl;
      std::cout << " Truth= "; truthDir1.Print();
      std::cout << " Reco = "; recoDir1.Print();
      std::cout << " - Photon 2 direction:" << std::endl;
      std::cout << " Truth= "; truthDir2.Print();
      std::cout << " Reco = "; recoDir2.Print();
    }


    TLorentzVector l1;
    double E1 = calibratedE(Signal1.getEtotalSR(2,true),recoDir1.dir().Eta());
    //additional smearing from constant term
    //E1 += lRndm.Gaus(0,0.01*E1);
    //std::cout << " Photon1: TrueE = " << truthE1 << " recoE = " << Signal1.getEtotalSR(2,true) << " calibE = " << E1 << std::endl;
    l1.SetPtEtaPhiE(recoDir1.dir().Pt()*E1,recoDir1.dir().Eta(),recoDir1.dir().Phi(),E1);
    TLorentzVector l2;
    double E2 = calibratedE(Signal2.getEtotalSR(2,true),recoDir2.dir().Eta());
    //E2 += lRndm.Gaus(0,0.01*E2);
    //std::cout << " Photon 2: TrueE = " << truthE2 << " recoE = " << Signal2.getEtotalSR(2,true)<< " calibE = " << E2 << std::endl;
    l2.SetPtEtaPhiE(recoDir2.dir().Pt()*E2,recoDir2.dir().Eta(),recoDir2.dir().Phi(),E2);


    //cut outside gen acceptance
    if (truthDir1.dir().Eta()<1.6 || truthDir1.dir().Eta()>2.7 || (truthE1/cosh(truthDir1.dir().Eta()))<=20 ||
	truthDir2.dir().Eta()<1.6 || truthDir2.dir().Eta()>2.7 || (truthE2/cosh(truthDir2.dir().Eta()))<=20) continue;

    //if (recoDir1.dir().Eta()<1.5 || recoDir1.dir().Eta()>2.8 || 
    //	recoDir2.dir().Eta()<1.5 || recoDir2.dir().Eta()>2.8 ||
    //	((E1/cosh(recoDir1.dir().Eta()))<40) || 
    //	((E2/cosh(recoDir2.dir().Eta())<40)) ) continue;

    nTwoPhotons++;

    hM.setRecoInfo(l1,l2,posFF1,posFF2);

    TLorentzVector l1nofit;
    double E1nofit = calibratedE(Signal1nofit.getEtotalSR(2,true),recoDir1nofit.dir().Eta());
    l1nofit.SetPtEtaPhiE(recoDir1nofit.dir().Pt()*E1nofit,recoDir1nofit.dir().Eta(),recoDir1nofit.dir().Phi(),E1nofit);
    TLorentzVector l2nofit;
    double E2nofit = calibratedE(Signal2nofit.getEtotalSR(2,true),recoDir2nofit.dir().Eta());
    l2nofit.SetPtEtaPhiE(recoDir2nofit.dir().Pt()*E2nofit,recoDir2nofit.dir().Eta(),recoDir2nofit.dir().Phi(),E2nofit);

    hMnofit.setRecoInfo(l1nofit,l2nofit,posFF1nofit,posFF2nofit);

    std::cout << " Check1: " << E1nofit << " " << truthE1 
	      << " (Enofit-Efit)/Efit " << (Signal1nofit.getEtotalSR(2,true)- Signal1.getEtotalSR(2,true))/Signal1.getEtotalSR(2,true)
	      << std::endl;
    std::cout << " Check2: " << E2nofit << " " << truthE2 
	      << " (Enofit-Efit)/Efit " << (Signal2nofit.getEtotalSR(2,true) - Signal2.getEtotalSR(2,true))/Signal2.getEtotalSR(2,true)
	      << std::endl;

    TLorentzVector t1;
    t1.SetPtEtaPhiE(truthDir1.dir().Pt()*truthE1,truthDir1.dir().Eta(),truthDir1.dir().Phi(),truthE1);
    TLorentzVector t2;
    t2.SetPtEtaPhiE(truthDir2.dir().Pt()*truthE2,truthDir2.dir().Eta(),truthDir2.dir().Phi(),truthE2);

    hM.setTruthInfo(t1,t2,
		    truthVtx1,truthVtx2);

    hMnofit.setTruthInfo(t1,t2,
			 truthVtx1,truthVtx2);

    hM.fillHistograms();
    hMnofit.fillHistograms();


    /*
    unsigned found = 0;
    for (unsigned iP(0); iP<(*genvec).size(); ++iP){//loop on gen particles    
      //if ((*genvec).size()!= 1) (*genvec)[iP].Print(std::cout);
      const HGCSSGenParticle lgen = (*genvec)[iP];
      if (lgen.trackID()<3){
	found++;
	//if (lgen.trackID()==1) g1 = ROOT::Math::PtEtaPhiEVector(lgen.pt(),lgen.eta(),lgen.phi(),lgen.E());
	//if (lgen.trackID()==2) g2 = ROOT::Math::PtEtaPhiEVector(lgen.pt(),lgen.eta(),lgen.phi(),lgen.E());
	//lgen.Print(std::cout);
      }
      if (found==2) break;
    }//loop on gen particles

    //std::cout << " -- Number of genparticles found: " << found << std::endl;
    if (found==2) {
      //std::cout << " -- Higgs mass = " << (g1+g2).M() << std::endl;
      nTwoPhotons++;
    }
    */



  }//loop on entries

  //finalise

  if (doFit1) lGamma1.finaliseFit();
  if (doFit2) lGamma2.finaliseFit();
  Signal1.finalise();
  Signal2.finalise();
  Signal1nofit.finalise();
  Signal2nofit.finalise();

  outputFile->Write();
  //outputFile->Close();
  
  std::cout << " -- Number of two photon events found: " << nTwoPhotons << std::endl;
  std::cout << " -- Number of matrix not found for photon 1: " << nMatrixNotFound1 << std::endl;
  std::cout << " -- Number of matrix not found for photon 2: " << nMatrixNotFound2 << std::endl;
  std::cout << " -- Number of events with closest cluster 1 away from truth : " << nTooFar1 << std::endl;
  std::cout << " -- Number of events with closest cluster 2 away from truth : " << nTooFar2 << std::endl;
  std::cout << " - End of higgsResolution program." << std::endl;

  return 0;
  

}//main
