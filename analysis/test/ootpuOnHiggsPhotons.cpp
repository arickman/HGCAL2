#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include <boost/algorithm/string.hpp>
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "boost/format.hpp"
#include "boost/function.hpp"

// helpful tools
#include "KDTreeLinkerAlgoT.h"
#include <unordered_map>
#include <unordered_set>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TProfile.h"
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

double absWeight(const unsigned layer){
  if (layer == 0) return 0.0378011;
  if (layer == 1) return 1;
  if (layer == 2) return 0.646989;
  if (layer == 3) return 0.617619;
  if (layer == 4) return 0.646989;
  if (layer == 5) return 0.617619;
  if (layer == 6) return 0.646989;
  if (layer == 7) return 0.617619;
  if (layer == 8) return 0.646989;
  if (layer == 9) return 0.617619;
  if (layer == 10) return 0.646989;
  if (layer == 11) return 0.942829;
  if (layer == 12) return 0.859702;
  if (layer == 13) return 0.942829;
  if (layer == 14) return 0.859702;
  if (layer == 15) return 0.942829;
  if (layer == 16) return 0.859702;
  if (layer == 17) return 0.942829;
  if (layer == 18) return 0.859702;
  if (layer == 19) return 0.942829;
  if (layer == 20) return 0.859702;
  if (layer == 21) return 1.37644;
  if (layer == 22) return 1.30447;
  if (layer == 23) return 1.37644;
  if (layer == 24) return 1.30447;
  if (layer == 25) return 1.37644;
  if (layer == 26) return 1.30447;
  if (layer == 27) return 1.37644;
  if (layer == 28) return 1.30447;
  if (layer == 29) return 1.79662;
  return 1;
};

double getSlope(const unsigned region){
  //return 73.4;
  if (region==3) return 80.4;//79.6;
  else if (region==2) return 81.1;
  else return 75.3;
};

double getOffset(const unsigned region){
  //return 1219;
  if (region==3) return 757;
  else if (region==2) return 727;
  else return 660;
};

double getEtotal(const std::vector<double> & Evec,
		 //const std::vector<unsigned> & region,
		 const double & eta){
  double Etot = 0;
  unsigned nL = Evec.size();
  for (unsigned iL(0); iL<nL;++iL){
    Etot += Evec[iL]*absWeight(iL);///getSlope(region[iL]);
  }
  double etacor = fabs(1./tanh(eta)); 
  return Etot*etacor;
}

double getCalibratedE(const std::vector<double> & Evec,
		      const std::vector<unsigned> & region,
		      const double eta){
  double Etot = getEtotal(Evec,eta);//region);
  double offset = getOffset(region[0]);//-77;//   +/-   86.7909     
  double slope = getSlope(region[0]);//79.6;//   +/-   0.455393
  return (Etot-offset)/slope;
  //calibration for signal region 2: 3*3 cm^2
  //return calibratedE(Etot,eta);
};

double radius(const double & aEta, const double & aZ){
  double theta = 2*atan(exp(-1.*aEta));
  return aZ*tan(theta);
};

double areaInCm2(const double & aR){
  return TMath::Pi()*pow(aR/10.,2);
};

bool getZpositions(const unsigned versionNumber, std::vector<double> & posz){
  std::ifstream fin;
  std::ostringstream finname;
  //finname << outFolder_ << "/zPositions.dat";
  finname << "data/zPositions_v" << versionNumber << ".dat";
  fin.open(finname.str());
  if (!fin.is_open()){
    std::cout << " Cannot open input file " << finname.str() << "! Refilling now..." << std::endl;
    return false;
  }
  
  std::cout << " Reading z position per layer from input file " << finname.str() << std::endl;
  
  while (!fin.eof()){
    unsigned l=posz.size();
    double z=0;
    fin>>l>>z;
    if (l<posz.size()){
      posz[l]=z;
      std::cout << " Layer " << l << ", z = " << z << std::endl;
    }
  }
  
  fin.close();
  return true;
  
}

void getMaximumCellFromGeom(const double & phimax,const double & etamax,std::vector<double> & xmax,std::vector<double> & ymax, const std::vector<double> & posz){

  for (unsigned iL(0); iL<xmax.size();++iL){
    double theta = 2*atan(exp(-1.*etamax));
    double rho = posz[iL]/cos(theta);
    xmax[iL] = rho*sin(theta)*cos(phimax);
    ymax[iL] = rho*sin(theta)*sin(phimax);

    if (xmax[iL]>0) xmax[iL]=static_cast<int>((xmax[iL]+4.999999)/10.)*10;
    else xmax[iL]=static_cast<int>((xmax[iL]-4.999999)/10.)*10;
    if (ymax[iL]>0) ymax[iL]=static_cast<int>((ymax[iL]+4.999999)/10.)*10;
    else ymax[iL]=static_cast<int>((ymax[iL]-4.999999)/10.)*10;

  }//loop on layers

}

struct MyHit{
  double x;
  double y;
  double z;
  double e;
  unsigned layer;
  bool hasSignal;
  ROOT::Math::XYZPoint position() const{
    return ROOT::Math::XYZPoint(x/10.,y/10.,z/10.);//in cm
  };
  double E() const {
    return e;
  };
};

void fillKDTree(KDTree & _hit_kdtree,
		std::vector<KDNode> & _hit_nodes,
		const std::vector<MyHit> & rechits){


  _hit_kdtree.clear();
  std::vector<bool> usable_rechits(rechits.size(),true);

  _hit_nodes.reserve(rechits.size());

  KDTreeCube kd_boundingregion = 
    fill_and_bound_kd_tree(rechits,usable_rechits,_hit_nodes);
  _hit_kdtree.build(_hit_nodes,kd_boundingregion);
  _hit_nodes.clear();
};

void findNeighbours(const double & cell_size_mm,
		    KDTree & _hit_kdtree,
		    const MyHit & current_cell, 
		    std::vector<KDNode> & found, 
		    const unsigned layerRange=1){
  const ROOT::Math::XYZPoint pos = current_cell.position();


  double cell_size = cell_size_mm/10.+0.01;//cm +0.01safety margin for double precision...
  unsigned layer = current_cell.layer;

  //std::cout << " -- x,y,z,l = " << pos.X() << " " << pos.Y() << " " << pos.Z() << " " << layer << std::endl;
  //std::cout << " -- kdtree size: " << _hit_kdtree.size() << std::endl;

  auto x_rh = minmax(pos.x()+cell_size,pos.x()-cell_size);
  auto y_rh = minmax(pos.y()+cell_size,pos.y()-cell_size);
  auto z_rh = minmax(zPos(layer-layerRange)/10.,zPos(layer+layerRange)/10.);//in cm!!!

  KDTreeCube hit_searchcube((float)x_rh.first,(float)x_rh.second,
			    (float)y_rh.first,(float)y_rh.second,
			    (float)z_rh.first,(float)z_rh.second);

  //std::cout << " -- min and max X = " << (float)x_rh.first << " " << (float)x_rh.second << std::endl
  //<< " -- min and max Y = " << (float)y_rh.first << " " <<(float)y_rh.second  << std::endl
  //<< " -- min and max Z = " << (float)z_rh.first << " " << (float)z_rh.second << std::endl;


  _hit_kdtree.search(hit_searchcube,found);

  //std::cout << " -- Found " << found.size() << " nodes." << std::endl;

};

double getFactor(const double & radius){
  if (radius<750) return 2;
  else if (radius > 1200) return 2./3;
  return 1;
};

unsigned getRegion(const double & radius){
  if (radius<750) return 1;
  else if (radius > 1200) return 3;
  return 2;
};

double getCellSize(const double & radius){
  if (radius<750) return 7.5;
  return 10;
};

int main(int argc, char** argv){//main  

  //Input output and config options
  std::string cfg;
  bool concept;
  unsigned pNevts;
  std::string inFilePath;
  std::string outFilePath;
  unsigned debug;
  unsigned nPu;
  unsigned run;

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
    ("debug,d",        po::value<unsigned>(&debug)->default_value(0))
    ("nPu,p",          po::value<unsigned>(&nPu)->default_value(140))
    ("run,r",          po::value<unsigned>(&run)->default_value(0))
    ;

  po::store(po::command_line_parser(argc, argv).options(config).allow_unregistered().run(), vm);
  po::store(po::parse_config_file<char>(cfg.c_str(), config), vm);
  po::notify(vm);


  std::cout << " -- Input parameters: " << std::endl
	    << " -- Input file path: " << inFilePath << std::endl
	    << " -- Output file path: " << outFilePath << std::endl
	    << " -- Requiring " << nPu << " pu events." << std::endl
	    << " -- Processing ";
  if (pNevts == 0) std::cout << "all events." << std::endl;
  else std::cout << pNevts << " events." << std::endl;


  /////////////////////////////////////////////////////////////
  //hardcoded
  /////////////////////////////////////////////////////////////

  const double mipE = 2.35;//fC/mip
  //const double cell_size = 10;//mm

  bool doHgg = true;

  /////////////////////////////////////////////////////////////
  //input
  /////////////////////////////////////////////////////////////

  //DigiPu140_version12_model2_BOFF_et70_alpha0.164_phi0.250pi_run0.root
  std::ostringstream inputbase;
  inputbase << inFilePath; //<< "/MipStudy.root";
  TChain *lTree = new TChain("RecoTree");

  for (unsigned i(1); i<11;++i){
    std::ostringstream input;
    input << inputbase.str() << "_" << i
	  << ".root";
    lTree->AddFile(input.str().c_str());
    std::cout << "Adding file:" << input.str().c_str() << std::endl;
  }
  
  //TFile *simFile = 0;
  //if (!testInputFile(input.str(),simFile)) return 1;
  //TTree *lTree = (TTree*)simFile->Get("RecoTree");
  if (!lTree){
    std::cout << " -- Error, tree RecoTree cannot be opened. Exiting..." << std::endl;
    return 1;
  }

  /////////////////////////////////////////////////////////////
  //Info
  /////////////////////////////////////////////////////////////

  //HGCSSInfo * info=(HGCSSInfo*)simFile->Get("Info");
  const double cellSize = 2.5;
  const unsigned versionNumber = 12;//info->version();
  const unsigned model = 2;//info->model();
  
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
 
  myDetector.buildDetector(versionNumber,concept,isCaliceHcal);

  HGCSSCalibration mycalib(inFilePath);

  const unsigned nLayers = myDetector.nLayers();
  const unsigned nSections = myDetector.nSections();

  std::cout << " -- N layers = " << nLayers << std::endl
	    << " -- N sections = " << nSections << std::endl;


  TRandom3 lRndm(0);
  std::cout << " -- Random number seed: " << lRndm.GetSeed() << std::endl;

  TFile *inputSig = 0;
  std::ostringstream input;
  input << "/afs/cern.ch/work/a/amagnan/PFCalEEAna/gitV00-02-14/version12/Hgg/pu0/run" << run <<".root";
  if (doHgg) inputSig = TFile::Open(input.str().c_str());
  else inputSig = TFile::Open("PLOTS/gitV00-02-12/version12/gamma/200um/eta19_et60_pu0_new_logw3.root");
  if (!inputSig) {
    std::cout << " Signal file not found." << std::endl;
    return 1;
  }
  if (doHgg) inputSig->cd("Gamma1Fit");
  TTree *sigTree1 = (TTree*)gDirectory->Get("EcellsSR2");
  std::vector<double> init;
  init.resize(25,0);
  std::vector<double> truthPosXg1;
  std::vector<double> truthPosYg1;
  std::vector<std::vector<double> > Exyg1;
  double Etruth1 = 0;
  if (doHgg) sigTree1->SetBranchAddress("truthE",&Etruth1);
  else Etruth1 = 60.*cosh(1.9);
  Exyg1.resize(nLayers,init);
  truthPosXg1.resize(nLayers,0);
  truthPosYg1.resize(nLayers,0);
  for (unsigned iL(0);iL<nLayers;++iL){
    std::ostringstream label;
    label.str("");     
    label << "TruthPosX_" << iL;
    sigTree1->SetBranchAddress(label.str().c_str(),&truthPosXg1[iL]);
    label.str("");     
    label << "TruthPosY_" << iL;
    sigTree1->SetBranchAddress(label.str().c_str(),&truthPosYg1[iL]);
    for (unsigned iy(0);iy<3;++iy){
      for (unsigned ix(0);ix<3;++ix){
	unsigned idx = 3*iy+ix;
	label.str("");     
	label << "E_" << iL << "_" << idx;
	sigTree1->SetBranchAddress(label.str().c_str(),&Exyg1[iL][idx]);
      }
    }
  }

  if (doHgg) inputSig->cd("Gamma2Fit");
  TTree *sigTree2 = 0;
  if (doHgg) sigTree2 = (TTree*)gDirectory->Get("EcellsSR2");

  double Etruth2 = 0;
  if (doHgg) sigTree2->SetBranchAddress("truthE",&Etruth2);
  std::vector<double> truthPosXg2;
  std::vector<double> truthPosYg2;
  std::vector<std::vector<double> > Exyg2;
  Exyg2.resize(nLayers,init);
  truthPosXg2.resize(nLayers,0);
  truthPosYg2.resize(nLayers,0);
  if (doHgg) {
    for (unsigned iL(0);iL<nLayers;++iL){
      std::ostringstream label;
      label.str("");     
      label << "TruthPosX_" << iL;
      sigTree2->SetBranchAddress(label.str().c_str(),&truthPosXg2[iL]);
      label.str("");     
      label << "TruthPosY_" << iL;
      sigTree2->SetBranchAddress(label.str().c_str(),&truthPosYg2[iL]);
      for (unsigned iy(0);iy<3;++iy){
	for (unsigned ix(0);ix<3;++ix){
	  unsigned idx = 3*iy+ix;
	  label.str("");     
	  label << "E_" << iL << "_" << idx;
	  sigTree2->SetBranchAddress(label.str().c_str(),&Exyg2[iL][idx]);
	}
      }
    }
  }

  //////////////////////////////////////////////////
  //////////////////////////////////////////////////
  ///////// Output File // /////////////////////////
  //////////////////////////////////////////////////
  //////////////////////////////////////////////////


  std::ostringstream output;
  output << outFilePath << "_pu" << nPu << "_run" << run << ".root";
  TFile *outputFile = TFile::Open(output.str().c_str(),"RECREATE");
  
  if (!outputFile) {
    std::cout << " -- Error, output file " << output.str().c_str() << " cannot be opened. Please create output directory. Exiting..." << std::endl;
    return 1;
  }
  else {
    std::cout << " -- output file " << outputFile->GetName() << " successfully opened." << std::endl;
  }
  outputFile->cd();

  TH1F *hphi = new TH1F("hphi",";#phi;showers",100,-3.1416,3.1416);
  TH1F *heta = new TH1F("heta",";#eta;showers",100,1.5,3.0);
  TH1F *hpt = new TH1F("hpt",";p_{T} (GeV);showers",180,20,200);
  TH1F *hE = new TH1F("hE",";E (GeV);showers",100,0,1000);

  TH1F *hphi1 = new TH1F("hphi1",";#phi;showers",100,-3.1416,3.1416);
  TH1F *heta1 = new TH1F("heta1",";#eta;showers",100,1.5,3.0);
  TH1F *hpt1 = new TH1F("hpt1",";p_{T} (GeV);showers",180,20,200);
  TH1F *hE1 = new TH1F("hE1",";E (GeV);showers",100,0,1000);

  TH1F *hphi2 = new TH1F("hphi2",";#phi;showers",100,-3.1416,3.1416);
  TH1F *heta2 = new TH1F("heta2",";#eta;showers",100,1.5,3.0);
  TH1F *hpt2 = new TH1F("hpt2",";p_{T} (GeV);showers",180,20,200);
  TH1F *hE2 = new TH1F("hE2",";E (GeV);showers",100,0,1000);

  TH1F *E_noOOTPU = new TH1F("E_noOOT",";E_{#gamma}-E_{truth} (GeV);showers",1000,-100,100);
  TH1F *E_mips = new TH1F("E_mips",";E_{#gamma} (mips);showers",10000,0,100000);
  TH1F *E_OOTPU = new TH1F("E_OOT",";E_{#gamma}-E_{truth} (GeV);showers",1000,-100,100);
  TH1F *E_OOTpastPU = new TH1F("E_OOTpast",";E_{#gamma}-E_{truth} (GeV);showers",1000,-100,100);
  TH2F *EcalibvsTruth = new TH2F("EcalibvsTruth",";E_{truth} (GeV);E_{#gamma} (GeV);showers",5000,0,5000,5000,0,5000);
  TH2F *ErawvsTruth = new TH2F("ErawvsTruth",";E_{truth} (GeV);E_{#gamma} (Mips);showers",1000,0,1000,10000,0,50000);
  TH2F *ErawvsTruth1 = new TH2F("ErawvsTruth1",";E_{truth} (GeV);E_{#gamma} (Mips);showers",1000,0,1000,10000,0,50000);
  TH2F *ErawvsTruth2 = new TH2F("ErawvsTruth2",";E_{truth} (GeV);E_{#gamma} (Mips);showers",1000,0,1000,10000,0,50000);
  TH2F *ErawvsTruth3 = new TH2F("ErawvsTruth3",";E_{truth} (GeV);E_{#gamma} (Mips);showers",1000,0,1000,10000,0,50000);
  TH2F *EratiovsEta = new TH2F("EratiovsEta",";#eta;E_{#gamma}/E_{truth};showers",100,1.5,3.0,400,0,1.5);



  TH1F *ErecooverEtrue = new TH1F("ErecooverEtrue",";E_{#gamma}/E_{truth};showers",400,0,1.5);
  TH1F *ErecooverEtrueOOT = new TH1F("ErecooverEtrueOOT",";E_{#gamma}/E_{truth};showers",400,0,1.5);
  TH1F *ErecooverEtrueOOTpast = new TH1F("ErecooverEtrueOOTpast",";E_{#gamma}/E_{truth};showers",400,0,1.5);

  TH1F *EpuoverErecoFutur = new TH1F("EpuoverErecoFutur",";E_{PU}/E_{#gamma};showers",1000,0,0.1);
  TH1F *EpuoverErecoPast = new TH1F("EpuoverErecoPast",";E_{PU}/E_{#gamma};showers",1000,0,0.1);
  TH1F *EpuoverErecoAll = new TH1F("EpuoverErecoAll",";E_{PU}/E_{#gamma};showers",1000,0,0.1);

  TH1F *EpuFutur = new TH1F("EpuFutur",";E_{PU} (mips);showers",1000,0,500);
  TH1F *EpuPast = new TH1F("EpuPast",";E_{PU} (mips);showers",1000,0,500);
  TH1F *EpuAll = new TH1F("EpuAll",";E_{PU} (mips);showers",1000,0,500);

  const unsigned nbx = 11;
  TH1F *Epu_perbx[nbx];
  TH1F *nAbove_perbx[nbx];
  TH2F *nAbovevsEta = new TH2F("nAbovevsEta",";#eta;n_{cells}(E>100 fC);showers",
			       14,1.5,2.9,
			       200,0,200);

  double bxthresh[nbx] = {250,250,1000,1930,2910,3870,4870,5880,6930,7990,9080};

  for (unsigned ibx(0); ibx<nbx;++ibx){
    bxthresh[ibx] = bxthresh[ibx]/mipE;
    std::ostringstream label;
    label.str("");
    label << "Epu_bx_" << ibx;
    Epu_perbx[ibx] = new TH1F(label.str().c_str(),";E_{PU} (mips);showers",500,0,500);
    label.str("");
    label << "nAbove_bx_" << ibx;
    nAbove_perbx[ibx] = new TH1F(label.str().c_str(),";n(E>thresh);showers",271,0,271);
  }


  const unsigned nEvts = ((pNevts > sigTree1->GetEntries() || pNevts==0) ? static_cast<unsigned>(sigTree1->GetEntries()) : pNevts) ;

  std::vector<HGCSSRecoHit> * rechitvec = 0;
  unsigned nPuVtx = 0;
  lTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
  if (lTree->GetBranch("nPuVtx")) lTree->SetBranchAddress("nPuVtx",&nPuVtx);

  std::vector<double> posz;
  posz.resize(nLayers,0);
  if (!getZpositions(versionNumber,posz)) return 1;

  //loop on events
  const unsigned nPuEvts = lTree->GetEntries();
  std::cout << "- Number of pu events available: " << lTree->GetEntries()  << std::endl;
  std::cout << "- Number of signal events available: " << sigTree1->GetEntries();
  if (doHgg) std::cout << " " << sigTree2->GetEntries();
  std::cout << std::endl;

  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
    if (debug) std::cout << "... Processing entry: " << ievt << std::endl;
    else if (ievt%10 == 0) std::cout << "... Processing entry: " << ievt << std::endl;

    sigTree1->GetEntry(ievt);
    if (doHgg) sigTree2->GetEntry(ievt);

    std::vector<double> EperLayer1;
    EperLayer1.resize(nLayers,0);
    std::vector<double> EperLayer1OOT;
    EperLayer1OOT.resize(nLayers,0);
    std::vector<double> EperLayer2;
    EperLayer2.resize(nLayers,0);
    std::vector<double> EperLayer2OOT;
    EperLayer2OOT.resize(nLayers,0);

    std::set<unsigned> xyid1[nbx][nLayers];
    std::set<unsigned> xyid2[nbx][nLayers];
    std::vector<unsigned> nCellsAbove1;
    nCellsAbove1.resize(nbx,0);
    std::vector<unsigned> nCellsAbove2;
    nCellsAbove2.resize(nbx,0);
    unsigned nover100_1 = 0;
    unsigned nover100_2 = 0;
    Direction dir1(truthPosXg1[0]/posz[0],truthPosYg1[0]/posz[0]);
    Direction dir2(truthPosXg2[0]/posz[0],truthPosYg2[0]/posz[0]);

    std::vector<double> xmax1;
    xmax1.resize(nLayers,0);
    std::vector<double> ymax1;
    ymax1.resize(nLayers,0);
    double phi1 = dir1.phi();
    double eta1 = dir1.eta();
    double theta1 = 2*atan(exp(-1.*eta1));
    hphi1->Fill(phi1);
    heta1->Fill(eta1);
    hpt1->Fill(Etruth1/cosh(eta1));
    hE1->Fill(Etruth1);

    if (debug) std::cout << " Truth photon: x=" << truthPosXg1[0] << " y=" << truthPosYg1[0] << " E=" << Etruth1 << " eta=" << eta1 << " phi=" << phi1 << std::endl;

    getMaximumCellFromGeom(phi1,eta1,xmax1,ymax1,posz);

    std::vector<double> xmax2;
    xmax2.resize(nLayers,0);
    std::vector<double> ymax2;
    ymax2.resize(nLayers,0);
    double phi2 = dir2.phi();
    double eta2 = dir2.eta();
    double theta2 = 2*atan(exp(-1.*eta2));
    getMaximumCellFromGeom(phi2,eta2,xmax2,ymax2,posz);
    hphi2->Fill(phi2);
    heta2->Fill(eta2);
    hpt2->Fill(Etruth2/cosh(eta2));
    hE2->Fill(Etruth2);

    bool fid1 = eta1>1.5 && eta1<2.9;
    fid1 = fid1 && (Etruth1/cosh(eta1))>40;
    bool fid2 = eta2>1.5 && eta2<2.9;
    fid2 = fid2 && (Etruth2/cosh(eta2))>40;

    if (!fid1 && !fid2) continue;

    bool region1[3] = {false,false,false};
    bool region2[3] = {false,false,false};
    std::vector<unsigned> regiong1;
    regiong1.resize(nLayers,0);
    std::vector<unsigned> regiong2;
    regiong2.resize(nLayers,0);

    for (unsigned iL(0);iL<nLayers;++iL){
      double r1 = sqrt(pow(xmax1[iL],2)+pow(ymax1[iL],2));
      double f1 = getFactor(r1);
      double r2 = sqrt(pow(xmax2[iL],2)+pow(ymax2[iL],2));
      double f2 = getFactor(r2);
      region1[getRegion(r1)-1] = true;
      region2[getRegion(r2)-1] = true;
      regiong1[iL] = getRegion(r1);
      regiong2[iL] = getRegion(r2);
      //double cs1 = getCellSize(r1);
      //double cs2 = getCellSize(r2);

      for (unsigned iy(0);iy<3;++iy){
	for (unsigned ix(0);ix<3;++ix){
	  unsigned idx = 3*iy+ix;
	  if (debug>1) std::cout << iL << " " << idx << " " << Exyg1[iL][idx] << std::endl;

	  //if (fabs(ix-1)<1.5 && fabs(iy-2)<1.5) {
	  EperLayer1[iL] += Exyg1[iL][idx];
	  EperLayer2[iL] += Exyg2[iL][idx];
	    //}

	  if (Exyg1[iL][idx]>(100./mipE*f1)) nover100_1++;	    
	  if (Exyg2[iL][idx]>(100./mipE*f2)) nover100_2++;	    
	  //fill ids of cells above thresh
	  for (unsigned ibx(1); ibx<nbx;++ibx){//loop on bx
	    if (Exyg1[iL][idx]>(bxthresh[ibx]*f1)) {
	      xyid1[ibx][iL].insert(idx);
	      nCellsAbove1[ibx]++;
	    }
	    if (Exyg2[iL][idx]>(bxthresh[ibx]*f2)) {
	      xyid2[ibx][iL].insert(idx);
	      nCellsAbove2[ibx]++;
	    }
	  }

	}
      }
      EperLayer1OOT[iL] = EperLayer1[iL];
      EperLayer2OOT[iL] = EperLayer2[iL];
    }//loop on layers


    fid1 = fid1 && (
    		    (region1[0] && !region1[1] && !region1[2]) ||
		    (!region1[0] && region1[1] && !region1[2]) ||
    		    (!region1[0] && !region1[1] && region1[2])
		    );
    fid2 = fid2 && (
    		    (region2[0] && !region2[1] && !region2[2]) ||
		    (!region2[0] && region2[1] && !region2[2]) ||
    		    (!region2[0] && !region2[1] && region2[2])
		    );

    if (!fid1 && !fid2) continue;

    for (unsigned ibx(1); ibx<nbx;++ibx){//loop on bx
      if (fid1) nAbove_perbx[ibx]->Fill(nCellsAbove1[ibx]);
      if (fid2) nAbove_perbx[ibx]->Fill(nCellsAbove2[ibx]);
    }

    if (fid1){
      hphi->Fill(phi1);
      heta->Fill(eta1);
      hpt->Fill(Etruth1/cosh(eta1));
      hE->Fill(Etruth1);
      nAbovevsEta->Fill(eta1,nover100_1);
    }
    if (fid2){
      hphi->Fill(phi2);
      heta->Fill(eta2);
      hpt->Fill(Etruth2/cosh(eta2));
      hE->Fill(Etruth2);
      nAbovevsEta->Fill(eta2,nover100_2);
    }

    std::set<unsigned> lidxSet;

    double EtotPu1 = 0;
    double EtotPuFutur1 = 0;
    double EtotPu2 = 0;
    double EtotPuFutur2 = 0;

    for (unsigned ibx(0); ibx<nbx;++ibx){//loop on bx

      if (ibx>0 && nCellsAbove1[ibx]==0 && nCellsAbove2[ibx]==0) continue; 

      unsigned ipu = 0;
      while (1){
	ipu = lRndm.Integer(nPuEvts);
	if (lidxSet.find(ipu)==lidxSet.end()){
	  lidxSet.insert(ipu);
	  break;
	}
	else {
	  std::cout << " -- Found duplicate ! Taking another shot." << std::endl;
	}
      }
      lTree->GetEntry(ipu);
      
      unsigned nover1 = 0;
      double Ebx1 = 0;
      unsigned nover2 = 0;
      double Ebx2 = 0;
      for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on hits
	const HGCSSRecoHit & lHit = (*rechitvec)[iH];
	
	double posx = lHit.get_x();
	double posy = lHit.get_y();
	double lR = sqrt(pow(posx,2)+pow(posy,2));
	unsigned layer = lHit.layer();
	//double leta = lHit.eta();
	double energy = lHit.energy();
	if (fabs(posx-xmax1[layer]) < 15.1 && 
	    fabs(posy-ymax1[layer]) < 15.1){
	  int ix = (posx-xmax1[layer])/getCellSize(lR);
	  int iy = (posy-ymax1[layer])/getCellSize(lR);
	  unsigned idx = 3*(iy+1)+(ix+1);
	  double r1 = sqrt(pow(xmax1[layer],2)+pow(ymax1[layer],2));
	  double f1 = getFactor(r1);
	  if ((ibx==0 && energy>(bxthresh[0]*f1)) || 
	      (ibx>0 && xyid1[ibx][layer].find(idx)!=xyid1[ibx][layer].end())){
	    EperLayer1OOT[layer] += energy;
	    nover1++;
	    Ebx1 += energy;
	  }

	}
	if (doHgg && fabs(posx-xmax2[layer]) < 15.1 && 
	    fabs(posy-ymax2[layer]) < 15.1){
	  int ix = (posx-xmax2[layer])/getCellSize(lR);
	  int iy = (posy-ymax2[layer])/getCellSize(lR);
	  unsigned idx = 3*(iy+1)+(ix+1);
	  double r2 = sqrt(pow(xmax2[layer],2)+pow(ymax2[layer],2));
	  double f2 = getFactor(r2);

	  if ((ibx==0 && energy>(bxthresh[0]*f2)) ||
	      (ibx>0 && xyid2[ibx][layer].find(idx)!=xyid2[ibx][layer].end())
	      ){
	    EperLayer2OOT[layer] += energy;
	    nover2++;
	    Ebx2 += energy;
	  }
	}
      }//loop on hits
      if (fid1) {
	Epu_perbx[ibx]->Fill(Ebx1);
      }
      EtotPu1 += Ebx1;
      if (doHgg){
	if (fid2){
	  Epu_perbx[ibx]->Fill(Ebx2);
	}
	EtotPu2 += Ebx2;
      }
      if (ibx==0){
	if (fid1){
	  E_OOTpastPU->Fill(getCalibratedE(EperLayer1OOT,regiong1,eta1)-Etruth1);
	  ErecooverEtrueOOTpast->Fill(getCalibratedE(EperLayer1OOT,regiong1,eta1)/Etruth1);
	  EpuoverErecoPast->Fill(Ebx1/getEtotal(EperLayer1,eta1));
	  EpuPast->Fill(Ebx1);
	}
	if (doHgg && fid2){
	  E_OOTpastPU->Fill(getCalibratedE(EperLayer2OOT,regiong2,eta2)-Etruth2);
	  ErecooverEtrueOOTpast->Fill(getCalibratedE(EperLayer2OOT,regiong2,eta2)/Etruth2);
	  EpuoverErecoPast->Fill(Ebx2/getEtotal(EperLayer2,eta2));
	  EpuPast->Fill(Ebx2);
	}
      } else {
	EtotPuFutur1 += Ebx1;
	if (doHgg) EtotPuFutur2 += Ebx2;
      }
    }//loop on bx
    if (fid1){
      EpuoverErecoFutur->Fill(EtotPuFutur1/getEtotal(EperLayer1,eta1));
      EpuoverErecoAll->Fill(EtotPu1/getEtotal(EperLayer1,eta1));
      EpuFutur->Fill(EtotPuFutur1);
      EpuAll->Fill(EtotPu1);
      E_noOOTPU->Fill(getCalibratedE(EperLayer1,regiong1,eta1)-Etruth1);
      E_mips->Fill(getEtotal(EperLayer1,eta1));
      EcalibvsTruth->Fill(Etruth1,getCalibratedE(EperLayer1,regiong1,eta1));
      ErawvsTruth->Fill(Etruth1,getEtotal(EperLayer1,eta1));
      if (region1[0] && !region1[1] && !region1[2]) ErawvsTruth1->Fill(Etruth1,getEtotal(EperLayer1,eta1));
      if (!region1[0] && region1[1] && !region1[2]) ErawvsTruth2->Fill(Etruth1,getEtotal(EperLayer1,eta1));
      if (!region1[0] && !region1[1] && region1[2]) ErawvsTruth3->Fill(Etruth1,getEtotal(EperLayer1,eta1));
      E_OOTPU->Fill(getCalibratedE(EperLayer1OOT,regiong1,eta1)-Etruth1);
      ErecooverEtrue->Fill(getCalibratedE(EperLayer1,regiong1,eta1)/Etruth1);
      EratiovsEta->Fill(eta1,getCalibratedE(EperLayer1,regiong1,eta1)/Etruth1);
      ErecooverEtrueOOT->Fill(getCalibratedE(EperLayer1OOT,regiong1,eta1)/Etruth1);
    }
    if (doHgg && fid2) {
      EpuoverErecoFutur->Fill(EtotPuFutur2/getEtotal(EperLayer2,eta2));
      EpuoverErecoAll->Fill(EtotPu2/getEtotal(EperLayer2,eta2));
      EpuFutur->Fill(EtotPuFutur2);
      EpuAll->Fill(EtotPu2);
      E_noOOTPU->Fill(getCalibratedE(EperLayer2,regiong2,eta2)-Etruth2);
      EcalibvsTruth->Fill(Etruth2,getCalibratedE(EperLayer2,regiong2,eta2));
      ErawvsTruth->Fill(Etruth2,getEtotal(EperLayer2,eta2));
      if (region2[0] && !region2[1] && !region2[2]) ErawvsTruth1->Fill(Etruth2,getEtotal(EperLayer2,eta2));
      if (!region2[0] && region2[1] && !region2[2]) ErawvsTruth2->Fill(Etruth2,getEtotal(EperLayer2,eta2));
      if (!region2[0] && !region2[1] && region2[2]) ErawvsTruth3->Fill(Etruth2,getEtotal(EperLayer2,eta2));
      EratiovsEta->Fill(eta2,getCalibratedE(EperLayer2,regiong2,eta2)/Etruth2);
      E_OOTPU->Fill(getCalibratedE(EperLayer2OOT,regiong2,eta2)-Etruth2);
      ErecooverEtrue->Fill(getCalibratedE(EperLayer2,regiong2,eta2)/Etruth2);
      ErecooverEtrueOOT->Fill(getCalibratedE(EperLayer2OOT,regiong2,eta2)/Etruth2);
    }
  }//loop on entries

  outputFile->cd();
  outputFile->Write();
  outputFile->Close();

  return 0;
  
}//main
