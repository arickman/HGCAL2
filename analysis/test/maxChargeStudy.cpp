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

int main(int argc, char** argv){//main  

  //Input output and config options
  std::string cfg;
  bool concept;
  unsigned pNevts;
  std::string inFilePath;
  std::string outFilePath;
  unsigned debug;
  unsigned nPu;

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

  const unsigned nEta = 1;

  const double eta[nEta] = {2.0};//1.7,2.0,2.5};

  const double deta = 0.1;
  const double Ethresh = 106;//106;

  const double cell_size = 10;//mm
  //const double Ethresh = 0.6;

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

  TFile *inputSig = TFile::Open("/afs/cern.ch/work/a/amagnan/PFCalEEAna/PLOTS/gitV00-02-14/version12/Hgg/200um/E3by3pu0.root");
  inputSig->cd("Gamma1Fit");
  TTree *sigTree1 = (TTree*)gDirectory->Get("EcellsSR2");
  inputSig->cd("Gamma2Fit");
  TTree *sigTree2 = (TTree*)gDirectory->Get("EcellsSR2");


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

  //////////////////////////////////////////////////
  //////////////////////////////////////////////////
  ///////// Output File // /////////////////////////
  //////////////////////////////////////////////////
  //////////////////////////////////////////////////


  std::ostringstream output;
  output << outFilePath << "_pu" << nPu << ".root";
  TFile *outputFile = TFile::Open(output.str().c_str(),"RECREATE");
  
  if (!outputFile) {
    std::cout << " -- Error, output file " << output.str().c_str() << " cannot be opened. Please create output directory. Exiting..." << std::endl;
    return 1;
  }
  else {
    std::cout << " -- output file " << outputFile->GetName() << " successfully opened." << std::endl;
  }
  outputFile->cd();

  TH2F *hxy = new TH2F("hxy",";x (mm);y (mm); hits",
			 113,-1695,1695,
			113,-1695,1695);

  TH2F *hphieta = new TH2F("hphieta",";#phi;#eta; hits",
			   100,-3.1416,3.1416,			     
			   50,1.6,2.4);
  TH2F *hxyAbove = new TH2F("hxyAbove",";x (mm);y (mm); hits (E>250 fC)",
			    113,-1695,1695,
			    113,-1695,1695);

  TH2F *hphietaAbove = new TH2F("hphietaAbove",";#phi;#eta; hits (E>250 fC)",
				100,-3.1416,3.1416,			     
				50,1.6,2.4);

  TH1F *nAbove = new TH1F("nAbove",";n(E>250 fC)",50,0,50);
  TH1F *nTot = new TH1F("nTot",";n in SR",300,0,300);
  TH2F *nAbovevsLayer = new TH2F("nAbovevsLayer",";layer;n_{1#times 1 cm^{2}}(E>250 fC)/n(1.9 < #eta < 2.1)",
				 nLayers,0,nLayers,10000,0,0.005);
  TProfile *prof_nAbovevsLayer = new TProfile("prof_nAbovevsLayer",";layer;n_{1#times 1 cm^{2}}(E>250 fC)/n(1.9 < #eta < 2.1)",
				 nLayers,0,nLayers,0,1);

  TH1F *EmaxNeighbourSC3D = new TH1F("EmaxNeighbourSC3D",";E_{max}^{SC3D} (Mips);hits in 1.9 < #eta < 2.1",400,0,100);

  TH2D *EneighvsEcell = new TH2D("EneighvsEcell",
				 ";E_{cell};E_{max}^{SC3D}",
				 100,100,1000,
				 400,0,100);

  TProfile *hLayerProfile_above50mips = new TProfile("hLayerProfile_above50mips",";layer;E (mips) in 1.9 < #eta < 2.1",nLayers,0,nLayers,0,10000);

  const unsigned nEvts = ((pNevts > lTree->GetEntries() || pNevts==0) ? static_cast<unsigned>(lTree->GetEntries()) : pNevts) ;

  std::vector<HGCSSRecoHit> * rechitvec = 0;
  unsigned nPuVtx = 0;
  lTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
  if (lTree->GetBranch("nPuVtx")) lTree->SetBranchAddress("nPuVtx",&nPuVtx);

  std::vector<double> posz;
  posz.resize(nLayers,0);
  if (!getZpositions(versionNumber,posz)) return 1;

  //loop on events
  const unsigned nPuEvts = lTree->GetEntries();
  std::cout << "- Number of events available: " << nPuEvts  << std::endl;

  double lArea[nLayers];
  for (unsigned iL(0);iL<nLayers;++iL){
    lArea[iL] = areaInCm2(radius(eta[0]-deta,posz[iL]))-areaInCm2(radius(eta[0]+deta,posz[iL]));
    std::cout << " Layer " << iL << " area = " << lArea[iL] << " cm^2 r1=" 
	      << radius(eta[0]-deta,posz[iL]) << " r2=" << radius(eta[0]+deta,posz[iL])
	      << std::endl;
  }

  std::vector<KDNode> _hit_nodes;
  KDTree _hit_kdtree;


  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
    if (debug) std::cout << "... Processing entry: " << ievt << std::endl;
    else if (ievt%100 == 0) std::cout << "... Processing entry: " << ievt << std::endl;

    lTree->GetEntry(ievt);
    std::vector<MyHit> lHitVec;

    unsigned nRandomCones = 100;
    double phistep = 2*TMath::Pi()/nRandomCones;
    double phi0 = lRndm.Uniform(-1.*TMath::Pi(),TMath::Pi());
    double eta0 = lRndm.Uniform(eta[0]-deta,eta[0]+deta);
    //if (debug_) std::cout << "--- etamax = " << etamax << " phimax=" << phimax << " phistep = " << phistep << std::endl;

    std::vector<unsigned> nAboveTot;
    nAboveTot.resize(nLayers,0);
    std::vector<double> EperLayer;
    EperLayer.resize(nLayers,0);

    for (unsigned ipm(0);ipm<nRandomCones;++ipm){
      std::vector<double> xmax;
      xmax.resize(nLayers,0);
      std::vector<double> ymax;
      ymax.resize(nLayers,0);
      double phirc = phi0;
      if (ipm%2==0) phirc += ipm/2*phistep+phistep/2.;
      else  phirc = phirc - ipm/2*phistep-phistep/2.;
      if (phirc < -1.*TMath::Pi()) phirc+=2.*TMath::Pi();
      if (phirc > TMath::Pi()) phirc-=2.*TMath::Pi();
      getMaximumCellFromGeom(phirc,eta0,xmax,ymax,posz);

      unsigned nover = 0;
      unsigned ntot = 0;
      for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on hits
	const HGCSSRecoHit & lHit = (*rechitvec)[iH];

	double posx = lHit.get_x();
	double posy = lHit.get_y();
	unsigned layer = lHit.layer();
	double leta = lHit.eta();
	double energy = lHit.energy();

	if (ipm==0 && fabs(leta-eta[0])<deta){
	  if (energy>Ethresh) nAboveTot[layer]++;
	  EperLayer[layer]+=energy;
	  MyHit lHit;
	  lHit.e = energy;
	  if (lHit.e==0) continue;
	  lHit.layer = layer;
	  lHit.x = posx;
	  lHit.y = posy;
	  lHit.z = posz[layer];
	  lHit.hasSignal = false;
	  lHitVec.push_back(lHit);
	}

	if (fabs(posx-xmax[layer]) < 15.1 && 
	    fabs(posy-ymax[layer]) < 15.1){
	  ntot++;
	  hxy->Fill(posx,posy);
	  hphieta->Fill(lHit.phi(),leta);
	  if (energy>Ethresh){
	    hxyAbove->Fill(posx,posy);
	    hphietaAbove->Fill(lHit.phi(),leta);
	    nover++;
	  }
	}
      }

      nAbove->Fill(nover);
      nTot->Fill(ntot);

    }//loop on random cones

    fillKDTree(_hit_kdtree,
	       _hit_nodes,
	       lHitVec);
	
    if (debug) std::cout << " ----- Number of hits in vector : " << lHitVec.size() << std::endl;
    bool fillShowerProfile = false;
    for (unsigned iH(0); iH<lHitVec.size();++iH){
      const MyHit& current = lHitVec[iH];
      unsigned iL = current.layer;
      
      std::vector<KDNode> neighbours;
      findNeighbours(cell_size,_hit_kdtree,current,neighbours,2);
      
      if (debug>1) std::cout << " -- Number of closest neighbours found: " << neighbours.size() << std::endl; 
      //if (neighbours.size()>nMaxNeigh) nMaxNeigh=neighbours.size();
      
      double EmaxSC3D = 0;
      for( const KDNode& nbourpoint :neighbours ) {
	unsigned index = nbourpoint.data;
	if (index == iH) {
	  //neighbours.erase(neighbours.begin()+iH);
	  continue;
	}
	const MyHit& nbour = lHitVec[index];
	if ((nbour.layer==iL  && (nbour.x==current.x || nbour.y==current.y )) || (abs(nbour.layer-iL)<2 && nbour.x==current.x && nbour.y==current.y )){
	  if (nbour.e>EmaxSC3D) EmaxSC3D = nbour.e;
	  //if (nbour.e>=EmaxCut) break;
	}
      }
      if (current.e>Ethresh) {
	EmaxNeighbourSC3D->Fill(EmaxSC3D);
	EneighvsEcell->Fill(current.e,EmaxSC3D);
	if (EmaxSC3D>50) fillShowerProfile = true;
      }
    }

    for (unsigned iL(0);iL<nLayers;++iL){
      if (fillShowerProfile) hLayerProfile_above50mips->Fill(iL,EperLayer[iL]);
      nAbovevsLayer->Fill(iL,nAboveTot[iL]/lArea[iL]);
      prof_nAbovevsLayer->Fill(iL,nAboveTot[iL]/lArea[iL]);
    }
  }//loop on entries

  outputFile->cd();
  outputFile->Write();
  outputFile->Close();

  return 0;
  
}//main
