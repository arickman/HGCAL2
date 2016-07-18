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
  unsigned debug;
  unsigned nPu;
  double etamean;

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
    ("etamean,e",      po::value<double>(&etamean)->default_value(2.85))
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
  const unsigned nNoise = 1;//5;//10;

  //const double noise[nNoise] = {0,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5};
  const double noise[nNoise] = {0.0};//0.6,0.7,0.8,0.9,1.0};
  double eta[nEta];// = {2.85};//1.7,2.0,2.5};

  eta[0] = etamean;

  const double deta = 0.05;

  //const double Ethresh = 0.6;

  /////////////////////////////////////////////////////////////
  //input
  /////////////////////////////////////////////////////////////

  std::ostringstream inputbase;
  inputbase << inFilePath //<< "/MipStudy.root";
	    << "/MipStudy_eta" << (eta[0]-0.15) << "_" << (eta[0]+0.15);
  TChain *lTree = new TChain("RecoTree");
  bool first = true;
  HGCSSInfo * info = 0;

  for (unsigned i(0); i<11;++i){
    unsigned start = i*50+1;
    unsigned end = start+49;
    if (i==0) start = 26;
    if (i==10) end = 525;
    std::ostringstream input;
    input << inputbase.str() << "_MinBias" << start << "_" << end 
	  << ".root";
    lTree->AddFile(input.str().c_str());
    std::cout << "Adding MinBias file:" << input.str().c_str() << std::endl;
    if (first){
      TFile *simFile = TFile::Open(input.str().c_str());
      if (simFile) {
	info =(HGCSSInfo*)simFile->Get("Info");
	if (info) first=false;
      }
    }
  }
  
  //TFile *simFile = 0;
  //if (!testInputFile(input.str(),simFile)) return 1;
  //TTree *lTree = (TTree*)simFile->Get("RecoTree");
  if (!lTree){
    std::cout << " -- Error, tree RecoTree cannot be opened. Exiting..." << std::endl;
    return 1;
  }

  std::cout << " Trees added." << std::endl;

  /////////////////////////////////////////////////////////////
  //Info
  /////////////////////////////////////////////////////////////

  //HGCSSInfo * info=(HGCSSInfo*)simFile->Get("Info");
  const double cellSize = info->cellSize();
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


  HGCSSGeometryConversion geomConv(inFilePath,model,cellSize);
  //set granularity to get cellsize for PU subtraction
  std::vector<unsigned> granularity;
  granularity.resize(nLayers,4);//0.5*0.5 cells
  geomConv.setGranularity(granularity);
  geomConv.initialiseHistos();

  TRandom3 lRndm(0);
  std::cout << " -- Random number seed: " << lRndm.GetSeed() << std::endl;

  //////////////////////////////////////////////////
  //////////////////////////////////////////////////
  ///////// Output File // /////////////////////////
  //////////////////////////////////////////////////
  //////////////////////////////////////////////////

  TFile *outputFile = TFile::Open(outFilePath.c_str(),"RECREATE");
  
  if (!outputFile) {
    std::cout << " -- Error, output file " << outFilePath << " cannot be opened. Please create output directory. Exiting..." << std::endl;
    return 1;
  }
  else {
    std::cout << " -- output file " << outputFile->GetName() << " successfully opened." << std::endl;
  }
  outputFile->cd();

  /*  TH1F *p_EmaxNeighbour[nNoise];

  for (unsigned in(0); in<nNoise;++in){
    std::ostringstream label;
    label << "EmaxNeighbour_noise" << noise[in]*100;
    p_EmaxNeighbour[in] = new TH1F(label.str().c_str(),";E_{max}^{neighbour};hits",100,0,5);
  }

  TH1F *hitSpectrum[nEta][nNoise][nLayers];
  TH1F *p_nHits[nEta][nNoise][nLayers];
  for (unsigned ie(0);ie<nEta;++ie){
    for (unsigned in(0); in<nNoise;++in){
      for (unsigned il(0); il<nLayers;++il){
	std::ostringstream label;
	label << "hitSpectrum_eta" << eta[ie]*10 << "_noise" << noise[in]*100 << "_layer" << il;
	hitSpectrum[ie][in][il] = new TH1F(label.str().c_str(),";E_{1#times1 cm^{2}} (mips);N_{cells}",100,0,5);
	label.str("");
	label << "Nhits_eta" << eta[ie]*10 << "_noise" << noise[in]*100 << "_layer" << il;
	p_nHits[ie][in][il] = new TH1F(label.str().c_str(),";N_{1#times1 cm^{2}};N_{events}",500,0,5000);
      }
    }
  }
  */
  TTree *outtree = new TTree("MipStudy","HGC standalone simulation mip study variables");
  std::vector<unsigned> cellIdsX[nEta][nNoise];
  std::vector<unsigned> cellIdsY[nEta][nNoise];
  std::vector<unsigned> cellIdsZ[nEta][nNoise];
  std::vector<double> energies[nEta][nNoise];
  std::vector<bool> signal[nEta][nNoise];
  for (unsigned ie(0);ie<nEta;++ie){
    for (unsigned in(0); in<nNoise;++in){
      std::ostringstream label;
      label << "X_" << ie << "_" << in;
      outtree->Branch(label.str().c_str(),&cellIdsX[ie][in]);
      label.str("");
      label << "Y_" << ie << "_" << in;
      outtree->Branch(label.str().c_str(),&cellIdsY[ie][in]);
      label.str("");
      label << "Z_" << ie << "_" << in;
      outtree->Branch(label.str().c_str(),&cellIdsZ[ie][in]);
      label.str("");
      label << "E_" << ie << "_" << in;
      outtree->Branch(label.str().c_str(),&energies[ie][in]);
      label.str("");
      label << "signal_" << ie << "_" << in;
      outtree->Branch(label.str().c_str(),&signal[ie][in]);
    }
  }

  const unsigned nEvts = ((pNevts > lTree->GetEntries() || pNevts==0) ? static_cast<unsigned>(lTree->GetEntries()) : pNevts) ;

  std::vector<HGCSSRecoHit> * rechitvec = 0;
  lTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);

  //loop on events
  const unsigned nPuEvts = lTree->GetEntries();
  std::cout << "- Number of PU events available: " << nPuEvts  << std::endl;

  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
    if (debug) std::cout << "... Processing entry: " << ievt << std::endl;
    else if (ievt%10 == 0) std::cout << "... Processing entry: " << ievt << std::endl;

    //lTree->GetEntry(ievt);
    //get poisson <140>
    unsigned nPuVtx = lRndm.Poisson(nPu);
    //for single interactions
    if (nPu==0) nPuVtx = 1;
 
    if (debug) std::cout << " -- Adding " << nPuVtx << " events to signal event: " << ievt << std::endl;
    for (unsigned iV(0); iV<nPuVtx; ++iV){//loop on interactions
      unsigned ipuevt = lRndm.Integer(nPuEvts);
      if (nPu==0) ipuevt = ievt;
      lTree->GetEntry(ipuevt);

      unsigned prevLayer = 10000;
      DetectorEnum type = DetectorEnum::FECAL;
      unsigned subdetLayer=0;
      for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on hits
	HGCSSRecoHit lHit = (*rechitvec)[iH];
	//double posz = lHit.get_z();
	double leta = lHit.eta();
	if (fabs(leta-eta[0])>= deta) continue;
	// && posz>minZ && posz<maxZ;
	unsigned layer = lHit.layer();
	if (layer != prevLayer){
	  const HGCSSSubDetector & subdet = myDetector.subDetectorByLayer(layer);
	  type = subdet.type;
	  subdetLayer = layer-subdet.layerIdMin;
	  prevLayer = layer;
	  //std::cout << " - layer " << layer << " " << subdet.name << " " << subdetLayer << std::endl;
	}      
	double energy = lHit.energy();
	double posx = lHit.get_x();
	double posy = lHit.get_y();
	double posz = lHit.get_z();
	geomConv.fill(type,subdetLayer,energy,0,posx,posy,posz);
	
      }//loop on hits
    }//loop on interactions

    unsigned nTotBins = 0;
    std::map<std::pair<int,int>,double> lMap[nEta][nNoise][nLayers];
  
    for (unsigned iL(0); iL<nLayers; ++iL){//loop on layers
      TH2D *histE = geomConv.get2DHist(iL,"E");
      nTotBins += histE->GetNbinsX()*histE->GetNbinsY();
      
      if (debug>1) std::cout << " Nhits in layer " << iL << " : " << histE->GetNbinsX()*histE->GetNbinsY() << std::endl;

      for (int iX(1); iX<histE->GetNbinsX()+1;++iX){
	for (int iY(1); iY<histE->GetNbinsY()+1;++iY){
	    double x = histE->GetXaxis()->GetBinCenter(iX);
	    double y = histE->GetYaxis()->GetBinCenter(iY);
	    double z = zPos(iL);
	    ROOT::Math::XYZPoint position(x,y,z);
	    double leta = position.Eta();
	    //consider only 3 eta rings: add noise and digitise
	    //if (!(fabs(leta-eta[ieta])<deta ||
	    //fabs(leta-2.0)<deta ||
	    //fabs(leta-2.5)<deta)) continue;
	    double simE = histE->GetBinContent(iX,iY);
	    bool isSig = false;
	    if (simE>0) isSig = true;
	    //add noise
	    for (unsigned ie(0);ie<nEta;++ie){
	      if (fabs(leta-eta[ie])>= deta) continue;
	      for (unsigned in(0); in<nNoise;++in){
		double lNoise = lRndm.Gaus(0,noise[in]);
		//double digiE = static_cast<unsigned>((simE+lNoise)*4)/4.;
		double digiE = simE+lNoise;
		if (digiE>0) {
		  lMap[ie][in][iL].insert(std::pair<std::pair<int,int>,double>(std::pair<int,int>(iX,iY),digiE));
		  cellIdsX[ie][in].push_back(iX);
		  cellIdsY[ie][in].push_back(iY);
		  cellIdsZ[ie][in].push_back(iL);
		  energies[ie][in].push_back(digiE);
		  signal[ie][in].push_back(isSig);
		}
	      }
	    }//loop on eta
	}//loop on y
      }//loop on x
    }//loop on layers

    outtree->Fill();
    
    //unsigned nHits[nEta][nNoise][nLayers];
    for (unsigned ie(0);ie<nEta;++ie){
      for (unsigned in(0); in<nNoise;++in){
	cellIdsX[ie][in].clear();
	cellIdsY[ie][in].clear();
	cellIdsZ[ie][in].clear();
	energies[ie][in].clear();
	signal[ie][in].clear();
	/*
	  for (unsigned iL(0); iL<nLayers; ++iL){//loop on layers
	  nHits[ie][in][iL] = 0;
	  TH2D *histE = geomConv.get2DHist(iL,"E");
	  if (debug>1) std::cout << " -- Number of hits in map[" << in << "] = " << lMap[ie][in][iL].size() << std::endl;
	  std::map<std::pair<int,int>,double>::iterator iter = lMap[ie][in][iL].begin();
	  for (;iter != lMap[ie][in][iL].end();++iter){
	    double Ehit = iter->second;
	    if (Ehit<Ethresh) continue;
	    int iX = iter->first.first;
	    int iY = iter->first.second;

	    //look at swiss-cross neighbours
	    double E0=0,E1=0,E2=0,E3=0;
	    if (lMap[ie][in][iL].count(std::pair<int,int>(iX,iY-1))>0) E0 = lMap[ie][in][iL][std::pair<int,int>(iX,iY-1)];
	    if (lMap[ie][in][iL].count(std::pair<int,int>(iX,iY+1))>0) E1 = lMap[ie][in][iL][std::pair<int,int>(iX,iY+1)];
	    if (lMap[ie][in][iL].count(std::pair<int,int>(iX-1,iY))>0) E2 = lMap[ie][in][iL][std::pair<int,int>(iX-1,iY)];
	    if (lMap[ie][in][iL].count(std::pair<int,int>(iX+1,iY))>0) E3 = lMap[ie][in][iL][std::pair<int,int>(iX+1,iY)];
	    double Emax = std::max(std::max(E0,E1),std::max(E2,E3));
	    p_EmaxNeighbour[in]->Fill(Emax);

	    //look at before and after layers
	    double x = histE->GetXaxis()->GetBinCenter(iX);
	    double y = histE->GetYaxis()->GetBinCenter(iY);

	    if (debug>1) std::cout << " -- " 
				   << iL << ", "
				   << iX << "-" << x << ", "
				   << iY << "-" << y << ", "
				   << "Emax = " << Emax << std::endl;

	    double z = zPos(iL);
	    ROOT::Math::XYZPoint position(x,y,z);
	    double ltheta = position.Theta();
	    double lphi = position.Phi();

	    double dxp = iL < (nLayers-1) ? (zPos(iL+1)-z)*tan(ltheta)*cos(lphi) : 0;
	    double dxm = iL>1? (zPos(iL-1)-z)*tan(ltheta)*cos(lphi) : 0;
	    double dyp = iL < (nLayers-1) ? (zPos(iL+1)-z)*tan(ltheta)*sin(lphi) : 0;
	    double dym = iL>1 ? (zPos(iL-1)-z)*tan(ltheta)*sin(lphi) : 0;
	    //track in previous and next layers
	    bool track1 = true;
	    if (iL>1){
	      if (lMap[ie][in][iL-1].count(std::pair<int,int>(iX,iY))>0) track1 = lMap[ie][in][iL-1][std::pair<int,int>(iX,iY)]>Ethresh;
	      if ( dym<0 ) {
		if ( lMap[ie][in][iL-1].count(std::pair<int,int>(iX,iY-1))>0 ) track1 = track1 || lMap[ie][in][iL-1][std::pair<int,int>(iX,iY-1)]>Ethresh;
		if ( dxm<0 ){
		  if (lMap[ie][in][iL-1].count(std::pair<int,int>(iX-1,iY-1))>0 ) track1 = track1 || lMap[ie][in][iL-1][std::pair<int,int>(iX-1,iY-1)]>Ethresh;
		  if (lMap[ie][in][iL-1].count(std::pair<int,int>(iX-1,iY))>0 ) track1 = track1 || lMap[ie][in][iL-1][std::pair<int,int>(iX-1,iY)]>Ethresh;
		}
		else if ( dxm>0 ){
		  if (lMap[ie][in][iL-1].count(std::pair<int,int>(iX+1,iY-1))>0 )  track1 = track1 || lMap[ie][in][iL-1][std::pair<int,int>(iX+1,iY-1)]>Ethresh;
		  if (lMap[ie][in][iL-1].count(std::pair<int,int>(iX+1,iY))>0 )  track1 = track1 || lMap[ie][in][iL-1][std::pair<int,int>(iX+1,iY)]>Ethresh;
		}
	      }
	      else if ( dym>0 ){
		if ( lMap[ie][in][iL-1].count(std::pair<int,int>(iX,iY+1))>0 ) track1 = track1 || lMap[ie][in][iL-1][std::pair<int,int>(iX,iY+1)]>Ethresh;
		if ( dxm<0 ){
		  if (lMap[ie][in][iL-1].count(std::pair<int,int>(iX-1,iY+1))>0 ) track1 = track1 || lMap[ie][in][iL-1][std::pair<int,int>(iX-1,iY+1)]>Ethresh;
		  if (lMap[ie][in][iL-1].count(std::pair<int,int>(iX-1,iY))>0 ) track1 = track1 || lMap[ie][in][iL-1][std::pair<int,int>(iX-1,iY)]>Ethresh;
		}
		else if ( dxm>0 ){
		  if (lMap[ie][in][iL-1].count(std::pair<int,int>(iX+1,iY+1))>0 )  track1 = track1 || lMap[ie][in][iL-1][std::pair<int,int>(iX+1,iY+1)]>Ethresh;
		  if (lMap[ie][in][iL-1].count(std::pair<int,int>(iX+1,iY))>0 )  track1 = track1 || lMap[ie][in][iL-1][std::pair<int,int>(iX+1,iY)]>Ethresh;
		}
	      }
	    }
	    bool track2 = true;
	    if (iL < (nLayers-1) )
	      {
		if (lMap[ie][in][iL+1].count(std::pair<int,int>(iX,iY))>0) track2 = lMap[ie][in][iL+1][std::pair<int,int>(iX,iY)]>Ethresh;
		if ( dyp<0 ) {
		  if ( lMap[ie][in][iL+1].count(std::pair<int,int>(iX,iY-1))>0 ) track2 = track2 || lMap[ie][in][iL+1][std::pair<int,int>(iX,iY-1)]>Ethresh;
		  if ( dxp<0 ){
		    if (lMap[ie][in][iL+1].count(std::pair<int,int>(iX-1,iY-1))>0 ) track2 = track2 || lMap[ie][in][iL+1][std::pair<int,int>(iX-1,iY-1)]>Ethresh;
		    if (lMap[ie][in][iL+1].count(std::pair<int,int>(iX-1,iY))>0 ) track2 = track2 || lMap[ie][in][iL+1][std::pair<int,int>(iX-1,iY)]>Ethresh;
		  }
		  else if ( dxp>0 ){
		    if (lMap[ie][in][iL+1].count(std::pair<int,int>(iX+1,iY-1))>0 )  track2 = track2 || lMap[ie][in][iL+1][std::pair<int,int>(iX+1,iY-1)]>Ethresh;
		    if (lMap[ie][in][iL+1].count(std::pair<int,int>(iX+1,iY))>0 )  track2 = track2 || lMap[ie][in][iL+1][std::pair<int,int>(iX+1,iY)]>Ethresh;
		  }
		}
		else if ( dyp>0 ){
		  if ( lMap[ie][in][iL+1].count(std::pair<int,int>(iX,iY+1))>0 ) track2 = track2 || lMap[ie][in][iL+1][std::pair<int,int>(iX,iY+1)]>Ethresh;
		  if ( dxp<0 ){
		    if (lMap[ie][in][iL+1].count(std::pair<int,int>(iX-1,iY+1))>0 ) track2 = track2 || lMap[ie][in][iL+1][std::pair<int,int>(iX-1,iY+1)]>Ethresh;
		    if (lMap[ie][in][iL+1].count(std::pair<int,int>(iX-1,iY))>0 ) track2 = track2 || lMap[ie][in][iL+1][std::pair<int,int>(iX-1,iY)]>Ethresh;
		  }
		  else if ( dxp>0 ){
		    if (lMap[ie][in][iL+1].count(std::pair<int,int>(iX+1,iY+1))>0 )  track2 = track2 || lMap[ie][in][iL+1][std::pair<int,int>(iX+1,iY+1)]>Ethresh;
		    if (lMap[ie][in][iL+1].count(std::pair<int,int>(iX+1,iY))>0 )  track2 = track2 || lMap[ie][in][iL+1][std::pair<int,int>(iX+1,iY)]>Ethresh;
		  }
		}
		
	      }
	    if (track1 && track2){
	      hitSpectrum[ie][in][iL]->Fill(Ehit);
	      nHits[ie][in][iL]++;
	    }
	  }//loop on map elements
	  p_nHits[ie][in][iL]->Fill(nHits[ie][in][iL]);
	}//loop on layers
	*/
      }//loop on noise
    }//loop on eta

    geomConv.initialiseHistos();
    
  }//loop on entries

  outputFile->cd();
  outputFile->Write();
  outputFile->Close();

  return 0;
  
}//main
