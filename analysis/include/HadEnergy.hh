#ifndef HadEnergy_h
#define HadEnergy_h

#include<string>
#include<iostream>
#include<fstream>
#include<sstream>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TAxis.h"
#include "TMath.h"
#include "TProfile.h"

#include "HGCSSRecoHit.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSGeometryConversion.hh"
#include "HGCSSCalibration.hh"
#include "HGCSSDetector.hh"

class HadEnergy{

public:
  HadEnergy(const std::vector<unsigned> & dropLay,
	    HGCSSDetector & myDetector, 
            TChain* lSimTree,
            TChain* lRecTree,
            TFile *outputFile,
            const unsigned pNevts);
  ~HadEnergy();

  inline void addLimMIP(const double lim){
    if(lim >0)LimMIP_.push_back(lim);
  };

  void bookHist(TFile *outputFile);

  inline void setOutputFile(const std::string outputName){
       outputName_ = outputName;
       TFile *outputFile_ =  TFile::Open(outputName.c_str(),"RECREATE");
       bookHist(outputFile_); 
  };
  
  inline void setFHtoE(const double FHtoEslope,const double FHtoEoffset){
       FHtoEslope_ = FHtoEslope;
       FHtoEoffset_ = FHtoEoffset;
  };

  inline void setBHtoE(const double BHtoEslope,const double BHtoEoffset){
       BHtoEslope_ = BHtoEslope;
       BHtoEoffset_ = BHtoEoffset;
  };

  inline void setFHtoBH(const double FHtoBHslope){
    FHtoBHslope_ = FHtoBHslope;
  };

  inline void setEEtoH(const double EEtoHslope){
    EEtoHslope_ = EEtoHslope;
  };

  inline void setEEtoHPar(const double EEtoHslopePar0,const double EEtoHslopePar1){
    EEtoHslopePar0_ = EEtoHslopePar0;
    EEtoHslopePar1_ = EEtoHslopePar1;
  };

  inline void setEEcalib(const double ECALslope,const double ECALoffset){
       ECALslope_ = ECALslope;
       ECALoffset_ = ECALoffset;
  };

  bool fillEnergies();
 
  double calcGlobalC(const double LimMIP, const double EmipMean, TH1F* spectrum); 
  
  inline double getGlobalC(const unsigned iLim) { 
      if(iLim < Cglobal_.size()) return Cglobal_[iLim];
      else return 0;
  };

private:

  std::string outputName_;
  TFile *outputFile_;
  TTree *outtree_;
  HGCSSDetector & myDetector_;
  TChain* lSimTree_;
  TChain* lRecTree_;
  unsigned pNevts_; 
  std::vector<unsigned> dropLay_;

  unsigned debug_; 
  //for tree
  double nLayers_;
  unsigned version_;
  unsigned nSections_;
  std::vector<double> absweight_;
  unsigned ievt_;
  //for tree
  double wgttotalE_;
  std::vector<double> correctedtotalE_;
  std::vector<double> energy_;

  std::vector<double> LimMIP_;
  std::vector<double> Cglobal_; 
  double EE_, EFHCAL_, EBHCAL_;
  double EmipMeanFH_;
  unsigned nhitsFH_; 

  TH1F *p_spectrum;
  TH1F *p_spectrum_hightail;
  TH1F *p_spectrum_lowtail;
  TH2F *p_spectrumByLayer;

  bool isCalibed_;
  double FHtoEslope_; 
  double FHtoEoffset_;
  double BHtoEslope_;
  double BHtoEoffset_;
  double ECALslope_;
  double ECALoffset_;
  double FHtoBHslope_;
  double EEtoHslope_;
  double EEtoHslopePar0_;
  double EEtoHslopePar1_;
};

#endif
