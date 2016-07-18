#ifndef HiggsMass_h
#define HiggsMass_h



#include<string>
#include<iostream>
#include<fstream>
#include<sstream>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"

#include "HGCSSRecoHit.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSPUenergy.hh"
#include "SignalRegion.hh"
#include "PositionFit.hh"

#include "Math/Vector4D.h"
#include "Math/Vector4Dfwd.h"
#include "Math/Point3D.h"
#include "Math/Point3Dfwd.h"

#include "TLorentzVector.h"
#include "TRandom3.h"

bool CheckCosTheta(const TLorentzVector & v){
  double cosTheta = v.CosTheta();
  if (cosTheta*cosTheta < 1) return true;
  return false;
};

class HiggsMass{

public:
  HiggsMass(){
    rand_.SetSeed(0);
  };

  ~HiggsMass(){};

  void setRecoInfo(const TLorentzVector & g1,
		   const TLorentzVector & g2,
		   const ROOT::Math::XYZPoint & posFF1,
		   const ROOT::Math::XYZPoint & posFF2);
    
  void setTruthInfo(const TLorentzVector & g1,
		    const TLorentzVector & g2,
		    const ROOT::Math::XYZPoint & vtx1,
		    const ROOT::Math::XYZPoint & vtx2);
		    
  void initialiseHistograms(TFile *fout, 
			    const std::string folder);

  void fillHistograms();


private:

  TRandom3 rand_;

  //reco
  TLorentzVector g1_;
  TLorentzVector g2_;

  ROOT::Math::XYZPoint posFF1_;
  ROOT::Math::XYZPoint posFF2_;

  //truth
  TLorentzVector tg1_;
  TLorentzVector tg2_;
  ROOT::Math::XYZPoint tvtx1_;
  ROOT::Math::XYZPoint tvtx2_;


  //from true dir, true E
  //from true direction, recoE
  TH1F *p_trueDir_trueE;
  TH1F *p_trueDir_recoE;
  //from projected positionFF+true vertex, trueE
  //from projected positionFF+true vertex, recoE
  TH1F *p_position_trueE;
  TH1F *p_position_recoE;
  //from reco direction, trueE
  //from reco direction, recoE
  TH1F *p_angle_trueE;
  TH1F *p_angle_recoE;
  //from projected positionFF+smeared vertex, trueE
  //from projected positionFF+smeared vertex, recoE
  TH1F *p_position_vtxsmear_trueE;
  TH1F *p_position_vtxsmear_recoE;

  TH1F *p_ErecooverEtrue;

  TH1F *p_vtx_x;
  TH1F *p_vtx_y;
  TH1F *p_vtx_z;

  TH1F *p_dvtx_x;
  TH1F *p_dvtx_y;
  TH1F *p_dvtx_z;

  TH1F *p_pT_Higgs;
  TH1F *p_eta_Higgs;

  TH2F *p_pTvseta[2];
  TH2F *p_MvspT[2];

  TH1F *p_pT_gamma1[2];
  TH1F *p_eta_gamma1[2];
  TH1F *p_phi_gamma1[2];

  TH1F *p_pT_gamma2[2];
  TH1F *p_eta_gamma2[2];
  TH1F *p_phi_gamma2[2];

};//class


#endif
