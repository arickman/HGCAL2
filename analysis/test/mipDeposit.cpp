#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include <boost/algorithm/string.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TStyle.h"

#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TString.h"
#include "TLatex.h"

#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSParameters.hh"


int main(int argc, char** argv){//main


  TString lSuffix = "version23_e100";
  TString plotBase = "PLOTS/version23/mu-/";

  TFile *inputFile = TFile::Open("root://eoscms//eos/cms/store/user/amagnan/HGCalHEGeant4/run_0/mu-/HGcal_"+lSuffix+".root");
  if (!inputFile) {
    std::cout << " -- Error, input file cannot be opened. Exiting..." << std::endl;
    return 1;
  }
  TTree *lTree = (TTree*)inputFile->Get("HGCSSTree");
  if (!lTree){
    std::cout << " -- Error, tree cannot be opened. Exiting..." << std::endl;
    return 1;
  }
  
  TCanvas *myc = new TCanvas("myc","myc",500,500);

  const unsigned nLayers = 54;//33;
  const unsigned nHcalSiLayers = 0;//24;  

  TH2F *p_nHits = new TH2F("nHits","; layer; Number of hits; Events",nLayers,0,nLayers,20,0,20);
  TH1F *p_hitEnergy_si = new TH1F("hitEnergy_si",";E (MeV);SimHits",250,0,1);
  TH1F *p_hitEnergySel_si = new TH1F("hitEnergySel_si",";E (MeV);SimHits",250,0,1);
  TH1F *p_hitEnergy_scint = new TH1F("hitEnergy_scint",";E (MeV);SimHits",500,0,100);
  TH1F *p_hitEnergySel_scint = new TH1F("hitEnergySel_scint",";E (MeV);SimHits",1000,0,10);

  std::vector<HGCSSSimHit> * simhitvec = 0;
  float volNb = 0;
  lTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
  lTree->SetBranchAddress("volNb",&volNb);

  const unsigned nEvts = lTree->GetEntries();

  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries

    lTree->GetEntry(ievt);
    
    if (ievt%(nLayers*1000)==0) std::cout << "entry " << ievt << " volNb = " << volNb << std::endl;

    unsigned nHits = 0;
    double energySel = 0;
    for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){//loop on hits
      HGCSSSimHit lHit = (*simhitvec)[iH];
      double energy = lHit.energy();
      if (energy>0) {
	if (volNb < nHcalSiLayers) p_hitEnergy_si->Fill(energy);
	else p_hitEnergy_scint->Fill(energy);
	nHits++;
      }
    }//loop on hits

    p_nHits->Fill(volNb,nHits);

    if (nHits>0 && nHits<2) {
      for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){//loop on hits
	HGCSSSimHit lHit = (*simhitvec)[iH];
	double energy = lHit.energy();
	if (energy>0){
	  if (volNb < nHcalSiLayers) p_hitEnergySel_si->Fill(energy);
	  else  p_hitEnergySel_scint->Fill(energy);
	}
      }
    }

  }//loop on entries

  
  myc->cd();
  gPad->SetLogz(1);
  gStyle->SetOptStat(1111110);
  p_nHits->Draw("colz");

  myc->Update();
  myc->Print(plotBase+"/mipHits.png");
  myc->Print(plotBase+"/mipHits.pdf");
  myc->Print(plotBase+"/mipHits.C");


  myc->cd();
  gPad->SetLogy(1);
  gStyle->SetOptStat(1111110);
  gStyle->SetOptFit(1111);
  p_hitEnergy_si->Draw();
  p_hitEnergy_si->Fit("landau","R+","",0.035,1);

  myc->Update();
  myc->Print(plotBase+"/mipDepositAll_si.png");
  myc->Print(plotBase+"/mipDepositAll_si.pdf");
  myc->Print(plotBase+"/mipDepositAll_si.C");

  myc->cd();
  gPad->SetLogy(1);
  gStyle->SetOptStat(1111110);
  gStyle->SetOptFit(1111);
  p_hitEnergySel_si->Draw();
  p_hitEnergySel_si->Fit("landau","LR+","",0.02,1);

  myc->Update();
  myc->Print(plotBase+"/mipDepositSel_si.png");
  myc->Print(plotBase+"/mipDepositSel_si.pdf");
  myc->Print(plotBase+"/mipDepositSel_si.C");


  myc->cd();
  gPad->SetLogy(1);
  gStyle->SetOptStat(1111110);
  gStyle->SetOptFit(1111);
  p_hitEnergy_scint->Draw();
  p_hitEnergy_scint->Fit("landau","R+","",0.035,1);
  
  myc->Update();
  myc->Print(plotBase+"/mipDepositAll_scint.png");
  myc->Print(plotBase+"/mipDepositAll_scint.pdf");
  myc->Print(plotBase+"/mipDepositAll_scint.C");

  myc->cd();
  gPad->SetLogy(1);
  gStyle->SetOptStat(1111110);
  gStyle->SetOptFit(1111);
  p_hitEnergySel_scint->Draw();
  p_hitEnergySel_scint->Fit("landau","LR+","",0.02,1);
  
  myc->Update();
  myc->Print(plotBase+"/mipDepositSel_scint.png");
  myc->Print(plotBase+"/mipDepositSel_scint.pdf");
  myc->Print(plotBase+"/mipDepositSel_scint.C");



  return 0;

}//main
