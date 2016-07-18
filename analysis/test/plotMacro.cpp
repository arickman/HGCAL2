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

  const double Emip = 0.0548;//in MeV
  unsigned genEn[]={10};
  const unsigned nGenEn=sizeof(genEn)/sizeof(unsigned);
  
  TCanvas *myc = new TCanvas("myc","myc",1);//1000,500);
  //myc->Divide(2,1);
  
  for (unsigned iE(0); iE<nGenEn; ++iE){
    std::cout << "- Processing energy : " << genEn[iE] 
	      << std::endl;

    TString genEnStr = "";
    genEnStr += genEn[iE];
    
    TString lSuffix[2] = {"3x100","200"};

    TH1F *Etot[2];
    Etot[0] = new TH1F("Etotal_2x100",";E (MeV)",500,0,100);
    Etot[1] = new TH1F("Etotal_200",";E (MeV)",500,0,100);

    for (unsigned iS(0); iS<2;++iS){

      TFile *inputFile = TFile::Open("/afs/cern.ch/work/a/amagnan/public/HGCalEEGeant4/git_V00-02-05/version_8/model_1/e-/BOFF/e_"+genEnStr+"/PFcal_"+lSuffix[iS]+".root");
      if (!inputFile) {
	std::cout << " -- Error, input file cannot be opened. Exiting..." << std::endl;
	return 1;
      }
      TTree *lTree = (TTree*)inputFile->Get("HGCSSTree");
      if (!lTree){
	std::cout << " -- Error, tree RecoTree cannot be opened either. Exiting..." << std::endl;
	return 1;
      }
    
      TH1F *simhitEnergy = new TH1F("simhitEnergy",";E (MIPs)",750,0,1500);
    
      std::vector<HGCSSSimHit> * simhitvec = 0;
      lTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
      
      const unsigned nEvts = lTree->GetEntries();
      
      double maxEhit = 0;
      
      for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
	if (ievt%10==0) std::cout << "... Processing entry: " << ievt << std::endl;
	
	lTree->GetEntry(ievt);

	double Etotal = 0;
	for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){//loop on hits
	  HGCSSSimHit lHit = (*simhitvec)[iH];
	  if (lHit.energy()>maxEhit) maxEhit = lHit.energy();
	  if (lHit.silayer()<2){
	    simhitEnergy->Fill(lHit.energy());///Emip);
	    Etotal+=lHit.energy();
	  }
	}//loop on hits
      
	Etot[iS]->Fill(Etotal);

      }//loop on entries
      
      std::cout << " -- max hit energy = " << maxEhit << std::endl;
      
      myc->cd();
      TLatex lat;
      float yMax = 0;
      gPad->SetLogy(1);
      yMax = simhitEnergy->GetMaximum();
      simhitEnergy->GetXaxis()->SetRangeUser(0,maxEhit);
      simhitEnergy->Draw();
      lat.DrawLatex(0,yMax*2,"SIM e-, "+genEnStr+" GeV, 2.5#times2.5 mm^{2} cells");
      
      myc->Update();
      myc->Print("PLOTS/SimHitEnergy_"+genEnStr+"GeV"+lSuffix[iS]+".png");
      myc->Print("PLOTS/SimHitEnergy_"+genEnStr+"GeV"+lSuffix[iS]+".pdf");
      
    }//loop on scenario

    myc->cd();
    Etot[0]->SetLineColor(2);
    Etot[0]->SetMarkerColor(2);
    Etot[0]->SetMarkerStyle(21);
    Etot[0]->Draw("PE");
    Etot[1]->SetLineColor(9);
    Etot[1]->SetMarkerColor(9);
    Etot[1]->SetMarkerStyle(22);
    Etot[1]->Draw("PEsame");
    myc->Update();
    myc->Print("PLOTS/TotalHitEnergy_"+genEnStr+"GeV.png");
    myc->Print("PLOTS/TotalHitEnergy_"+genEnStr+"GeV.pdf");

  }//loop on genEn
  
  
  return 0;
  
}//main
