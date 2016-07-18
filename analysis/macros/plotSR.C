#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>

#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TString.h"
#include "TLatex.h"
#include "TGaxis.h"

#include "TDRStyle.h"

int plotSR(){

  SetTdrStyle();

  const unsigned nLayers = 30;

  TCanvas *myc = new TCanvas("myc","",1);
  gStyle->SetOptStat(0);

  TH2F *h5 = new TH2F("h5","SR 5;width (pads);layer",
		     5,-2.5,2.5,nLayers,0,nLayers);
  TH2F *h6 = new TH2F("h6","SR 6;width (pads);layer",
		     5,-2.5,2.5,nLayers,0,nLayers);

  const unsigned nSR=8;
  double noise = 0.12;

  TF1 *fgaus = new TF1("fgaus","TMath::Gaus(x,[0],[1],0)",-2,2);
  fgaus->SetParameters(0,0.12);
  double proba = fgaus->Integral(0.5,2);

  std::cout << " -- Proba to pass 0.5 MIP cut just from gaussian noise : " << proba << std::endl;

  // myc->cd();
  //fgaus->Draw();

  //return 1;

  unsigned nCells[nSR] = {1,4,9,16,25,400,518,140*140};
  double calib = 268.;
  TLine *line = new TLine(0.2,0,0.2,nLayers);
      
  TLatex lat;

  myc->Clear();
  myc->Divide(2,1);

  for (unsigned iSR(0); iSR<nSR;++iSR){

    if (iSR==5){
      for (unsigned iL(0);iL<nLayers;++iL){	      
	if (iL==0) h5->Fill(0.,iL);
	else if (iL<5) h5->Fill(0.,iL);
	else if (iL<10) {h5->Fill(0.,iL);h5->Fill(1.,iL);}
	else if (iL<15) {h5->Fill(-1.,iL);h5->Fill(0.,iL);h5->Fill(1.,iL);}
	else if (iL<20) {h5->Fill(-1.,iL);h5->Fill(0.,iL);h5->Fill(1.,iL);h5->Fill(2.,iL);}
	else {h5->Fill(-2.,iL);h5->Fill(-1.,iL);h5->Fill(0.,iL);h5->Fill(1.,iL);h5->Fill(2.,iL);}
      }
      
      myc->cd(1);
      h5->Draw("col");
      
      line->Draw();
      lat.DrawLatex(0.3,15,"photon");
      lat.DrawLatexNDC(0.4,0.96,"SR 5");

    }
    else if (iSR == 6){
      for (unsigned iL(0);iL<nLayers;++iL){	      
	if (iL==0) h6->Fill(0.,iL);
	else if (iL<5) {h6->Fill(0.,iL);}
	else if (iL<12) {h6->Fill(-1.,iL);h6->Fill(0.,iL);h6->Fill(1.,iL);}
	else {h6->Fill(-2.,iL);h6->Fill(-1.,iL);h6->Fill(0.,iL);h6->Fill(1.,iL);h6->Fill(2.,iL);}
      }
      
      myc->cd(2);
      h6->Draw("col");
      line->Draw();
      lat.DrawLatex(0.3,15,"photon");
      lat.DrawLatexNDC(0.4,0.96,"SR 6");

      myc->Print("PLOTS/SignalRegions5and6.pdf");

    }
    else {
      nCells[iSR] *= nLayers;
    }
    double noiseMIP = noise*nCells[iSR]*proba;
    double noiseGeV = noiseMIP/calib;

    std::cout << " SR " << iSR << " " << nCells[iSR] << " " << noiseMIP << " MIPs, " << noiseGeV << " GeV." << std::endl;
    
    
  }
  
  return 0;
}//main
