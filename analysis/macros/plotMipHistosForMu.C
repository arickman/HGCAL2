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
#include "TMath.h"
#include "TVectorD.h"
#include "Math/Vector3D.h"
#include "Math/Vector3Dfwd.h"
#include "Math/Point2D.h"
#include "Math/Point2Dfwd.h"
#include "TSystem.h"

#include "TDRStyle.h"


Double_t langaufun(Double_t *x, Double_t *par) {

   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation), 
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

      // Numeric constants
      Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      Double_t mpshift  = 0;//-0.22278298;       // Landau maximum location

      // Control constants
      Double_t np = 100.0;      // number of convolution steps
      Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      Double_t xx;
      Double_t mpc;
      Double_t fland;
      Double_t sum = 0.0;
      Double_t xlow,xupp;
      Double_t step;
      Double_t i;


      // MP shift correction
      mpc = par[1] - mpshift * par[0]; 

      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];

      step = (xupp-xlow) / np;

      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);

         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
      }

      return (par[2] * step * sum * invsq2pi / par[3]);
};


Double_t langaufungaus(Double_t *x, Double_t *par) {
  
  return langaufun(x,par)+par[4]*TMath::Gaus(x[0],par[5],par[6],true);

};

/*bool getZpositions(const unsigned versionNumber, std::vector<double> & posz){
  std::ifstream fin;
  std::ostringstream finname;
  //finname << outFolder_ << "/zPositions.dat";
  finname << "../data/zPositions_v" << versionNumber << ".dat";
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
  
  }*/

int plotMipHistosForMu(){

  ///////////////////////////////////////////////////////////////////
  ////////// Hardcoded config ///////////////////////////////////////
  ///////////////////////////////////////////////////////////////////

  const double eta = 2.8;
  const double deta = 0.15;

  const unsigned nLayers = 30;

  ///////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////

  //std::vector<double> posz;
  //posz.resize(nLayers,0);
  //if (!getZpositions(12,posz)) return 1;


  std::ostringstream outFilePath;
  outFilePath << "/afs/cern.ch/work/a/amagnan/PFCalEEAna/PLOTS/gitV00-02-12/version12/MinBias/MipStudy_eta" << eta-deta << "_" << eta+deta << "_MinBias5_5.root";
  TFile *file = TFile::Open(outFilePath.str().c_str());
  file->cd();

  TTree *tree = (TTree*)gDirectory->Get("RecoTree"); 

  //TF1 *landaugaus = new TF1("landaugaus",langaufungaus,0,5,7);

  double mpshift  = 0;//-0.22278298;

  TH1F *hitSpectrum[nLayers];

  std::ostringstream label;

  SetTdrStyle();
  TCanvas *myc = new TCanvas("myc","myc",2000,1000);
  TCanvas *mycL = new TCanvas("mycL","mycL",2000,1000);
  mycL->Divide(5,2);
  TCanvas *mycS = new TCanvas("mycS","mycS",2000,1000);
  TCanvas *mycM = new TCanvas("mycM","mycM",2000,1000);


  // TCanvas *mycL[nLayers];
  // for (unsigned il(0); il<nLayers;++il){
  //   label.str("");
  //   label << "mycL" << il ;
  //   mycL[il] = new TCanvas(label.str().c_str(),label.str().c_str(),1500,1000);
  //   mycL[il]->Divide(3,2);
  // }


  gStyle->SetOptStat("eMRuo");
  gStyle->SetOptFit(1111);
  TLatex lat;
  char buf[500];


  myc->cd();

  //TH1F *hitSpectrum[nEta][nNoise][nLayers];
  TH1F *hitSpectrumAll = new TH1F("hitSpectrumAll",";E_{1#times1cm^{2}}(mips);N_{cells}",100,0,5);
  std::ostringstream lvar;
  lvar << "2*HGCSSRecoHitVec.energy_>>hitSpectrumAll";
  //lvar << "2*HGCSSRecoHitVec.energy_/(fabs(posz[HGCSSRecoHitVec.layer_])/sqrt(posz[HGCSSRecoHitVec.layer_]*posz[HGCSSRecoHitVec.layer_]+HGCSSRecoHitVec.xpos_*HGCSSRecoHitVec.xpos_+HGCSSRecoHitVec.ypos_*HGCSSRecoHitVec.ypos_))>>hitSpectrumAll";

  tree->Draw(lvar.str().c_str(),"(HGCSSRecoHitVec.layer_>=0)");

  TGraphErrors *grMPV = new TGraphErrors();
  grMPV->SetName("grMPV");
  TGraphErrors *grSigma = new TGraphErrors();
  grSigma->SetName("grSigma");


  for (unsigned il(0); il<nLayers;++il){
    label.str("");
    label << "hitSpectrum_layer" << il;
    hitSpectrum[il] = new TH1F(label.str().c_str(),";E_{1#times1cm^{2}}(mips);N_{cells}",100,0,5);
    lvar.str("");
    lvar << "2*HGCSSRecoHitVec.energy_>>" << label.str();
    //lvar << "2*HGCSSRecoHitVec.energy_/(fabs(posz[HGCSSRecoHitVec.layer_])/sqrt(posz[HGCSSRecoHitVec.layer_]*posz[HGCSSRecoHitVec.layer_]+HGCSSRecoHitVec.xpos_*HGCSSRecoHitVec.xpos_+HGCSSRecoHitVec.ypos_*HGCSSRecoHitVec.ypos_))>>" << label.str();
    std::ostringstream lcut;
    lcut << "(HGCSSRecoHitVec.layer_==" << il << ")";
    tree->Draw(lvar.str().c_str(),lcut.str().c_str());
  }

  myc->Clear();
  

  for (unsigned il(0); il<nLayers;++il){
    mycL->cd(il/3+1);
    //hitSpectrum[il]->Rebin(2);
    hitSpectrum[il]->SetLineColor(1+il%3);
    hitSpectrum[il]->Draw(il%3==0?"":"same");
    if (hitSpectrum[il]->GetEntries()==0) continue;
    TF1 *fit = 0;
    double mpc = 0;
    double minfit = 0.7;
    double maxfit = 1.7;
    hitSpectrum[il]->Fit("landau","R+","same",minfit,maxfit);
    fit = (TF1*)hitSpectrum[il]->GetFunction("landau");
    if (!fit) continue;
    //landaugaus->SetParameters(fit->GetParameter(2),fit->GetParameter(1),fit->GetParameter(0),0,fit->GetParameter(0),0,0);
    mpc = fit->GetParameter(1) - mpshift * fit->GetParameter(2); 
    if (il%3==0){
      sprintf(buf,"Layer %d",il);
      lat.DrawLatexNDC(0.2,0.85,buf);
      sprintf(buf,"%3.2f < #eta < %3.2f",eta-deta,eta+deta);
      lat.DrawLatexNDC(0.2,0.75,buf);
    }
    else {
      lat.SetTextColor(il%3+1);
      sprintf(buf,"Layer %d",il);
      lat.DrawLatexNDC(0.35,0.75-0.1*(il%3),buf);
      lat.SetTextColor(1);
    }
    grMPV->SetPoint(il,il,mpc);
    grMPV->SetPointError(il,0,fit->GetParError(1));
    grSigma->SetPoint(il,il,fit->GetParameter(2));
    grSigma->SetPointError(il,0,fit->GetParError(2));

  }
  label.str("");
  label << "PLOTS/Muons/MipSpectrum_eta" << eta*10 ;
  //label << "_layer" << il;
  //label << suffix.str();
  label << "_perlayer.pdf";
  mycL->Print(label.str().c_str());

  gStyle->SetOptFit(1111);
  gStyle->SetOptStat("eMRuo");

  myc->cd();
  gPad->SetLogy(1);
  hitSpectrumAll->Draw();
  double minfit = 0.5;
  double maxfit = 1.5;
  minfit = 0.75;
  maxfit = 1.4;
  hitSpectrumAll->Fit("landau","R+","same",minfit,maxfit);
  TF1 *fit = (TF1*)hitSpectrumAll->GetFunction("landau");
  sprintf(buf,"%3.2f < #eta < %3.2f",eta-deta,eta+deta);
  lat.DrawLatexNDC(0.2,0.75,buf);

  myc->Update();
  label.str("");
  label << "PLOTS/Muons/MipSpectrum_eta" << eta*10 ;
  label << ".pdf";
  myc->Print(label.str().c_str());

  myc->cd();
  gPad->SetLogy(0);
  hitSpectrumAll->Draw();
  sprintf(buf,"%3.2f < #eta < %3.2f",eta-deta,eta+deta);
  lat.DrawLatexNDC(0.2,0.75,buf);

  myc->Update();
  label.str("");
  label << "PLOTS/Muons/MipSpectrum_eta" << eta*10 ;
  label << "_nolog.pdf";
  myc->Print(label.str().c_str());


  mycM->cd();
  gPad->SetGridy(1);
  grMPV->SetTitle(";layer;Landau MPV (MIPs)");
  grMPV->SetMarkerStyle(20);
  grMPV->SetMarkerColor(1);
  grMPV->SetLineColor(1);
  grMPV->SetMinimum(0.8);
  grMPV->SetMaximum(1.2);
  grMPV->Draw("APL");
  sprintf(buf,"#eta = %3.1f",eta);
  lat.DrawLatexNDC(0.2,0.85,buf);

  mycM->Update();
  label.str("");
  label << "PLOTS/Muons/LandauMPV_eta" << eta*10 ;
  //label << suffix.str();
  label << ".pdf";
  mycM->Print(label.str().c_str());

  mycS->cd();
  gPad->SetGridy(1);
  grSigma->SetTitle(";layer;Landau #sigma (MIPs)");
  grSigma->SetMarkerStyle(20);
  grSigma->SetMarkerColor(1);
  grSigma->SetLineColor(1);
  grSigma->SetMinimum(0.);
  grSigma->SetMaximum(0.3);
  grSigma->Draw("APL");
  sprintf(buf,"#eta = %3.1f",eta);
  lat.DrawLatexNDC(0.2,0.85,buf);

  mycS->Update();
  label.str("");
  label << "PLOTS/Muons/LandauSigma_eta" << eta*10 ;
  //label << suffix.str();
  label << ".pdf";
  mycS->Print(label.str().c_str());


  return 0;
}//main
