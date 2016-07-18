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
#include "TRandom3.h"
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

int toyStoN2() {

  SetTdrStyle();
  TCanvas *myc = new TCanvas("myc","myc",1500,1000);

  const double noise = 1.;
  const unsigned nEvts = 100000;
  const unsigned nPts = 10;
  TRandom3 lRndm;

  for (unsigned ip(0);ip<nPts+1;++ip){
    const unsigned nEvtsNoise = ip*nEvts/nPts;
    
    TH1F* mipSignal = new TH1F("mipSignal",";E (mips);cells",100,0,5);
    TH1F* mip = new TH1F("mip",";E (mips);cells",100,0,5);
    TH1F* mipLangaus = new TH1F("mipLangaus",";E (mips);cells",100,0,5);
    TH1F* mipNoise = new TH1F("mipNoise",";E (mips);cells",100,0,5);
    TF1 *landaugaus = new TF1("landaugaus",langaufungaus,0,5,7);
    Double_t mpshift  = 0;//-0.22278298;       // Landau maximum location
    
    for (unsigned ie(0); ie<nEvts;++ie){
      double landauMean = 1;
      double lSignal = lRndm.Landau(landauMean,0.12);
      mipSignal->Fill(lSignal);
      double lNoise = lRndm.Gaus(0,noise);
      mip->Fill(lSignal+lNoise);
      mipLangaus->Fill(lSignal+lNoise);
    }
    for (unsigned ie(0); ie<nEvtsNoise;++ie){
      double lNoise = lRndm.Gaus(0,noise);
      mip->Fill(lNoise);
      mipNoise->Fill(lNoise);
    }
    
    double lmax = std::max(mipSignal->GetMaximum(),mip->GetMaximum());
    
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    myc->cd();
    mipSignal->SetMaximum(lmax);
    mipSignal->SetLineColor(2);
    //mipSignal->SetMarkerStyle(22);
    //mipSignal->SetMarkerColor(2);
    mipSignal->Draw("");
    mipSignal->Fit("landau","0");
    TF1 *fit = (TF1*)mipSignal->GetFunction("landau");
    fit->SetLineColor(2);
    fit->Draw("same");
    // MP shift correction
    double mpc = fit->GetParameter(1);// - mpshift * fit->GetParameter(2); 
    double mpcerr = fit->GetParError(1);//sqrt(pow(fit->GetParError(1),2)+pow(mpshift*fit->GetParError(2),2));
    //landauCor->SetParameters(fit->GetParameter(2),fit->GetParameter(1),fit->GetParameter(0),0);
    //mipSignal->Fit("landauCor","+","same",0.8,1.6);
    mip->SetLineColor(4);
    //mip->SetMarkerStyle(20);
    //mip->SetMarkerColor(1);
    mip->Draw("same");
    landaugaus->SetParameters(fit->GetParameter(2),mpc,fit->GetParameter(0),noise,nEvtsNoise,0,noise);
    landaugaus->FixParameter(3,noise);
    landaugaus->FixParameter(5,0);
    landaugaus->FixParameter(6,noise);
    landaugaus->SetLineColor(4);
    mip->Fit("landaugaus","B+","same");
    
    mipNoise->SetLineColor(3);
    mipNoise->Draw("same");
    mipLangaus->SetLineColor(6);
    mipLangaus->Draw("same");

    TLatex lat;
    char buf[500];
    sprintf(buf,"signal MPV = %3.3f #pm %3.3f",mpc,mpcerr);
    lat.DrawLatexNDC(0.4,0.85,buf);
    sprintf(buf,"signal+noise MPV = %3.3f #pm %3.3f",landaugaus->GetParameter(1),landaugaus->GetParError(1));
    lat.DrawLatexNDC(0.4,0.75,buf);
    sprintf(buf,"purity = %3.2f",1.0*nEvts/(nEvts+nEvtsNoise));
    lat.DrawLatexNDC(0.4,0.6,buf);

    TLegend *leg = new TLegend(0.55,0.16,0.94,0.5);
    leg->SetFillColor(10);
    leg->AddEntry(mipSignal,"Signal","L");
    leg->AddEntry(mipNoise,"Noise","L");
    leg->AddEntry(mipLangaus,"Signal#otimesNoise","L");
    leg->AddEntry(mip,"Noise+Signal#otimesNoise","L");
    leg->Draw("same");

    std::ostringstream label;
    label << "PLOTS/ToyStoN2_purity" << std::setprecision(2) << 1.0*nEvts/(nEvts+nEvtsNoise) << "_noise" << noise << ".pdf";
    myc->Print(label.str().c_str());
  }

  return 0;

}//main
