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
#include "TPaveStats.h"


int plotOOTpu(){//main


  SetTdrStyle();
  gStyle->SetPadTopMargin(0.07);

  //TFile *input = TFile::Open("../PLOTS/pT60eta19Photons_MB_pu140.root");
  TFile *input = TFile::Open("../PLOTS/HiggsPhotons_MBtest_pu140.root");
  input->cd();
  TCanvas *myc = new TCanvas("myc","myc",1500,1000);
  TCanvas *myc2 = new TCanvas("myc2","myc2",1500,1000);
  TCanvas *myc3 = new TCanvas("myc3","myc3",1500,1000);
  TCanvas *myc4 = new TCanvas("myc4","myc4",1500,1000);
  TCanvas *myc5 = new TCanvas("myc5","myc5",1500,1000);

  myc->cd();

  const unsigned nH = 3;
  TH1F *h[nH];

  std::string histname[nH] = {
    "ErecooverEtrue",
    "ErecooverEtrueOOTpast",
    "ErecooverEtrueOOT"
  };  

  std::string labelStr[nH] = {
    "no PU BX 0",
    "PU BX -1",
    "PU BX -1+1-10"
  };

  gStyle->SetOptStat(0);//"eMRo");
  gStyle->SetOptFit(0);//11);
  TPaveStats *st[nH];

  TF1 *f[nH];
  TLegend *leg = new TLegend(0.63,0.7,0.93,0.9);
  leg->SetFillStyle(0);

  unsigned lcolor[nH] = {1,4,2};

  for (unsigned ih(0);ih<nH;++ih){
    h[ih] = (TH1F*)gDirectory->Get(histname[ih].c_str());
    if (!h[ih]) continue;
    h[ih]->Rebin(2);
    h[ih]->GetXaxis()->SetRangeUser(0.9,1.2);
    h[ih]->SetLineColor(lcolor[ih]);
    h[ih]->SetMarkerColor(lcolor[ih]);
    if (ih>0) h[ih]->SetMarkerStyle(20+ih);
    h[ih]->Draw(ih==0?"":"PEsame");
    //h[ih]->Fit("gaus","","same",0.9,1.02);
    //f[ih] = (TF1*)h[ih]->GetFunction("gaus");
    //f[ih]->SetLineColor(lcolor[ih]);
    //f[ih]->Draw("same");
    gPad->Update();
    //st[ih] = (TPaveStats*)h[ih]->FindObject("stats");
    //if (!st[ih]) continue;
    //st[ih]->SetLineColor(lcolor[ih]);
    //st[ih]->SetTextColor(lcolor[ih]);
    //st[ih]->SetX1NDC(0.63);
    //st[ih]->SetX2NDC(0.93);
    //st[ih]->SetY1NDC(0.72-0.23*ih);
    //st[ih]->SetY2NDC(0.92-0.23*ih);
    //st[ih]->SetLabel(labelStr[ih].c_str());
    leg->AddEntry(h[ih],labelStr[ih].c_str(),ih==0?"L":"PL");
  }

  TLatex lat;
  //lat.DrawLatexNDC(0.3,0.96,"Single #gamma p_{T}=60 GeV, #eta=1.9");
  lat.DrawLatexNDC(0.3,0.94,"Pythia H#rightarrow#gamma#gamma, 1.5 < #eta^{#gamma} < 2.9");
  lat.DrawLatexNDC(0.02,0.02,"HGCAL Standalone Simulation");

  leg->Draw("same");

  myc->Update();
  myc->Print("OOTPUImpact_Hgg.pdf");

  TLegend *leg2 = new TLegend(0.13,0.72,0.43,0.92);
  leg2->SetFillStyle(0);
  labelStr[0] = "PU BX -1";
  labelStr[1] = "PU BX 1-10";
  labelStr[2] = "PU BX -1+1-10";

  myc4->cd();
  gPad->SetLogy(1);
  histname[0] = "EpuPast";
  histname[1] = "EpuFutur";
  histname[2] = "EpuAll";
  for (unsigned ih(0);ih<nH;++ih){
    h[ih] = (TH1F*)gDirectory->Get(histname[ih].c_str());
    if (!h[ih]) continue;
    h[ih]->Rebin(10);
    h[ih]->GetXaxis()->SetRangeUser(0,1000);
    h[ih]->GetYaxis()->SetRangeUser(0.3,50000);
    h[ih]->SetLineColor(lcolor[ih]);
    h[ih]->SetMarkerColor(lcolor[ih]);
    h[ih]->SetMarkerStyle(20+ih);
    h[ih]->Draw(ih==0?"PE":"PEsames");
    gPad->Update();
    st[ih] = (TPaveStats*)h[ih]->FindObject("stats");
    if (!st[ih]) continue;
    st[ih]->SetLineColor(lcolor[ih]);
    st[ih]->SetTextColor(lcolor[ih]);
    st[ih]->SetX1NDC(0.63);
    st[ih]->SetX2NDC(0.93);
    st[ih]->SetY1NDC(0.72-0.21*ih);
    st[ih]->SetY2NDC(0.92-0.21*ih);
    //st[ih]->SetLabel(labelStr[ih].c_str());
    leg2->AddEntry(h[ih],labelStr[ih].c_str(),"PL");
  }

  //lat.DrawLatexNDC(0.3,0.96,"Single #gamma p_{T}=60 GeV, #eta=1.9");
  lat.DrawLatexNDC(0.3,0.94,"Pythia H#rightarrow#gamma#gamma, 1.5 < #eta^{#gamma} < 2.9");
  lat.DrawLatexNDC(0.02,0.02,"HGCAL Standalone Simulation");

  leg2->Draw("same");

  myc4->Update();
  myc4->Print("OOTPUEnergy_Hgg.pdf");

  myc5->cd();
  gPad->SetLogy(1);
  histname[0] = "EpuoverErecoPast";
  histname[1] = "EpuoverErecoFutur";
  histname[2] = "EpuoverErecoAll";
  for (unsigned ih(0);ih<nH;++ih){
    h[ih] = (TH1F*)gDirectory->Get(histname[ih].c_str());
    if (!h[ih]) continue;
    h[ih]->Rebin(5);
    h[ih]->GetXaxis()->SetRangeUser(0.,0.05);
    h[ih]->SetLineColor(lcolor[ih]);
    h[ih]->SetMarkerColor(lcolor[ih]);
    h[ih]->SetMarkerStyle(20+ih);
    h[ih]->Draw(ih==0?"PE":"PEsames");
    gPad->Update();
    st[ih] = (TPaveStats*)h[ih]->FindObject("stats");
    if (!st[ih]) continue;
    st[ih]->SetLineColor(lcolor[ih]);
    st[ih]->SetTextColor(lcolor[ih]);
    st[ih]->SetX1NDC(0.63);
    st[ih]->SetX2NDC(0.93);
    st[ih]->SetY1NDC(0.72-0.21*ih);
    st[ih]->SetY2NDC(0.92-0.21*ih);
  }

  //lat.DrawLatexNDC(0.3,0.96,"Single #gamma p_{T}=60 GeV, #eta=1.9");
  lat.DrawLatexNDC(0.3,0.94,"Pythia H#rightarrow#gamma#gamma, 1.5 < #eta^{#gamma} < 2.9");
  lat.DrawLatexNDC(0.02,0.02,"HGCAL Standalone Simulation");

  leg2->Draw("same");

  myc5->Update();
  myc5->Print("OOTPUEOverPhotonE_Hgg.pdf");


  myc2->Clear();
  myc2->Divide(5,2);
  myc3->Divide(5,2);
  const unsigned nbx = 11;
  TH1F *Epu_perbx[nbx];
  TH1F *nAbove_perbx[nbx];

  gStyle->SetOptStat("eMRuo");
  gStyle->SetStatW(0.6);
  gStyle->SetStatH(0.4);

  for (unsigned ibx(1); ibx<nbx;++ibx){
    std::ostringstream label;
    label.str("");
    label << "Epu_bx_" << ibx;
    Epu_perbx[ibx] = (TH1F*)gDirectory->Get(label.str().c_str());
    label.str("");
    label << "nAbove_bx_" << ibx;
    nAbove_perbx[ibx] = (TH1F*)gDirectory->Get(label.str().c_str());

    myc2->cd(ibx);
    gPad->SetLogy(1);
    Epu_perbx[ibx]->GetXaxis()->SetRangeUser(0,200);
    Epu_perbx[ibx]->Draw();

    label.str("");
    label <<  "BX " << ibx;
    lat.DrawLatexNDC(0.5,0.96,label.str().c_str());

    myc3->cd(ibx);
    gPad->SetLogy(1);
    nAbove_perbx[ibx]->GetXaxis()->SetRangeUser(0,200);
    if (ibx>3) nAbove_perbx[ibx]->GetXaxis()->SetRangeUser(0,50);
    if (ibx>8) nAbove_perbx[ibx]->GetXaxis()->SetRangeUser(0,20);
    nAbove_perbx[ibx]->Draw();
    lat.DrawLatexNDC(0.5,0.96,label.str().c_str());

  }

  myc2->Update();
  myc2->Print("futurOOTPUenergy.pdf");
  myc3->Update();
  myc3->Print("futurOOTPUnHitsAbove.pdf");

  return 0;
}//main
