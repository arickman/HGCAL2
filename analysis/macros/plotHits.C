#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>
#include<cmath>

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

int plotHits() {

  //TFile *in1 = TFile::Open("../PLOTS/gitV00-02-03/version23/pi-/validation_20.root");
  //TFile *in2 = TFile::Open("../PLOTS/gitV00-02-03/version23/pi-/validation_20_1cm.root");
  //TFile *in3 = TFile::Open("../PLOTS/gitV00-02-03/version23/pi-/validation_20_3cm_noXTnonoise.root");
  TFile *in3 = TFile::Open("../PLOTS/gitV00-02-03/version23/pi-/validation_e30.root");
  in3->cd();

  TCanvas *myc = new TCanvas("myc","",1);
  myc->cd();
  myc->SetGridx(1);
  myc->SetGridy(1);
  gStyle->SetOptStat(0);

  TH1F *simnoise = (TH1F*)gDirectory->Get("p_simnoisehitEnergy");
  TH1F *xtalknoise = (TH1F*)gDirectory->Get("p_xtalknoisehitEnergy");
  TH1F *sim = (TH1F*)gDirectory->Get("p_simhitEnergy");
  TH1F *rec = (TH1F*)gDirectory->Get("p_rechitEnergy");
  sim->Sumw2();
  simnoise->Sumw2();
  xtalknoise->Sumw2();
  rec->Sumw2();

  rec->SetLineColor(4);
  rec->SetLineWidth(2);
  rec->SetMarkerColor(4);
  rec->SetMarkerStyle(2);
  sim->SetLineColor(2);
  sim->SetLineWidth(2);
  sim->SetMarkerColor(2);
  sim->SetMarkerStyle(3);
  //sim->Rebin(5);
  simnoise->SetLineColor(kViolet);
  simnoise->SetLineWidth(2);
  simnoise->SetMarkerColor(kViolet);
  simnoise->SetMarkerStyle(1);
  xtalknoise->SetLineColor(kGreen);
  xtalknoise->SetLineWidth(2);
  xtalknoise->SetMarkerColor(kGreen);
  xtalknoise->SetMarkerStyle(20);
  //simnoise->Rebin(5);
  //rec->Rebin(5);

  int binMin = sim->FindBin(0.5);
  int binMax = sim->GetNbinsX()+1;//FindBin(4);

  std::cout << " -- entries: " << sim->GetEntries() << " " 
	    << simnoise->GetEntries() << " "
	    << xtalknoise->GetEntries() << " "
	    << rec->GetEntries() << std::endl;
 std::cout << " -- integral 0.5-4: " << sim->Integral(binMin,binMax) << " " 
	    << simnoise->Integral(binMin,binMax) << " "
	    << xtalknoise->Integral(binMin,binMax) << " "
	    << rec->Integral(binMin,binMax) << std::endl;


  sim->Scale(1./sim->Integral(binMin,binMax));
  simnoise->Scale(1./simnoise->Integral(binMin,binMax));
  xtalknoise->Scale(1./xtalknoise->Integral(binMin,binMax));
  rec->Scale(1./rec->Integral(binMin,binMax));

  sim->GetXaxis()->SetRangeUser(0.,4);
  simnoise->GetXaxis()->SetRangeUser(0.,4);
  xtalknoise->GetXaxis()->SetRangeUser(0.,4);
  rec->GetXaxis()->SetRangeUser(0.,4);

  sim->GetYaxis()->SetTitle("Normalised hit energy / 0.05 MIP");


  //sim->Scale(505278./sim->GetEntries());
  //simnoise->Scale(505278./simnoise->GetEntries());
  //rec->Scale(375960./rec->GetEntries());
  sim->SetTitle("30 GeV pi-");
  //simnoise->GetYaxis()->SetRangeUser(0,24000);

  sim->Draw("hist");
  sim->GetYaxis()->SetRangeUser(0,0.1);
  simnoise->Draw("histsame");
  xtalknoise->Draw("histsame");
  rec->Draw("histsame");
  //sim->Draw("same");

  TLegend *leg = new TLegend(0.4,0.6,0.89,0.89);
  leg->SetFillColor(10);
  leg->AddEntry(sim,"SimHits","P");
  leg->AddEntry(simnoise,"SimHits + noise","L");
  leg->AddEntry(xtalknoise,"SimHits + xtalk + noise","L");
  leg->AddEntry(rec,"RecHits","L");
  leg->Draw("same");

  TLatex lat;
  lat.SetTextColor(2);
  char buf[500];
  sprintf(buf,"<E>=%3.3f #pm %3.3f MIPs",(sim->GetMean()),(sim->GetMeanError()));
  //lat.DrawLatexNDC(0.15,0.8,buf);
  sprintf(buf,"#sigma(E)=%3.3f #pm %3.3f MIPs",(sim->GetRMS()),(sim->GetRMSError()));
  //lat.DrawLatexNDC(0.15,0.75,buf);
  lat.SetTextColor(4);
  sprintf(buf,"<E>=%3.3f #pm %3.3f MIPs",(rec->GetMean()),(rec->GetMeanError()));
  //lat.DrawLatexNDC(0.15,0.7,buf);
  sprintf(buf,"#sigma(E)=%3.3f #pm %3.3f MIPs",(rec->GetRMS()),(rec->GetRMSError()));
  //lat.DrawLatexNDC(0.15,0.65,buf);


  myc->Print("HitSpectrum_pi-_30.pdf");

  return 0;
}//main
