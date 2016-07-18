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

int plotMaxCharge(){//main

  SetTdrStyle();

  const unsigned nLayers = 30;

  const unsigned nPhi = 100;

  const unsigned npu = 2;
  const unsigned pu[npu] = {140,200};

  TFile *input[npu];
  input[0] = TFile::Open("../PLOTS/MaxChargeStudy_MB_pu140.root");
  input[1] = TFile::Open("../PLOTS/MaxChargeStudy_MB_pu200.root");

  std::ostringstream label;
  const unsigned nC = 10;
  TCanvas *myc[nC];
  for (unsigned ic(0); ic<nC;++ic){
    label.str("");
    label << "myc" << ic;
    myc[ic] = new TCanvas(label.str().c_str(),label.str().c_str(),1500,750);
    myc[ic]->Divide(2,1);
  }

  gStyle->SetOptStat("eMRuo");
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);
  gStyle->SetStatW(0.4);
  gStyle->SetStatH(0.4);

  std::string histname[nC] = {
    "hxy",
    "hphieta",
    "hxyAbove",
    "hphietaAbove",
    "nAbove",
    "nTot",
    "prof_nAbovevsLayer",
    "EmaxNeighbourSC3D",
    "hLayerProfile_above50mips",
    "EmaxNeighbourSC3DInt"
  };

  TH1 *hist[nC];

  TLatex lat;
  char buf[500];
  TH2F *nAbovevsLayer[npu];
  TH1F *EmaxNeighbourSC3DInt[npu];
  EmaxNeighbourSC3DInt[0] = new TH1F("EmaxNeighbourSC3DInt0",";E_{max}^{SC3D} (Mips);(hits E<E_{max}^{SC3D})/(hits in 1.9 < #eta < 2.1)",400,0,100);
  EmaxNeighbourSC3DInt[1] = new TH1F("EmaxNeighbourSC3DInt1",";E_{max}^{SC3D} (Mips);(hits E<E_{max}^{SC3D})/(hits in 1.9 < #eta < 2.1)",400,0,100);


  for (unsigned ih(0); ih<nC; ++ih){
    for (unsigned iF(0); iF<npu;++iF){
      input[iF]->cd();
      myc[ih]->cd(iF+1);
      if (ih<4) {
	gPad->SetRightMargin(0.2);
	gStyle->SetStatX(0.78);
	gStyle->SetStatY(0.94);
	gStyle->SetStatW(0.4);
	if (ih%2==0) {
	  gStyle->SetStatH(0.05);
	  gStyle->SetOptStat("e");
	}
	else {
	  gStyle->SetStatH(0.2);
	  gStyle->SetOptStat("eMR");
	}
	hist[ih] = (TH2F*)gDirectory->Get(histname[ih].c_str());
	nAbovevsLayer[iF] = (TH2F*)gDirectory->Get("nAbovevsLayer");
	hist[ih]->Draw("colz");
 	sprintf(buf,"Pythia Minbias Pu %d",pu[iF]);
	lat.DrawLatexNDC(0.3,0.96,buf);
      }
      else {
	gStyle->SetStatX(0.9);
	gStyle->SetStatY(0.9);
	gStyle->SetStatW(0.4);
	gStyle->SetStatH(0.4);
	gStyle->SetOptStat("eMRuo");
	if (ih< (nC-1)) hist[ih] = (TH1F*)gDirectory->Get(histname[ih].c_str());
	else hist[ih] = EmaxNeighbourSC3DInt[iF];
	if (!hist[ih]) return 1;

	hist[ih]->SetMarkerStyle(20);
	hist[ih]->SetMarkerColor(1);
	hist[ih]->SetLineColor(1);
	unsigned total = 0;
	std::cout << " -- n in bin 0 = " << hist[ih]->GetBinContent(1) << std::endl;
	for (int b(1);b<hist[ih]->GetNbinsX()+1;++b){
	  total += hist[ih]->GetBinContent(b)*hist[ih]->GetXaxis()->GetBinLowEdge(b);
	}
	std::cout << " File " << iF << " hist " << histname[ih] << " nTotalCells = " << total << std::endl;
	if (ih==4) {
	  hist[ih]->GetXaxis()->SetRangeUser(1,50);
	  hist[ih]->GetYaxis()->SetRangeUser(1,8000);

	  gPad->SetLogy(1);
	  gPad->SetGridx(1);
	  gPad->SetGridy(1);
	}
	else if (ih<6) hist[ih]->GetXaxis()->SetRangeUser(0,130);
	else if (ih==6) {
	  gPad->SetLogy(1);
	  gStyle->SetOptStat(0);
	  hist[ih]->GetYaxis()->SetRangeUser(0.000001,0.001);
	}
	else if (ih==7) {
	  gPad->SetLogy(1);
	  for (int b(1);b<hist[ih]->GetNbinsX()+1;++b){
	    double error = 0;
	    double integral = hist[ih]->IntegralAndError(0,b,error);
	    EmaxNeighbourSC3DInt[iF]->SetBinContent(b,integral/hist[ih]->Integral(0,hist[ih]->GetNbinsX()+2));
	    EmaxNeighbourSC3DInt[iF]->SetBinError(b,error/hist[ih]->Integral(0,hist[ih]->GetNbinsX()+2));
	  }
	  std::cout << " -- Check integral : " << hist[ih]->Integral(0,hist[ih]->GetNbinsX()+2) << std::endl;
	}
	else if (ih==8) gStyle->SetOptStat("e");
	else if (ih==9) {
	  gPad->SetGridx(1);
	  gPad->SetGridy(1);
	  hist[ih]->GetYaxis()->SetRangeUser(0,1);
	  gStyle->SetOptStat(0);
	}
	hist[ih]->Draw("PE");
	gPad->Modified();
	gPad->Update();
	if (ih==9) {
	  gStyle->SetOptFit(0);//1111);
	  TF1 *f1 = new TF1("f1","[0]+[1]*x",0,100);
	  hist[ih]->Fit("f1","S","0",0,5);
	  TF1 *f2 = new TF1("f2","[0]+[1]*x",0,100);
	  hist[ih]->Fit("f2","S","0",40,100);
	  //f1 = (TF1*)hist[ih]->GetFunction("f1");
	  f1->SetLineColor(4);
	  f1->SetLineWidth(2);
	  f1->Draw("same");
	  lat.SetTextColor(4);
	  lat.DrawLatexNDC(0.2,0.6,"HIPs");
	  gPad->Modified();
	  gPad->Update();
	  //f2 = (TF1*)hist[ih]->GetFunction("f2");
	  f2->SetLineColor(6);
	  f2->SetLineWidth(2);
	  f2->Draw("same");
	  lat.SetTextColor(6);
	  lat.DrawLatexNDC(0.7,0.8,"Showers");
	  gPad->Modified();
	  gPad->Update();
	  lat.SetTextColor(1);
	}
	//if (ih==6) nAbovevsLayer[iF]->Draw("colzsame");
	sprintf(buf,"Pythia Minbias Pu %d",pu[iF]);
	lat.DrawLatexNDC(0.3,0.96,buf);
	if (ih==4){
	  sprintf(buf,"n_{tot}(E>250fC)=%d",total);
	  lat.DrawLatexNDC(0.4,0.3,buf);
	  double rate = total/(1000*nPhi*270.);
	  unsigned p=0;
	  for (;p<10;++p){
	    std::cout << p << " " << rate*pow(10,p) << " " << pow(10,p) << std::endl;
	    if (rate*pow(10,p)>=1) break;
	  }
	  unsigned denom = pow(10,p);
	  sprintf(buf,"rate = %3.1f/%d",rate*pow(10,p),denom);
	  lat.DrawLatexNDC(0.4,0.2,buf);
	}
      }
      myc[ih]->Update();
      gPad->Modified();
      gPad->Update();
    }//loop on files
    myc[ih]->Update();
    label.str("");
    label << "PLOTS/" << histname[ih] << "_MB.pdf";
    myc[ih]->Print(label.str().c_str());

  }//loop on hists





  return 0;
}//main
