#include<string>
#include<iostream>
#include<fstream>
#include<sstream>

#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TString.h"
#include "TLatex.h"
#include "TPaletteAxis.h"

double xLayer(const unsigned layer) {
  double w1b = 1.6+3+0.1;
  double w2b = 3.3+3+0.1;
  double w3b = 5.6+3+0.1;
  double wa = 0.1+2;
  double W1 = 10*(w1b+wa);
  double W2 = 10*(w2b+wa);
  double W3 = 10*(w3b+wa);

  double totalWidth = 10*(w1b+w2b+w3b)+30*wa;
  unsigned il = layer%10+1;

  double offset = layer/10==0 ? (il*w1b + (il-1)*wa) :
    ( layer/10==1 ? W1 + il*w2b + (il-1)*wa :
      W1 + W2 + il*w3b + (il-1)*wa );

  return -totalWidth/2. + offset; //in mm
  
}

double getEta(double ypos, std::string etaStr, const double xLaymm){
  double eta0 = 0;
  if (etaStr.find("20")!=etaStr.npos) eta0 = 2.0;
  else if (etaStr.find("25")!=etaStr.npos) eta0 = 2.5;
  else if (etaStr.find("30")!=etaStr.npos) eta0 = 3.0;
  else if (etaStr.find("35")!=etaStr.npos) eta0 = 3.5;
  else {
    //std::cout << " ERROR: new eta point ! Please implement value to convert y position... Setting it to 10..." << std::endl;
     //exit(1);
     eta0=10;
  }
  //find height using center of detector
  double xDist0 = 3100;
  double theta0 = 2*atan(exp(-eta0));
  double y0 = tan(theta0)*xDist0;//mm
  //find theta at center of layer.
  double xDist = 3100+xLaymm;//distance from beam spot to layer in mm
  double tanTheta = (y0+ypos)/xDist;
  double theta = atan(tanTheta);
  double eta = -log(tan(theta/2.));
  //std::cout << " - etaBin " << etaStr << " y0 = " << y0 << " y = " << ypos << " eta = " << eta << std::endl;

  return eta;

}

int plotOcc(){//main  

  //const unsigned nS = 7;
  //TString scenario[nS] = {"0","1","2","3","4","5","6"};
  const unsigned nS = 4;
  std::string scenario[nS] = {
    "scenario_0/PedroPU/eta20/",
    "scenario_0/PedroPU/eta25/",
    "scenario_0/PedroPU/eta30/",
    "scenario_0/PedroPU/eta35/"
  };

  bool doOnlySummary = true;

  TString etaStr[nS] = {
    "#eta = 2.0",
    "#eta = 2.5",
    "#eta = 3.0",
    "#eta = 3.5",
  };

  const double eta[nS] = {2.0,2.5,3.0,3.5};


  const unsigned nV = 1;
  TString version[nV] = {"20"};
  const double Emip = 0.0548;//in MeV

  const unsigned nLayers = 30;
  const unsigned nEcalLayers = 30;

  const unsigned nOcc = 4;
  const unsigned occThreshold[nOcc] = {1,5,10,20};

  double xLay[nLayers];
  for (unsigned iL(0); iL<nLayers; ++iL){
    xLay[iL] = xLayer(iL);
  }

  const unsigned nPlots = 5;
  const unsigned nPads = 6;
  const unsigned nCanvas = nPlots*nS;  
  TCanvas *myc[nCanvas];
  for (unsigned iC(0);iC<nCanvas;++iC){
    std::ostringstream lName;
    lName << "myc" << iC;
    myc[iC] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1500,1000);
    myc[iC]->Divide(3,2);
  }
  
  std::ostringstream saveName;

  for (unsigned iV(0); iV<nV;++iV){//loop on versions
    for (unsigned iS(0); iS<nS;++iS){//loop on scenarios
    
      TString plotDir = "../PLOTS/version_"+version[iV]+"/"+scenario[iS]+"/";
      std::ostringstream lName;
      lName << plotDir << "CalibHistos.root";
      TFile *inputFile = TFile::Open(lName.str().c_str());
      if (!inputFile) {
	std::cout << " -- Error, input file " << lName.str() << " cannot be opened. Exiting..." << std::endl;
	return 1;
      }

      double maxOcc[nOcc][nLayers];
      double maxOccErr[nOcc][nLayers];

      double lay[nLayers];
      double layerr[nLayers];
      for (unsigned iL(0); iL<nLayers; ++iL){
	lay[iL] = iL;
	layerr[iL] = 0;
      }
      double Emax[nOcc];
      for (unsigned iO(0); iO<nOcc; ++iO){

	TH2F *p_occupancy[nLayers];
	for (unsigned iL(0); iL<nLayers; ++iL){
	  maxOcc[iO][iL] = 0;
	  maxOccErr[iO][iL]  = 0;
	}

	Emax[iO] = 0;
	unsigned counter = 0;

	for (unsigned iL(0); iL<nLayers; ++iL){
	  lName.str("");
	  lName << "p_occupancy_" << occThreshold[iO] << "_" << iL;
	  p_occupancy[iL] = (TH2F*)gDirectory->Get(lName.str().c_str());
	  if (!p_occupancy[iL]) {
	    std::cout << " -- ERROR, pointer for histogram is null for layer: " << iL << ". Exiting..." << std::endl;
	    return 1;
	  }
	  p_occupancy[iL]->Sumw2();
	  p_occupancy[iL]->Scale(100./1000.);

	  int binNominal = 1;
	  double minDiff = 100;
	  for (int iB(1); iB<p_occupancy[iL]->GetNbinsY()+1;++iB){
	    char buf[500];
	    double ypos = p_occupancy[iL]->GetYaxis()->GetBinCenter(iB);
	    double lEta = getEta(ypos,scenario[iS],xLay[iL]);
	    sprintf(buf,"#eta %1.2f",lEta);
	    p_occupancy[iL]->GetYaxis()->SetBinLabel(iB,buf);
	    //find eta closest to nominal
	    if (fabs(lEta-eta[iS])<minDiff) {
	      minDiff = fabs(lEta-eta[iS]);
	      binNominal = iB;
	    }
	  }

	  int bin0Xp = p_occupancy[iL]->GetXaxis()->FindBin(0.5);
	  int bin0Xn = p_occupancy[iL]->GetXaxis()->FindBin(-0.5);
	  double content0p = p_occupancy[iL]->GetBinContent(bin0Xp,binNominal);
	  double content0n = p_occupancy[iL]->GetBinContent(bin0Xn,binNominal);
	  int bin0X = content0p>content0n ? bin0Xp : bin0Xn;
	  maxOcc[iO][iL] = p_occupancy[iL]->GetBinContent(bin0X,binNominal);
	  maxOccErr[iO][iL] = p_occupancy[iL]->GetBinError(bin0X,binNominal);
	  double ypos = p_occupancy[iL]->GetYaxis()->GetBinCenter(binNominal);
	  double lEta = getEta(ypos,scenario[iS],xLay[iL]);
	  std::cout << iO << " " << iL << " " 
		    << p_occupancy[iL]->GetXaxis()->GetBinCenter(bin0X) << " "
		    << ypos << " " << lEta << " " 
		    << maxOcc[iO][iL] << " +/- " << maxOccErr[iO][iL] << std::endl;
	  
	  double Etot = p_occupancy[iL]->GetMaximum();
	  if (Etot > Emax[iO]) Emax[iO] = Etot;
	}//loop on layers

	if (doOnlySummary) continue;

	gStyle->SetOptStat(0);
	for (unsigned iL(0); iL<nLayers; ++iL){
	  myc[nPlots*iS+counter]->cd(iL+1-nPads*counter);
	  gPad->SetRightMargin(0.16);
	  p_occupancy[iL]->GetXaxis()->SetLabelSize(0.05);
	  p_occupancy[iL]->GetYaxis()->SetLabelSize(0.05);
	  p_occupancy[iL]->GetZaxis()->SetLabelSize(0.05);
	  p_occupancy[iL]->GetXaxis()->SetTitleSize(0.05);
	  p_occupancy[iL]->GetYaxis()->SetTitleSize(0.05);
	  p_occupancy[iL]->GetZaxis()->SetTitleSize(0.05);
	  p_occupancy[iL]->GetZaxis()->SetRangeUser(0,Emax[iO]);
	  p_occupancy[iL]->GetZaxis()->SetTitle("Occupancy(%)");
	  p_occupancy[iL]->GetZaxis()->SetTitleOffset(-0.5);
	  p_occupancy[iL]->Draw("colz");
	  gPad->Update();
	  TPaletteAxis *palette = (TPaletteAxis*)p_occupancy[iL]->GetListOfFunctions()->FindObject("palette");
	  palette->SetX1NDC(0.84);
	  palette->SetX2NDC(0.89);
	  palette->SetY1NDC(0.10);
	  palette->SetY2NDC(0.89);
	  TLatex lat;
	  lat.SetTextSize(0.07);
	  char buf[100];
	  sprintf(buf,"Layer %d",iL);
	  lat.DrawLatex(-80,110,buf);
	  if (iL%nPads==nPads-1){
	    myc[nPlots*iS+counter]->Update();
	    saveName.str("");
	    saveName << plotDir << "/occupancy_" << occThreshold[iO] << "_" << counter;
	    myc[nPlots*iS+counter]->Print((saveName.str()+".png").c_str());
	    myc[nPlots*iS+counter]->Print((saveName.str()+".pdf").c_str());
	    counter++;
	  }
	  
	}//loop on layers

      }//loop on occupancies

      myc[iS]->Clear();
      myc[iS]->cd();

      TGraphErrors *gr[nOcc];
      TLegend *leg = new TLegend(0.7,0.75,1.,1.);
      leg->SetFillColor(10);

      for (unsigned iO(0); iO<nOcc; ++iO){
	gr[iO] = new TGraphErrors(nLayers,lay,maxOcc[iO],layerr,maxOccErr[iO]);
	//gr[iO]->SetMaximum(Emax[0]);
	gr[iO]->SetMinimum(0);
	gr[iO]->SetMaximum(100);
	gr[iO]->GetXaxis()->SetTitle("Layer");
	gr[iO]->GetYaxis()->SetTitle("max occupancy (%)");
	gr[iO]->SetTitle(etaStr[iS]);
	gr[iO]->SetMarkerColor(4-iO);
	gr[iO]->SetMarkerStyle(23-iO);
	gr[iO]->SetLineColor(4-iO);
	gr[iO]->Draw(iO==0 ? "AP" : "P");
	char buf[100];
	sprintf(buf,"Threshold = %d MIPs",occThreshold[iO]);
	leg->AddEntry(gr[iO],buf,"P");
      }
      leg->Draw("same");
      myc[iS]->Update();
      saveName.str("");
      saveName << plotDir << "/maxOccupancy";
      myc[iS]->Print((saveName.str()+".png").c_str());
      myc[iS]->Print((saveName.str()+".pdf").c_str());



    }//loop on scenarios

  }//loop on versions
  
  return 0;


}//main
