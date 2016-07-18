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

int plotLateralSize(){//main

  //const unsigned nS = 7;
  //TString scenario[nS] = {"0","1","2","3","4","5","6"};
  const unsigned nS = 1;
  std::string scenario[nS] = {
    //"pi-/twiceSampling/GeVCal/EarlyDecay/",
    "pi-/concept/GeVCal/"
  };
  
  const unsigned nV = 1;
  TString version[nV] = {"21"};//,"0"};
  
  //const unsigned nLayers = 54;//9;//33;

  TString pDetector = "ECAL";
  
  const unsigned nLat = 20;
  const double latSize[nLat] = {10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200};//+/- in mm

  std::ostringstream saveName;

  unsigned genEn[]={5,10,25,40,50,60,80,100,150,200,300,400,500};//,1000};//,1000,2000};
  //unsigned genEn[]={40,50,60,80,100,200,300,400,500,1000,2000};
  //unsigned genEn[]={5,10,20,25,50,75,100,125,150,175,200,300,500};
  //unsigned genEn[]={5,20,50,100,150,200};
  //unsigned genEn[]={50};
  const unsigned nGenEn=sizeof(genEn)/sizeof(unsigned);
  TCanvas *mycF = new TCanvas("mycF","mycF",1500,1000);
  TCanvas *mycI = new TCanvas("mycI","mycI",1500,750);
  TCanvas *mycGr = new TCanvas("mycGr","mycGr",1500,1000);
  mycGr->Divide(5,3);
  TCanvas *mycE = new TCanvas("mycE","mycE",1000,1000);

  for (unsigned iV(0); iV<nV;++iV){//loop on versions
    for (unsigned iS(0); iS<nS;++iS){//loop on scenarios
      
      TString plotDir = "../PLOTS/version"+version[iV]+"/"+scenario[iS]+"/";
  
      TFile *inputFile = TFile::Open(plotDir+"CalibHistos.root");
      if (!inputFile) {
	std::cout << " -- Error, input file " << plotDir << "/CalibHistos.root cannot be opened. Trying next..." << std::endl;
	continue;
      }
      
      TH1F *p_EfracLateral[nGenEn][nLat];

      TH2F *p_meanFrac = new TH2F("p_meanFrac",";radius (mm);gen E (GeV);<E_{x#times x}>/E_{tot}",nLat,5,latSize[nLat-1]+5,nGenEn,0,nGenEn);
      TH2F *p_rmsFrac = new TH2F("p_rmsFrac",";radius (mm);gen E (GeV);#sigma(E_{x#times x}/E_{tot})",nLat,5,latSize[nLat-1]+5,nGenEn,0,nGenEn);
      double meanFrac[nGenEn][nLat];
      double rmsFrac[nGenEn][nLat];
      double latErr[nLat];

      double f90[nGenEn];
      double f90err[nGenEn];
      double E[nGenEn];
      double Eerr[nGenEn];

      TGraphErrors *gr[nGenEn];

      for (unsigned iE(0); iE<nGenEn; ++iE){
	
	std::cout << "- Processing energy : " << genEn[iE] 
		  << std::endl;
	std::ostringstream lName;
	
	E[iE] = genEn[iE];
	Eerr[iE] = 0;

	TString eStr;
	eStr += genEn[iE];
	p_meanFrac->GetYaxis()->SetBinLabel(iE+1,eStr);
	p_rmsFrac->GetYaxis()->SetBinLabel(iE+1,eStr);

	mycF->Clear();
	mycF->Divide(5,4);
	
	for (unsigned iR(0); iR<nLat; ++iR){
	  latErr[iR] = 0;
	  lName.str("");
	  lName << "p_EfracLateral_" << genEn[iE] << "_" << pDetector << "_" << iR;
	  p_EfracLateral[iE][iR] = (TH1F*)gDirectory->Get(lName.str().c_str());
	  if (!p_EfracLateral[iE][iR]) {
	    std::cout << " -- ERROR, pointer for histogram Efrac is null for energy " << genEn[iE] << " radius: " << iR << ". Exiting..." << std::endl;
	    return 1;
	  }

	  double overflows = p_EfracLateral[iE][iR]->GetBinContent(p_EfracLateral[iE][iR]->GetNbinsX()+1);
	  if (overflows>0) {
	    std::cout << " -- ERROR, overflows are not 0 : " << overflows << " entries: " << p_EfracLateral[iE][iR]->GetEntries() << std::endl;
	  } 

	  p_meanFrac->SetBinContent(iR+1,iE+1,p_EfracLateral[iE][iR]->GetMean());
	  p_rmsFrac->SetBinContent(iR+1,iE+1,p_EfracLateral[iE][iR]->GetRMS());

	  meanFrac[iE][iR] = p_EfracLateral[iE][iR]->GetMean();
	  rmsFrac[iE][iR] = p_EfracLateral[iE][iR]->GetRMS();

	  p_EfracLateral[iE][iR]->Sumw2();
	  //p_EfracLateral[iE][iR]->Rebin(4);
	  //gStyle->SetOptStat("eMRuoi");
	  gStyle->SetOptStat(0);
	  mycF->cd(iR+1);
	  p_EfracLateral[iE][iR]->Draw("PE");
	  double maxy = p_EfracLateral[iE][iR]->GetMaximum();

	  TLatex lat;
	  char buf[500];
	  sprintf(buf,"pi^{-}, E = %d GeV",genEn[iE]);
	  lat.DrawLatex(0.1,maxy*1.2,buf);
	  sprintf(buf,"Lateral size: %3.0f #times %3.0f mm^{2}",2*latSize[iR],2*latSize[iR]);
	  lat.DrawLatex(0.1,maxy*1,buf);
  

	}//loop on radius

	saveName.str("");
	saveName << plotDir << "/EtotLateralFraction_" << genEn[iE] << "GeV";
	mycF->Update();
	mycF->Print((saveName.str()+".png").c_str());
	mycF->Print((saveName.str()+".pdf").c_str());

	gStyle->SetOptStat(0);

	mycGr->cd(iE+1);
	gr[iE] = new TGraphErrors(nLat,latSize,meanFrac[iE],latErr,rmsFrac[iE]);
	gr[iE]->SetMarkerStyle(20);
	gr[iE]->SetMarkerColor(1);
	gr[iE]->SetLineColor(1);
	gr[iE]->GetXaxis()->SetTitle("Radius x (mm)");
	gr[iE]->GetYaxis()->SetTitle("<E_{x #times x}/E_{tot}>");
	char buf[500];
	sprintf(buf,"pi^{-}, E = %d GeV",genEn[iE]);
	gr[iE]->SetTitle(buf);
	gr[iE]->SetMinimum(0);
	gr[iE]->SetMaximum(1.1);
	gr[iE]->Draw("APL");
      
	for (unsigned iR(0); iR<nLat; ++iR){
	  if (meanFrac[iE][iR] > 0.9){
	    f90[iE] = latSize[iR];
	    f90err[iE] = 5;
	    break;
	  }
	}

      }//loop on energies

      mycGr->Update();
      saveName.str("");
      saveName << plotDir << "/meanLatFraction_Graph";
      mycGr->Update();
      mycGr->Print((saveName.str()+".png").c_str());
      mycGr->Print((saveName.str()+".pdf").c_str());
  
      mycE->cd();
      TGraphErrors *gr90 = new TGraphErrors(nGenEn,E,f90,Eerr,f90err);
      gr90->SetMarkerStyle(20);
      gr90->SetMarkerColor(1);
      gr90->SetLineColor(1);
      gr90->GetXaxis()->SetTitle("Energy (GeV)");
      gr90->GetYaxis()->SetTitle("90% radius (mm)");
      gr90->GetYaxis()->SetTitleOffset(1.3);
      gr90->SetTitle("");
      gr90->SetMinimum(0);
      //gr90->SetMaximum(150);
      gr90->Draw("APL");

      saveName.str("");
      saveName << plotDir << "/f90vsE";
      mycE->Update();
      mycE->Print((saveName.str()+".png").c_str());
      mycE->Print((saveName.str()+".pdf").c_str());
  
      //draw energy fractions
      gStyle->SetOptStat(0);
      mycI->Clear();
      mycI->Divide(2,1);
      mycI->cd(1);
      p_meanFrac->Draw("colz");
      mycI->cd(2);
      p_rmsFrac->Draw("colz");
      
      saveName.str("");
      saveName << plotDir << "/meanLatFraction";
      mycI->Update();
      mycI->Print((saveName.str()+".png").c_str());
      mycI->Print((saveName.str()+".pdf").c_str());
  

    }//loop on scenarios
    
  }//loop on versions

  return 0;

}//main
