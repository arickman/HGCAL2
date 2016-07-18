#include<string>
#include<iostream>
#include<fstream>
#include<sstream>

#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
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


int plotXY(){//main  

  bool doXYplots = false;
  const double Emip = 0.0548;//in MeV

  //const unsigned nS = 7;
  //TString scenario[nS] = {"0","1","2","3","4","5","6"};
  const unsigned nS = 1;
  std::string scenario[nS] = {
    "pi-/twiceSampling/GeVCal/EarlyDecay/"
    //"pi+/PU/eta20/",
    //"pi+/PU/eta25/",
    //"pi+/PU/eta30/"
    //"e-/noWeights/"
    //"scenario_0/GammaSmeared/eta20/",
    //"scenario_0/GammaSmeared/eta25/",
    //"scenario_0/GammaSmeared/eta30/",
    //"scenario_0/GammaSmeared/eta35/",
    //"scenario_0/GammaPU/eta20/",
    //"scenario_0/GammaPU/eta25/",
    //"scenario_0/GammaPU/eta30/",
    //"scenario_0/GammaPU/eta35/",
    //"scenario_0/PU/eta20/",
    //"scenario_0/PU/eta25/",
    //"scenario_0/PU/eta30/",
    //"scenario_0/PU/eta35/"
};
  
  const unsigned nV = 1;
  TString version[nV] = {"23"};
  bool isPU = false;
  const unsigned nLayers = 54;
  const unsigned nEcalLayers = 38;

  TCanvas *mycECAL = new TCanvas("mycECAL","mycECAL",1500,1000);
  TCanvas *mycHCAL = new TCanvas("mycHCAL","mycHCAL",1500,1000);
  TCanvas *mycAll = new TCanvas("mycAll","mycAll",1500,1000);
  TCanvas *myc = new TCanvas("myc","myc",1);
  
  for (unsigned iV(0); iV<nV;++iV){//loop on versions
    for (unsigned iS(0); iS<nS;++iS){//loop on scenarios
    
      //if (scenario[iS].find("PU") != scenario[iS].npos) isPU = true;

      TString plotDir = "../PLOTS/version"+version[iV]+"/"+scenario[iS]+"/";
      //plotDir += "noWeights/";
      //TString plotDir = "../PLOTS/version_"+version[iV]+"/scenario_"+scenario[iS]+"/";
      //TString plotDir = "../PLOTS/version_"+version[iV]+"/";

      TFile *inputFile = TFile::Open(plotDir+"CalibHistos.root");
      if (!inputFile) {
	std::cout << " -- Error, input file " << plotDir << "/CalibHistos.root cannot be opened. Exiting..." << std::endl;
	return 1;
      }
      

      //unsigned genEn[]={5,10,25,50,75,100,150,200,300,500};
      unsigned genEn[]={50};
      const unsigned nGenEn=sizeof(genEn)/sizeof(unsigned);
  
      TH2F *p_xy[nGenEn][nLayers];
      TH1F *p_Etot[nGenEn][nLayers];
      TH2F *p_EvsLayer[nGenEn];
      TProfile *prof_EvsLayer[nGenEn];
      double Emax[nGenEn];
      std::ostringstream saveName;

      Double_t X0[nEcalLayers];
      double X0tot = 0;
      for (unsigned iL(0); iL<nEcalLayers; ++iL){
	X0tot += iL<10 ? 0.5 : (iL<20 ? 0.8 : 1.2);
	X0[iL] = X0tot;
      }
      TH2F *p_EvsX0 = 0;
      TProfile *prof_EvsX0 = 0;

      for (unsigned iE(0); iE<nGenEn; ++iE){

	std::cout << "- Processing energy : " << genEn[iE] << std::endl;
	Emax[iE] = 0;

	double EmaxLayer = 0;
    
	bool stopProcessing = false;

	std::ostringstream lName;

	for (unsigned iL(0); iL<nLayers; ++iL){
	  lName.str("");
	  lName << "p_xy_" << genEn[iE] << "_" << iL;
	  p_xy[iE][iL] = (TH2F*)gDirectory->Get(lName.str().c_str());
	  if (!p_xy[iE][iL]) {
	    std::cout << " -- ERROR, pointer for histogram is null for layer: " << iL << ". Going to next one..." << std::endl;
	    stopProcessing=true;
	    break;
	  }
	  double Etot = p_xy[iE][iL]->GetMaximum();
	  if (Etot > Emax[iE]) Emax[iE] = Etot;

	  lName.str("");
	  lName << "p_Etot_" << genEn[iE] << "_" << iL;
	  p_Etot[iE][iL] = (TH1F*)gDirectory->Get(lName.str().c_str());
	  if (!p_Etot[iE][iL]) {
	    std::cout << " -- ERROR, pointer for Etot histogram is null for layer: " << iL << ". Exiting..." << std::endl;
	    return 1;
	  }
	  if (p_Etot[iE][iL]->GetMean() > EmaxLayer) EmaxLayer = p_Etot[iE][iL]->GetMean()+5*p_Etot[iE][iL]->GetRMS();
	}
	if (stopProcessing) continue;

	std::cout << " -- max energy " << Emax[iE] << std::endl;

	lName.str("");
	lName << "p_EvsLayer_" << genEn[iE];
	p_EvsLayer[iE] = (TH2F*)gDirectory->Get(lName.str().c_str());

	prof_EvsLayer[iE] = p_EvsLayer[iE]->ProfileX();

	mycECAL->Clear();
	unsigned counter = 1;
	if (doXYplots) mycECAL->Divide(6,5);
	else mycECAL->Divide(3,3);

	mycHCAL->Clear();
	if (doXYplots) mycHCAL->Divide(6,5);
	else mycHCAL->Divide(3,3);

	gStyle->SetOptStat(0);
    
	for (unsigned iL(0); iL<nLayers; ++iL){//loop on layers
	  std::cout << " -- Processing layer " << iL << std::endl;
	  if (doXYplots) {
	    if (iL<nEcalLayers) mycECAL->cd(iL+1);
	    else mycHCAL->cd(iL+1);
	  }
	  else {
	    if (iL%3==2 && counter < 10 && iL<nEcalLayers) {
	      mycECAL->cd(counter);
	      counter++;
	    }
	    else if (iL%2==1 && counter<19 && iL>=nEcalLayers){
	      std::cout << " -- Hcal layer " << iL << " counter = " << counter << std::endl;
	      mycHCAL->cd(counter-9);
	      counter++;
	    }
	    else continue;
	  }
	  gPad->SetLogz(1);
	  p_xy[iE][iL]->SetMaximum(Emax[iE]);
	  p_xy[iE][iL]->Scale(1./Emip);
	  p_xy[iE][iL]->GetXaxis()->SetLabelSize(0.06);
	  p_xy[iE][iL]->GetYaxis()->SetLabelSize(0.06);
	  p_xy[iE][iL]->GetXaxis()->SetTitleSize(0.05);
	  p_xy[iE][iL]->GetYaxis()->SetTitleSize(0.05);
	  p_xy[iE][iL]->GetZaxis()->SetTitleOffset(-0.5);
	  p_xy[iE][iL]->GetZaxis()->SetTitle("E(MIP)");

	  p_xy[iE][iL]->Draw("colz");
	  TLatex lat;
	  lat.SetTextSize(0.07);
	  char buf[500];
	  sprintf(buf,"Layer %d",iL);
	  lat.DrawLatex(-40,105,buf);
	  if (doXYplots) {
	    myc->cd();
	    p_xy[iE][iL]->Draw("colz");
	    myc->Update();
	    saveName.str("");
	    saveName << plotDir << "/xySimHits/xySimHits_layer" << iL << "_" << genEn[iE] << "GeV";
	    myc->Print((saveName.str()+".png").c_str());
	    myc->Print((saveName.str()+".pdf").c_str());
	    myc->cd();
	    myc->SetLogz(1);
	    p_xy[iE][iL]->Draw("colz");
	    myc->Update();
	    saveName.str("");
	    saveName << plotDir << "/xySimHits/xySimHits_layer" << iL << "_" << genEn[iE] << "GeV_log";
	    myc->Print((saveName.str()+".png").c_str());
	    myc->Print((saveName.str()+".pdf").c_str());
	  }
	}//loop on layers
	saveName.str("");
	if (doXYplots) saveName << plotDir << "/xySimHits_ECAL_" << genEn[iE] << "GeV";
	else saveName << plotDir << "/xySimHits_ECAL_" << genEn[iE] << "GeV_subset";
	mycECAL->Update();
	mycECAL->Print((saveName.str()+".png").c_str());
	mycECAL->Print((saveName.str()+".pdf").c_str());
	saveName.str("");
	if (doXYplots) saveName << plotDir << "/xySimHits_HCAL_" << genEn[iE] << "GeV";
	else saveName << plotDir << "/xySimHits_HCAL_" << genEn[iE] << "GeV_subset";
	mycHCAL->Update();
	mycHCAL->Print((saveName.str()+".png").c_str());
	mycHCAL->Print((saveName.str()+".pdf").c_str());

	//return 1;

	myc->cd();
	//p_EvsLayer[iE]->GetYaxis()->SetRangeUser(0,1000);
	p_EvsLayer[iE]->Draw("colz");
	prof_EvsLayer[iE]->SetMarkerStyle(23);
	prof_EvsLayer[iE]->SetMarkerColor(1);
	prof_EvsLayer[iE]->Draw("PEsame");
	myc->Update();
	saveName.str("");
	saveName << plotDir << "/ElayervsLayer_" << genEn[iE] << "GeV";
	myc->Print((saveName.str()+".png").c_str());
	myc->Print((saveName.str()+".pdf").c_str());

	if (iE==5){
	  p_EvsX0 = new TH2F("EvsX0",";Absorber X_{0}; E (MeV)",50,0,25,
			     //30,X0,
			     100,0,200);
	  //p_EvsLayer[iE]->GetYaxis()->GetBinLowEdge(1),
	  //p_EvsLayer[iE]->GetYaxis()->GetBinLowEdge(p_EvsLayer[iE]->GetNbinsY()+1)
	  //);
	  p_EvsX0->Sumw2();
	  for (int iX(0); iX<nEcalLayers;++iX){
	    for (int iY(1); iY<p_EvsLayer[iE]->GetNbinsY()+1;++iY){
	      if (iY==200) std::cout << iX << " " << iY << " " << X0[iX] << " " << p_EvsLayer[iE]->GetYaxis()->GetBinCenter(iY) << " " << p_EvsLayer[iE]->GetBinContent(iX+1,iY) << std::endl;
	      p_EvsX0->Fill(X0[iX],p_EvsLayer[iE]->GetYaxis()->GetBinCenter(iY),p_EvsLayer[iE]->GetBinContent(iX+1,iY)<0?0:p_EvsLayer[iE]->GetBinContent(iX+1,iY));
	      //p_EvsX0->SetBinContent(iX+1,iY,p_EvsLayer[iE]->GetBinContent(iX+1,iY));
	      //p_EvsX0->SetBinError(iX+1,iY,p_EvsLayer[iE]->GetBinError(iX+1,iY));
	    }
	  }
	  //p_EvsX0->RebinY(10);
	  p_EvsX0->Draw("colz");
	  prof_EvsX0 = p_EvsX0->ProfileX();
	  prof_EvsX0->SetMarkerStyle(23);
	  prof_EvsX0->SetMarkerColor(1);
	  prof_EvsX0->Draw("PEsame");
	  myc->Update();
	  saveName.str("");
	  saveName << plotDir << "/ElayervsX0_" << genEn[iE] << "GeV";
	  myc->Print((saveName.str()+".png").c_str());
	  myc->Print((saveName.str()+".pdf").c_str());

	  //return 1;
	}


      }//loop on energies


      mycAll->Clear();
      mycAll->Divide(5,2);
      
      for (unsigned iE(0); iE<nGenEn; ++iE){
      	mycAll->cd(iE+1);
      	p_EvsLayer[iE]->Draw("colz");
      	prof_EvsLayer[iE]->SetMarkerStyle(23);
      	prof_EvsLayer[iE]->SetMarkerColor(1);
      	prof_EvsLayer[iE]->Draw("PEsame");
      	TLatex lat;
      	char buf[500];
      	sprintf(buf,"Egen = %d GeV",genEn[iE]);
      	lat.DrawLatex(1,p_EvsLayer[iE]->GetMaximum()*1.1,buf);
      }//loop on energies

      mycAll->Update();
      saveName.str("");
      saveName << plotDir << "/ElayervsLayer";
      mycAll->Print((saveName.str()+".png").c_str());
      mycAll->Print((saveName.str()+".pdf").c_str());


    }//loop on scenarios

  }//loop on versions

  return 0;


}//main
