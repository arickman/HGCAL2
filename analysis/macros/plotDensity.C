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

double getBinWidth(const double etaMin,const double etaMax, std::string etaStr, const double xLaymm){
  double eta0 = 0;
  if (etaStr.find("eta20")!=etaStr.npos) eta0 = 2.0;
  else if (etaStr.find("eta25")!=etaStr.npos) eta0 = 2.5;
  else if (etaStr.find("eta30")!=etaStr.npos) eta0 = 3.0;
  else if (etaStr.find("eta35")!=etaStr.npos) eta0 = 3.5;
  else {
     std::cout << " ERROR: new eta point ! Please implement value to convert y position... Exiting..." << std::endl;
    exit(1);
  }

  double xDist = 310+xLaymm/10.;//distance in cm from beam spot to layer.
  double theta0 = 2*atan(exp(-eta0));
  double y0 = tan(theta0)*xDist;//cm

  double theta1 = 2*atan(exp(-etaMin));
  double ypos1 = xDist*tan(theta1)-y0;

  double theta2 = 2*atan(exp(-etaMax));
  double ypos2 = xDist*tan(theta2)-y0;

  double dy = ypos1-ypos2;//cm

  // std::cout << " - etaBin " << etaStr 
  // 	    << " etaMin = " << etaMin 
  // 	    << " etaMax = " << etaMax
  // 	    << " dy = " << dy
  // 	    << std::endl;

  return dy;

}

int plotDensity(){//main

  //const unsigned nS = 7;
  //TString scenario[nS] = {"0","1","2","3","4","5","6"};
  const unsigned nS = 4;
  std::string scenario[nS] = {
    "scenario_0/PedroPU/eta20/",
    "scenario_0/PedroPU/eta25/",
    "scenario_0/PedroPU/eta30/",
    "scenario_0/PedroPU/eta35/"
  };

  const unsigned nV = 1;
  TString version[nV] = {"20"};
  const double signalRegionInX=2;//cm

  double eta[4] = {2.0,2.5,3.,3.5};
  
  const double Emip = 0.056;//in MeV

  std::ostringstream saveName;
  bool isPU = false;
  
  const unsigned nLayers = 30;
      
  unsigned genEn[]={5,10,25,50,75,100,150,200,300,500};
  //unsigned genEn[]={10};
  const unsigned nGenEn=sizeof(genEn)/sizeof(unsigned);

  TH1F *p_Edensity[nV][nS][nGenEn][nLayers];

  double xLay[nLayers];
  for (unsigned iL(0); iL<nLayers; ++iL){
    xLay[iL] = xLayer(iL);
  }

  //canvas so they are created only once
  TCanvas *mycE = new TCanvas("mycE","mycE",1500,1000);
  TCanvas *mycAll = new TCanvas("mycAll","mycAll",1500,1000);
  mycAll->Divide(3,3);

  const unsigned nCanvas = 4;  
  TCanvas *myc[nCanvas];
  for (unsigned iC(0);iC<nCanvas;++iC){
    std::ostringstream lName;
    lName << "myc" << iC;
    myc[iC] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1);
  }

  for (unsigned iV(0); iV<nV;++iV){//loop on versions
    for (unsigned iS(0); iS<nS;++iS){//loop on scenarios
    
      if (scenario[iS].find("PU") != scenario[iS].npos) isPU = true;
  
      TString plotDir = "../PLOTS/version_"+version[iV]+"/"+scenario[iS]+"/";
      //plotDir += "noWeights/";
      //TString plotDir = "../PLOTS/version_"+version[iV]+"/scenario_"+scenario[iS]+"/";
      //TString plotDir = "../PLOTS/version_"+version[iV]+"/";

      TFile *inputFile = TFile::Open(plotDir+"CalibHistos.root");
      if (!inputFile) {
	std::cout << " -- Error, input file " << plotDir << "/CalibHistos.root cannot be opened. Exiting..." << std::endl;
	return 1;
      }
      
      for (unsigned iE(0); iE<nGenEn; ++iE){

	if (scenario[iS].find("Pedro")!=scenario[iS].npos) genEn[iE]=iE;

	std::cout << "- Processing energy : " << genEn[iE] 
		  << std::endl;
	
	std::ostringstream lName;
	for (unsigned iL(0); iL<nLayers; ++iL){
	  lName.str("");
	  lName << "p_Edensity_" << genEn[iE] << "_" << iL;
	  p_Edensity[iV][iS][iE][iL] = (TH1F*)(gDirectory->Get(lName.str().c_str())->Clone());
	  if (!p_Edensity[iV][iS][iE][iL]) {
	    std::cout << " -- ERROR, pointer for histogram Edensity is null for layer: " << iL << ". Exiting..." << std::endl;
	    return 1;
	  }
	  //std::cout << " ---- Edensity = " << p_Edensity[iV][iS][iE][iL]->GetMean() << std::endl;
	}//loop on layers


      }//loop on energies


    }//loop on scenarios
    
  }//loop on versions


  for (unsigned iV(0); iV<nV;++iV){//loop on versions
 
    TH1F *p_allEtaDensity[nLayers];
    TH2F *p_allEtaDensityLayervsEta = new TH2F("p_allEtaDensityLayervsEta",";#eta;Layer;E density (MIP/cm^{2})",61,1.88,4.32,30,0,30);
    TH2F *p_allEtaDensityEtavsLayer = new TH2F("p_allEtaDensityEtavsLayer",";Layer;#eta;E density (MIP/cm^{2})",30,0,30,61,1.88,4.32);
    p_allEtaDensityLayervsEta->Sumw2();
    p_allEtaDensityEtavsLayer->Sumw2();

    TH1F *p_etavsLayer[4];
    p_etavsLayer[0] = new TH1F("p_eta20vsLayer",";layer;E density (MIP/cm^{2})",30,0,30);
    p_etavsLayer[1] = new TH1F("p_eta25vsLayer",";layer;E density (MIP/cm^{2})",30,0,30);
    p_etavsLayer[2] = new TH1F("p_eta30vsLayer",";layer;E density (MIP/cm^{2})",30,0,30);
    p_etavsLayer[3] = new TH1F("p_eta35vsLayer",";layer;E density (MIP/cm^{2})",30,0,30);

    unsigned counter = 1;
    TLine * line20;
    TLine * line25;
    TLine * line30;
    TLine * line35;
    TLine * linem20;
    TLine * linem25;
    TLine * linem30;
    TLine * linem35;

    unsigned eSignal = 9;

    for (unsigned iL(0); iL<nLayers; ++iL){//loop on layers
      //plot Edensity
      std::ostringstream lName;
      lName.str("");
      lName << "p_allEtaDensity_" << iL;
      p_allEtaDensity[iL] = new TH1F(lName.str().c_str(),";#eta",130,1.8,4.4);
      p_allEtaDensity[iL]->Sumw2();
      
      double minEta[nS];
      double maxEta[nS];
      
      //fill histos
      //merge all energy files
      for (unsigned iE(isPU?0:eSignal); iE<(isPU?nGenEn:eSignal+1); ++iE){//loop on energies
	//unsigned iE = 5; 
	for (unsigned iS(0); iS<nS;++iS){//loop on scenarios
	  //just plot for pure PU
	  //if (scenario[iS].find("/PU") == scenario[iS].npos) continue;
	  
	  minEta[iS] = 5;
	  maxEta[iS] = 0;
	  
	  //get min and max for each scenario
	  for (int iB(1); iB<p_Edensity[iV][iS][iE][iL]->GetNbinsX()+1;++iB){//loop on bins
	    double etaMin = p_Edensity[iV][iS][iE][iL]->GetBinLowEdge(iB);
	    double etaMax = p_Edensity[iV][iS][iE][iL]->GetBinLowEdge(iB+1);
	    if (p_Edensity[iV][iS][iE][iL]->GetBinContent(iB) > 0){
	      if (etaMin < minEta[iS]) minEta[iS] = etaMin;
	      if (etaMax > maxEta[iS]) maxEta[iS] = etaMax;
	    }
	  }
	}//loop on scenarios
	for (unsigned iS(0); iS<nS;++iS){//loop on scenarios
	  //just plot for pure PU
	  //if (scenario[iS].find("/PU") == scenario[iS].npos) continue;

	  int binEta = p_Edensity[iV][iS][iE][iL]->GetXaxis()->FindBin(eta[iS]);

	  for (int iB(1); iB<p_Edensity[iV][iS][iE][iL]->GetNbinsX()+1;++iB){//loop on bins
	    bool overlap = false;
	    double etaMin = p_Edensity[iV][iS][iE][iL]->GetBinLowEdge(iB);
	    double etaMax = p_Edensity[iV][iS][iE][iL]->GetBinLowEdge(iB+1);
	    double etaCenter = p_Edensity[iV][iS][iE][iL]->GetBinCenter(iB);
	    
	    if ( (iS < nS-1 && etaCenter > minEta[iS+1]) ||
		 (iS > 0 && etaCenter < maxEta[iS-1])) overlap = true;
	    
	    double norm = (isPU? 1/1000.: 1/100.)*1/signalRegionInX*1/getBinWidth(etaMin,etaMax,scenario[iS],xLay[iL]); 
	    double density = p_Edensity[iV][iS][iE][iL]->GetBinContent(iB)*norm;
	    if (overlap) density = density/2.;
	    //std::cout << density << std::endl;
	    p_allEtaDensity[iL]->Fill(etaCenter,density/Emip);
	    p_allEtaDensityLayervsEta->Fill(etaCenter,iL,density/Emip);
	    p_allEtaDensityEtavsLayer->Fill(iL,etaCenter,density/Emip);

	    if (iB==binEta) p_etavsLayer[iS]->Fill(iL,density/Emip);

	  }//loop on bins
	}//loop on scenarios
      }//loop on energies
      
      p_allEtaDensity[iL]->Rebin(2);
      p_allEtaDensity[iL]->GetYaxis()->SetTitle("E density (MIP/cm^{2})");
      p_allEtaDensity[iL]->GetXaxis()->SetLabelSize(0.05);
      p_allEtaDensity[iL]->GetYaxis()->SetLabelSize(0.05);
      p_allEtaDensity[iL]->GetXaxis()->SetTitleSize(0.05);
      p_allEtaDensity[iL]->GetYaxis()->SetTitleSize(0.05);
      double maxY = isPU?250:100;
      p_allEtaDensity[iL]->GetYaxis()->SetRangeUser(0,maxY);
      if (iL%3==2 && counter < 10) {
	mycAll->cd(counter);
	counter++;
	
	p_allEtaDensity[iL]->Draw("PE");
	TLatex lat;
	char buf[500];
	lat.SetTextSize(0.07);
	sprintf(buf,"Layer %d",iL);
	lat.DrawLatex(2,maxY*1.01,buf);

	line20 = new TLine(minEta[0],0,minEta[0],maxY);
	line20->SetLineColor(2);
	line20->Draw();
	lat.SetTextColor(2);
	lat.DrawLatex(minEta[0],maxY*0.3,"Eta 2.0");
	line25 = new TLine(minEta[1],0,minEta[1],maxY);
	line25->SetLineColor(3);
	line25->Draw();
	lat.SetTextColor(3);
	lat.DrawLatex(minEta[1],maxY*0.5,"Eta 2.5");
	line30 = new TLine(minEta[2],0,minEta[2],maxY);
	line30->SetLineColor(4);
	line30->Draw();
	lat.SetTextColor(4);
	lat.DrawLatex(minEta[2],maxY*0.7,"Eta 3.0");
	line35 = new TLine(minEta[3],0,minEta[3],maxY);
	line35->SetLineColor(6);
	line35->Draw();
	lat.SetTextColor(6);
	lat.DrawLatex(minEta[3],maxY*0.9,"Eta 3.5");

	linem20 = new TLine(maxEta[0],0,maxEta[0],maxY);
	linem20->SetLineColor(2);
	linem20->Draw();
	linem25 = new TLine(maxEta[1],0,maxEta[1],maxY);
	linem25->SetLineColor(3);
	linem25->Draw();
	linem30 = new TLine(maxEta[2],0,maxEta[2],maxY);
	linem30->SetLineColor(4);
	linem30->Draw();
	linem35 = new TLine(maxEta[3],0,maxEta[3],maxY);
	linem35->SetLineColor(6);
	linem35->Draw();

      }

    }//loop on layers

    TString plotDir = "../PLOTS/version_"+version[iV]+"/"+scenario[0]+"/";

    saveName.str("");
    saveName << plotDir << "/../Edensity_subset";
    if (!isPU) saveName << "_" << genEn[eSignal] << "GeV";
    mycAll->Update();
    mycAll->Print((saveName.str()+".png").c_str());
    mycAll->Print((saveName.str()+".pdf").c_str());
    //for (unsigned iL(0); iL<nLayers; ++iL){
    //p_allEtaDensity[iL]->Delete();
    //}
    mycE->Clear();
    mycE->cd();
    gPad->SetLogz(1);
    p_allEtaDensityLayervsEta->GetXaxis()->SetLabelSize(0.05);
    p_allEtaDensityLayervsEta->GetYaxis()->SetLabelSize(0.05);
    p_allEtaDensityLayervsEta->GetXaxis()->SetTitleSize(0.05);
    p_allEtaDensityLayervsEta->GetYaxis()->SetTitleSize(0.05);
    p_allEtaDensityLayervsEta->GetZaxis()->SetTitleOffset(-0.5);
    p_allEtaDensityLayervsEta->GetZaxis()->SetRangeUser(isPU?0.1:0.01,isPU?250:300);
    p_allEtaDensityLayervsEta->Draw("colz");
    TLatex lat;
    char buf[500];
    sprintf(buf,"E_{MIP} = %2.1f keV",Emip*1000);
    lat.DrawLatex(2,31,buf);

    saveName.str("");
    saveName << plotDir << "/../EdensityLayervsEta";
    if (!isPU) saveName << "_" << genEn[eSignal] << "GeV";
    mycE->Update();
    mycE->Print((saveName.str()+".png").c_str());
    mycE->Print((saveName.str()+".pdf").c_str());

    mycE->Clear();
    mycE->cd();
    gPad->SetLogz(1);
    p_allEtaDensityEtavsLayer->GetXaxis()->SetLabelSize(0.05);
    p_allEtaDensityEtavsLayer->GetYaxis()->SetLabelSize(0.05);
    p_allEtaDensityEtavsLayer->GetZaxis()->SetLabelSize(0.05);
    p_allEtaDensityEtavsLayer->GetXaxis()->SetTitleSize(0.05);
    p_allEtaDensityEtavsLayer->GetYaxis()->SetTitleSize(0.05);
    p_allEtaDensityEtavsLayer->GetZaxis()->SetTitleSize(0.05);
    p_allEtaDensityEtavsLayer->GetZaxis()->SetTitleOffset(-0.5);
    p_allEtaDensityEtavsLayer->GetZaxis()->SetRangeUser(isPU?0.1:0.01,isPU?250:300);
    p_allEtaDensityEtavsLayer->Draw("lego2z");
    sprintf(buf,"E_{MIP} = %2.1f keV",Emip*1000);
    lat.DrawLatex(3,4.35,buf);

    saveName.str("");
    saveName << plotDir << "/../EdensityEtavsLayer";
    if (!isPU) saveName << "_" << genEn[eSignal] << "GeV";
    mycE->Update();
    mycE->Print((saveName.str()+".png").c_str());
    mycE->Print((saveName.str()+".pdf").c_str());
    
    for (unsigned iS(0); iS<nS;++iS){//loop on scenarios
      myc[iS]->cd();
      p_etavsLayer[iS]->Draw();
      std::cout << "-------- Eta = " << eta[iS] << " --------" << std::endl;
      for (int iB(1); iB<p_etavsLayer[iS]->GetNbinsX()+1;++iB){
	std::cout << "-- Bin " << iB << " " << p_etavsLayer[iS]->GetBinContent(iB) << std::endl;
      }

      saveName.str("");
      saveName << "../PLOTS/version_" << version[iV] << "/DensityPerLayer_eta" << eta[iS] ;
      myc[iS]->Update();
      myc[iS]->Print((saveName.str()+".png").c_str());
      myc[iS]->Print((saveName.str()+".pdf").c_str());
    }

  }//loop on versions

  return 0;
  
  
}//main
