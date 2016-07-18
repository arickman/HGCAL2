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
	
bool checkHistPointer(TH1 *pointer, std::string name){
  if (!pointer) {
    std::cout << " -- ERROR, pointer for histogram " << name << " is null. Exiting..." << std::endl;
    return false;
  }
  return true;
}

int plotGlobalCor(){//main  

  //const unsigned nS = 7;
  //TString scenario[nS] = {"0","1","2","3","4","5","6"};
  const unsigned nS = 1;
  std::string scenario[nS] = {
    "pi-/"
  };
  
  const unsigned nV = 1;
  TString version[nV] = {"23"};

  TString pSuffix = "";

  const unsigned nLayers = 54;//64 //54
  //const unsigned nHcalLayers = 38;

  unsigned genEn[]={//10,30,80};
    10,15,18,20,25,
    30,35,40,45,50,
    60,80};//,100,200,300,400,500
  //,1000,2000};
  //unsigned genEn[]={10,15,20,25,30,40,50,60,80,100,150,200,300,400,500};//,1000,2000};
  //unsigned genEn[]={50,100,150,200};
  //unsigned genEn[]={10,80};
  unsigned nGenEn=sizeof(genEn)/sizeof(unsigned);

  //const unsigned nLimits = 1;//5;
  //const double pElim[nLimits] = {5};
  const unsigned nLimits = 9;//5;
  const double pElim[nLimits] = {5,7.5,10,15,20,25,30,35,40};//,45,50,60,70,80};

  const unsigned limRef = 2;

  std::ostringstream lName;

  const unsigned nC = 7;
  const unsigned nCtot = nC+1;
  TCanvas *myc[nCtot];
  for (unsigned iC(0); iC<nCtot;++iC){
    lName.str("");
    lName << "Canvas_" << iC;
    myc[iC] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1000,700);
  }
  myc[0]->Divide(3,2);
  myc[1]->Divide(3,3);
  myc[2]->Divide(5,3);
  myc[3]->Divide(5,3);
  myc[4]->Divide(2,1);


  std::ostringstream saveName;
  
  for (unsigned iV(0); iV<nV;++iV){//loop on versions
    for (unsigned iS(0); iS<nS;++iS){//loop on scenarios
    
      TString plotDir = "../PLOTS/gitV00-02-04/version"+version[iV]+"/"+scenario[iS]+"/";
      
      TH1F *p_hitSpectrum_lowTail[nGenEn];
      TH1F *p_hitSpectrum_highTail[nGenEn];
      TH1F *p_hitSpectrum_ratio[nGenEn];
      TH1F *p_meanHitSpectrum_lowTail[nGenEn];
      TH1F *p_meanHitSpectrum_highTail[nGenEn];
      TH1F *p_Cglobal[nGenEn][nLimits];
      TH2F *p_EvsCglobal[nGenEn][nLimits];
      TH1F *p_Etotal[nGenEn];
      TH1F *p_Eshower[nGenEn][nLimits];

      double genEnErr[nGenEn];
      double grgenEn[nGenEn];

      unsigned lColor[20] = {1,2,3,4,5,6,7,8,9,209,
			     210,211,212,213,214,215,216,217,218,219};

      TLegend *leg = new TLegend(0.4,0.75,0.89,0.89);
      leg->SetFillColor(10);
      TLegend *legE = new TLegend(0.7,0.5,0.89,0.89);
      legE->SetFillColor(10);
      TLegend *legR = new TLegend(0.905,0.1,1.0,0.9);
      legR->SetFillColor(10);

      //draw calibration curves
      TGraphErrors *calib = new TGraphErrors();
      calib->SetName("calib");
      calib->SetMarkerStyle(21);
      calib->SetTitle("");

      char buf[500];

      TGraphErrors *reso[nLimits+1];
      for (unsigned iLim(0); iLim<nLimits+1;++iLim){
	reso[iLim] = new TGraphErrors();
	lName.str("");
	lName << "reso_"<<iLim;
	reso[iLim]->SetName(lName.str().c_str());
	reso[iLim]->SetMarkerStyle(20+iLim);
	reso[iLim]->SetMarkerColor(lColor[iLim]);
	reso[iLim]->SetLineColor(lColor[iLim]);
	reso[iLim]->SetTitle("");
	if (iLim<nLimits) sprintf(buf,"%1.1f",pElim[iLim]);
	else sprintf(buf,"Ref");
	legR->AddEntry(reso[iLim],buf,"P");
      }
      reso[nLimits]->SetMarkerStyle(21);
      reso[nLimits]->SetMarkerColor(2);
      reso[nLimits]->SetLineColor(2);

      for (unsigned iE(0); iE<nGenEn; ++iE){

	std::cout << " -- Processing energy " << genEn[iE] << std::endl;
	genEnErr[iE] = 0;
	grgenEn[iE] = genEn[iE];

	lName.str("");
	lName << plotDir << "globalCompensation_e" << genEn[iE] << pSuffix << ".root";
	TFile *inputFile = TFile::Open(lName.str().c_str());
	if (!inputFile) {
	  std::cout << " -- Error, input file " << lName.str() << " cannot be opened. Going to next..." << std::endl;
	  continue;
	  //return 1;
	}

	lName.str("");
	lName << "p_hitSpectrum_lowTail";
	p_hitSpectrum_lowTail[iE]  = (TH1F*)gDirectory->Get(lName.str().c_str());
	if (!checkHistPointer(p_hitSpectrum_lowTail[iE],lName.str())) return 1;
	p_hitSpectrum_lowTail[iE]->Scale(1./p_hitSpectrum_lowTail[iE]->GetEntries());

	lName.str("");
	lName << "p_hitSpectrum_highTail";
	p_hitSpectrum_highTail[iE]  = (TH1F*)gDirectory->Get(lName.str().c_str());
	if (!checkHistPointer(p_hitSpectrum_highTail[iE],lName.str())) return 1;
	p_hitSpectrum_highTail[iE]->Scale(1./p_hitSpectrum_highTail[iE]->GetEntries());

	lName.str("");
	lName << "p_hitSpectrum_ratio";
	p_hitSpectrum_ratio[iE] = (TH1F*)p_hitSpectrum_highTail[iE]->Clone(lName.str().c_str());
	p_hitSpectrum_ratio[iE]->Divide(p_hitSpectrum_highTail[iE],p_hitSpectrum_lowTail[iE]);

	lName.str("");
	lName << "p_meanHitSpectrum_lowTail";
	p_meanHitSpectrum_lowTail[iE]  = (TH1F*)gDirectory->Get(lName.str().c_str());
	if (!checkHistPointer(p_meanHitSpectrum_lowTail[iE],lName.str())) return 1;
	p_meanHitSpectrum_lowTail[iE]->Scale(1./p_meanHitSpectrum_lowTail[iE]->GetEntries());

	lName.str("");
	lName << "p_meanHitSpectrum_highTail";
	p_meanHitSpectrum_highTail[iE]  = (TH1F*)gDirectory->Get(lName.str().c_str());
	if (!checkHistPointer(p_meanHitSpectrum_highTail[iE],lName.str())) return 1;
	p_meanHitSpectrum_highTail[iE]->Scale(1./p_meanHitSpectrum_highTail[iE]->GetEntries());

	for (unsigned iLim(0); iLim<nLimits;++iLim){
	  lName.str("");
	  lName << "p_Cglobal_" << iLim;
	  p_Cglobal[iE][iLim] =  (TH1F*)gDirectory->Get(lName.str().c_str());
	  if (!checkHistPointer(p_Cglobal[iE][iLim],lName.str())) return 1;
	  
	  lName.str("");
	  lName << "p_Eshower_" << iLim;
	  p_Eshower[iE][iLim] = (TH1F*)gDirectory->Get(lName.str().c_str());
	  if (!checkHistPointer(p_Eshower[iE][iLim],lName.str())) return 1;
	  
	  lName.str("");
	  lName << "p_EvsCglobal_" << iLim;
	  p_EvsCglobal[iE][iLim] = (TH2F*)gDirectory->Get(lName.str().c_str());
	  if (!checkHistPointer(p_EvsCglobal[iE][iLim],lName.str())) return 1;
	}
	
	lName.str("");
	lName << "p_ErecoTotal";
	p_Etotal[iE] = (TH1F*)gDirectory->Get(lName.str().c_str());
	if (!checkHistPointer(p_Etotal[iE],lName.str())) return 1;

	gStyle->SetOptStat(0);

	sprintf(buf,"#pi^{-} %d GeV",genEn[iE]);

	//myc[0]->cd(6);
	p_hitSpectrum_lowTail[iE]->SetLineColor(4);
	p_hitSpectrum_highTail[iE]->SetLineColor(2);
	//p_hitSpectrum_lowTail[iE]->SetLineWidth(iE%3+1);
	p_hitSpectrum_lowTail[iE]->GetYaxis()->SetTitle("bin Entries / Entries");
	p_hitSpectrum_lowTail[iE]->SetTitle(buf);
	p_hitSpectrum_lowTail[iE]->Rebin(4);
	p_hitSpectrum_highTail[iE]->Rebin(4);
	p_hitSpectrum_lowTail[iE]->GetXaxis()->SetRangeUser(0.5,20);
	p_hitSpectrum_highTail[iE]->GetXaxis()->SetRangeUser(0.5,20);
	//p_hitSpectrum_lowTail[iE]->GetYaxis()->SetRangeUser(0,0.025);

	p_hitSpectrum_ratio[iE]->GetYaxis()->SetTitle("Ratio highTail/LowTail");
	p_hitSpectrum_ratio[iE]->GetXaxis()->SetRangeUser(0.5,20);
	p_hitSpectrum_ratio[iE]->GetYaxis()->SetRangeUser(0.,2.);
	if (genEn[iE]==10) {
	  myc[0]->cd(1);
	  p_hitSpectrum_lowTail[iE]->Draw("");
	  p_hitSpectrum_highTail[iE]->Draw("same");
	  leg->AddEntry(p_hitSpectrum_lowTail[iE],"Eevent < Emean-sigma","L");
	  leg->AddEntry(p_hitSpectrum_highTail[iE],"Eevent > Emean+sigma","L");
	  leg->Draw("same");
	  myc[0]->cd(4);
	  gPad->SetGridx(1);
	  gPad->SetGridy(1);
	  p_hitSpectrum_ratio[iE]->Draw("");
	}
	else if (genEn[iE]==30){
	  myc[0]->cd(2);
	  p_hitSpectrum_lowTail[iE]->Draw("");
	  p_hitSpectrum_highTail[iE]->Draw("same");
	  myc[0]->cd(5);
	  gPad->SetGridx(1);
	  gPad->SetGridy(1);
	  p_hitSpectrum_ratio[iE]->Draw("");
	}
	else if (genEn[iE]==80){
	  myc[0]->cd(3);
	  p_hitSpectrum_lowTail[iE]->Draw("");
	  p_hitSpectrum_highTail[iE]->Draw("same");
	  myc[0]->cd(6);
	  gPad->SetGridx(1);
	  gPad->SetGridy(1);
	  p_hitSpectrum_ratio[iE]->Draw("");
	}

	myc[2]->cd(iE+1);
	p_Eshower[iE][limRef]->SetMarkerStyle(20);
	p_Eshower[iE][limRef]->SetMarkerColor(1);
	p_Eshower[iE][limRef]->Rebin(1);
	double minx = std::min(p_Eshower[iE][limRef]->GetMean()-5*p_Eshower[iE][limRef]->GetRMS(),
			       p_Etotal[iE]->GetMean()-5*p_Etotal[iE]->GetRMS());
	double maxx = std::max(p_Eshower[iE][limRef]->GetMean()+5*p_Eshower[iE][limRef]->GetRMS(),
			       p_Etotal[iE]->GetMean()+5*p_Etotal[iE]->GetRMS());
	p_Eshower[iE][limRef]->GetXaxis()->SetRangeUser(minx,maxx);
	p_Eshower[iE][limRef]->SetMaximum(std::max(p_Eshower[iE][limRef]->GetMaximum(),p_Etotal[iE]->GetMaximum()));
	p_Eshower[iE][limRef]->Draw("PE");
	p_Etotal[iE]->SetMarkerStyle(21);
	p_Etotal[iE]->SetMarkerColor(2);
	p_Etotal[iE]->Rebin(1);
	//p_Etotal[iE]->Scale(p_Eshower[iE][limRef]->GetEntries()/p_Etotal[iE]->GetEntries());
	p_Etotal[iE]->Draw("PEsame");
	std::cout << "shower " << p_Eshower[iE][limRef]->GetEntries() << " reco " << p_Etotal[iE]->GetEntries() << std::endl;

	myc[6]->cd();
	gStyle->SetOptStat("MRuoi");
	for (unsigned iLim(0); iLim<nLimits;++iLim){
	  
	  TH1F *hTmp = (TH1F*)p_Eshower[iE][iLim]->Clone();
	  hTmp->SetMarkerStyle(20);
	  hTmp->SetMarkerColor(1);
	  hTmp->GetXaxis()->SetRangeUser(hTmp->GetMean()-5*hTmp->GetRMS(),hTmp->GetMean()+5*hTmp->GetRMS());
	  if (iLim==limRef) {
	    myc[3]->cd(iE+1);
	    hTmp->Draw("PE");
	  }
	  else myc[6]->cd();
	  gStyle->SetOptFit(1111);
	  hTmp->Fit("gaus","LR+","",
		    hTmp->GetMean()-2*hTmp->GetRMS(),
		    hTmp->GetMean()+2*hTmp->GetRMS());
	  //sprintf(buf,"E = %d GeV",genEn[iE]);
	  hTmp->SetTitle(buf);
	  TF1 *fitResult = hTmp->GetFunction("gaus");
	  if (iLim==limRef) {
	    Int_t npc=calib->GetN();
	    calib->SetPoint(npc,fitResult->GetParameter(1),genEn[iE]);
	    calib->SetPointError(npc,fitResult->GetParError(1),0.01*genEn[iE]);
	  }
	  Int_t np = reso[iLim]->GetN();
	  double lreso=fitResult->GetParameter(2)/fitResult->GetParameter(1);
	  double lresoErr = lreso*sqrt(pow(fitResult->GetParError(1)/fitResult->GetParameter(1),2)+pow(fitResult->GetParError(2)/fitResult->GetParameter(2),2));
	  reso[iLim]->SetPoint(np,1/sqrt(genEn[iE]),lreso);
	  reso[iLim]->SetPointError(np,0,lresoErr);
	}

	if (genEn[iE] == 40){
	  for (unsigned iLim(0); iLim<nLimits;++iLim){
	    myc[1]->cd(iLim+1);
	    //p_EvsCglobal[iE][iLim]->RebinX(2);
	    //p_EvsCglobal[iE][iLim]->RebinY(2);
	    p_EvsCglobal[iE][iLim]->GetXaxis()->SetRangeUser(p_Cglobal[iE][iLim]->GetMean()-5*p_Cglobal[iE][iLim]->GetRMS(),p_Cglobal[iE][iLim]->GetMean()+5*p_Cglobal[iE][iLim]->GetRMS());
	    p_EvsCglobal[iE][iLim]->GetYaxis()->SetRangeUser(p_Etotal[iE]->GetMean()-5*p_Etotal[iE]->GetRMS(),p_Etotal[iE]->GetMean()+5*p_Etotal[iE]->GetRMS());
	    p_EvsCglobal[iE][iLim]->SetTitle(buf);
	    p_EvsCglobal[iE][iLim]->Draw("colz");
	  }
	}

	Int_t np = reso[nLimits]->GetN();
	reso[nLimits]->SetPoint(np,1/sqrt(genEn[iE]),p_Etotal[iE]->GetRMS()/p_Etotal[iE]->GetMean());
	reso[nLimits]->SetPointError(np,0,p_Etotal[iE]->GetRMSError()/p_Etotal[iE]->GetMean());

	myc[4]->cd(1);
	gStyle->SetOptStat(0);
	p_meanHitSpectrum_lowTail[iE]->SetLineColor(4);
	p_meanHitSpectrum_highTail[iE]->SetLineColor(2);
	//p_meanHitSpectrum_lowTail[iE]->SetLineWidth(iE%3+1);
	p_meanHitSpectrum_lowTail[iE]->Rebin(2);
	p_meanHitSpectrum_highTail[iE]->Rebin(2);
	p_meanHitSpectrum_lowTail[iE]->GetXaxis()->SetRangeUser(0,20);
	//p_meanHitSpectrum_lowTail[iE]->GetYaxis()->SetRangeUser(0,0.025);
	p_meanHitSpectrum_lowTail[iE]->GetYaxis()->SetTitle("Normalised entries");
	p_meanHitSpectrum_lowTail[iE]->SetTitle(buf);

	if (genEn[iE]==10) {
	  myc[4]->cd(1);
	  p_meanHitSpectrum_lowTail[iE]->Draw("");
	  p_meanHitSpectrum_highTail[iE]->Draw("same");
	  leg->Draw("same");
	}
	else if (genEn[iE]==80){
	  myc[4]->cd(2);
	  p_meanHitSpectrum_lowTail[iE]->Draw("");
	  p_meanHitSpectrum_highTail[iE]->Draw("same");
	  leg->Draw("same");
	}

	for (unsigned iLim(0); iLim<nLimits;++iLim){
	  if (nC+iLim<nCtot) myc[nC+iLim]->cd();
	  p_Cglobal[iE][iLim]->SetLineColor(lColor[iE]);
	  p_Cglobal[iE][iLim]->SetLineWidth(3);
	  //p_Cglobal[iE][iLim]->SetFillStyle();
	  //p_Cglobal[iE][iLim]->SetFillColor(lColor[iE]);
	  p_Cglobal[iE][iLim]->Scale(1./p_Cglobal[iE][iLim]->GetEntries());
	  p_Cglobal[iE][iLim]->Rebin(8);
	  p_Cglobal[iE][iLim]->GetXaxis()->SetRangeUser(0.8,1.7);
	  p_Cglobal[iE][iLim]->GetYaxis()->SetRangeUser(0,0.15);

	  sprintf(buf,"C_{global} (e_{lim} = %1.1f MIP)",pElim[iLim]);

	  p_Cglobal[iE][iLim]->GetXaxis()->SetTitle(buf);
	  p_Cglobal[iE][iLim]->GetYaxis()->SetTitle("Normalised entries");
	  if (iLim==limRef && (iE%3==0 || genEn[iE]==80)) {
	    p_Cglobal[iE][limRef]->Draw(iE==0?"":"same");
	    sprintf(buf,"%d GeV",genEn[iE]);
	    if (iLim==limRef) legE->AddEntry(p_Cglobal[iE][limRef],buf,"L");
	  }
	}
	//gPad->SetLogz(1);

      }//loop on energies
 
      saveName.str("");
      saveName << plotDir << "/HitSpectra";
      myc[0]->Update();
      myc[0]->Print((saveName.str()+".png").c_str());
      myc[0]->Print((saveName.str()+".pdf").c_str());

      saveName.str("");
      saveName << plotDir << "/EtotvsCglobal";
      myc[1]->Update();
      myc[1]->Print((saveName.str()+".png").c_str());
      myc[1]->Print((saveName.str()+".pdf").c_str());

      saveName.str("");
      saveName << plotDir << "/EuncorEshower";
      myc[2]->Update();
      myc[2]->Print((saveName.str()+".png").c_str());
      myc[2]->Print((saveName.str()+".pdf").c_str());

      saveName.str("");
      saveName << plotDir << "/Eshower";
      myc[3]->cd();
      gStyle->SetOptStat("MRuoi");
      myc[3]->Update();
      myc[3]->Print((saveName.str()+".png").c_str());
      myc[3]->Print((saveName.str()+".pdf").c_str());

      saveName.str("");
      saveName << plotDir << "/EmipMean";
      gStyle->SetOptStat(0);
      myc[4]->Update();
      myc[4]->Print((saveName.str()+".png").c_str());
      myc[4]->Print((saveName.str()+".pdf").c_str());

      myc[5]->cd();
      gStyle->SetOptStat(0);
      gStyle->SetOptFit(0);
      TGraphErrors * gr = calib;
      gr->GetXaxis()->SetLabelSize(0.06);
      gr->GetXaxis()->SetTitleSize(0.06);
      gr->GetYaxis()->SetLabelSize(0.06);
      gr->GetYaxis()->SetTitleSize(0.06);
      gr->GetXaxis()->SetTitleOffset(0.7);
      gr->GetYaxis()->SetTitleOffset(0.8);
      gr->SetTitle("");
      gr->GetYaxis()->SetTitle("Beam energy [GeV]");
      gr->GetXaxis()->SetTitle("Eshower [GeV]");
      gr->Draw("ap");
      TF1 *fitFunc=new TF1("calibFunc","[0]*x+[1]*x*x+[2]*x*x*x",gr->GetXaxis()->GetXmin(),gr->GetXaxis()->GetXmax());
      fitFunc->SetParameters(0.98,0.002,0);
      fitFunc->SetLineColor(6);
      gr->Fit(fitFunc,"RME");
      TLatex lat;
      lat.SetTextColor(6);
      sprintf(buf,"<E> #propto E(a + b #times E +c #times E^2)");
      lat.DrawLatex(genEn[0],gr->GetYaxis()->GetXmax()*0.9,buf);
      sprintf(buf,"a = %3.3f #pm %3.3f ",fitFunc->GetParameter(0),fitFunc->GetParError(0));
      lat.DrawLatex(genEn[0],gr->GetYaxis()->GetXmax()*(0.8),buf);
      sprintf(buf,"b = %3.3f #pm %3.3f",fitFunc->GetParameter(1),fitFunc->GetParError(1));
      lat.DrawLatex(genEn[0],gr->GetYaxis()->GetXmax()*(0.7),buf);
      sprintf(buf,"c = %3.1e #pm %3.1e",fitFunc->GetParameter(2),fitFunc->GetParError(2));
      lat.DrawLatex(genEn[0],gr->GetYaxis()->GetXmax()*(0.6),buf);
      sprintf(buf,"chi2/NDF = %3.3f/%d = %3.3f",fitFunc->GetChisquare(),fitFunc->GetNDF(),fitFunc->GetChisquare()/fitFunc->GetNDF());
      lat.DrawLatex(genEn[0],gr->GetYaxis()->GetXmax()*(0.5),buf);

      saveName.str("");
      saveName << plotDir << "/EbeamvsEshower";
      myc[5]->Update();
      myc[5]->Print((saveName.str()+".png").c_str());
      myc[5]->Print((saveName.str()+".pdf").c_str());

      myc[6]->cd();
      gr = reso[0];
      gr->GetXaxis()->SetLabelSize(0.06);
      gr->GetXaxis()->SetTitleSize(0.06);
      gr->GetYaxis()->SetLabelSize(0.06);
      gr->GetYaxis()->SetTitleSize(0.06);
      gr->GetXaxis()->SetTitleOffset(0.7);
      gr->GetYaxis()->SetTitleOffset(0.8);
      gr->SetTitle("");
      gr->GetXaxis()->SetTitle("1/#sqrt{Beam energy} [1/#sqrt{GeV}]");
      gr->GetYaxis()->SetTitle("#sigma/E");
      gr->SetMinimum(0);
      gr->SetMaximum(0.25);
      gr->Draw("ap");
      TF1 *fitFunc2 =new TF1("resoFunc","sqrt([0]*[0]*x*x+[1]*[1])",gr->GetXaxis()->GetXmin(),gr->GetXaxis()->GetXmax());

      for (unsigned iLim(0); iLim<nLimits;++iLim){
	if (nC+iLim<nCtot){
	  saveName.str("");
	  saveName << plotDir << "/Cglobal_Elim" << pElim[iLim];
	  myc[nC+iLim]->cd();
	  legE->Draw("same");
	  myc[nC+iLim]->Update();
	  myc[nC+iLim]->Print((saveName.str()+".png").c_str());
	  myc[nC+iLim]->Print((saveName.str()+".pdf").c_str());
	}

	myc[6]->cd();
	reso[iLim]->Draw("p");

      }
      reso[nLimits]->Draw("pL");

      saveName.str("");
      saveName << plotDir << "/resolutionVsElim";
      myc[6]->cd();

      legR->Draw("same");
      myc[6]->Update();
      myc[6]->Print((saveName.str()+".png").c_str());
      myc[6]->Print((saveName.str()+".pdf").c_str());



    }//loop on scenarios
    
  }//loop on versions
  
  return 0;


}//main
