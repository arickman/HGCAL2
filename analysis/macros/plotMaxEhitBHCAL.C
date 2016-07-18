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

int plotMaxEhitBHCAL() {

  std::string pSuffix = "";

  const unsigned nS = 1;
  std::string scenario[nS] = {
   "u_pt"
  };

  const unsigned nV = 1;
  TString version[nV] = {"33"};//,"0"};
 
  const unsigned nHists = 2;
  std::string etastr[nHists] = {"eta18","eta18"};//,"eta25","eta29"};
  std::string label[nHists] = {"u quark","single pi^{-}"};
  const double leta[nHists] = {1.8,1.8};//,2.5,2.9};

  unsigned genEn[]={1,2,3};
  const unsigned nGenEn=sizeof(genEn)/sizeof(unsigned);

  unsigned genEt[nGenEn][nHists] = 
    {{320,320},//,170,110},
     {650,650},//,320,220},
     {1000,1000}//,500,330}
    };
  double genEnErr[nGenEn];
  double genEnD[nGenEn];
  TCanvas *mycL = new TCanvas("mycL","mycL",1500,1000);
  mycL->Divide(nHists,nGenEn);
  TCanvas *myc = new TCanvas("myc","myc",1);
  
  const unsigned nP = 3;
  TCanvas *mycmax[nP];
  mycmax[0] = new TCanvas("mycmax1","mycmax1",1);
  mycmax[1] = new TCanvas("mycmax3","mycmax3",1);
  mycmax[2] = new TCanvas("mycmax95","mycmax95",1);

  double mean[nHists][nGenEn];
  double meanerr[nHists][nGenEn];
  double rms[nHists][nGenEn];
  double rmserr[nHists][nGenEn];
  double mean3rms[nHists][nGenEn];
  double mean3rmserr[nHists][nGenEn];
  double int95[nHists][nGenEn];
  double int95err[nHists][nGenEn];

  char buf[500];
  TLatex lat;

  for (unsigned iV(0); iV<nV;++iV){//loop on versions
    for (unsigned iS(0); iS<nS;++iS){//loop on scenarios

      std::ostringstream plotDir;
      plotDir << "../PLOTS/gitV06a-03-05/version" << version[iV] << "/";
 
      TH1F *maxEhit[nHists];

      for (unsigned iE(0); iE<nGenEn; ++iE){
	std::cout << "- Processing energy : " << genEn[iE] 
		  << std::endl;
	genEnErr[iE] = 0;
	genEnD[iE] = genEn[iE];

	for (unsigned iH(0); iH<nHists;++iH){
	  TFile *inputFile = 0;
	  std::ostringstream linputStr;
	  linputStr << plotDir.str();
	  if (iH==0) linputStr << "/" << scenario[iS] << genEt[iE][0];
	  else linputStr << "/pi-";
	  linputStr << "/validation_" << etastr[iH];
	  if (iH==1) linputStr << "_et" << genEt[iE][iH];
	  linputStr << ".root";
	  inputFile = TFile::Open(linputStr.str().c_str());
	  if (!inputFile) {
	    std::cout << " -- Error, input file " << linputStr.str() << " cannot be opened. Exiting..." << std::endl;
	    return 1;
	  }
	  else std::cout << " -- File " << inputFile->GetName() << " sucessfully opened." << std::endl;
	  
	  //maxEhit[0] = (TH1F*)gDirectory->Get("p_maxEhit_2d5");
	  maxEhit[iH] = (TH1F*)gDirectory->Get("p_maxEhit_bh");

	  gStyle->SetOptStat("eMR");
	  myc->cd();
	  gStyle->SetStatX(0.89);
	  gStyle->SetStatY(0.89);
	  gStyle->SetStatH(0.5);
	  gStyle->SetStatW(0.3);
	  maxEhit[iH]->Rebin(genEn[iE]>500? 20*iE : (iE>0?2*iE:1));
	  maxEhit[iH]->GetXaxis()->SetRange(maxEhit[iH]->FindFirstBinAbove(0),maxEhit[iH]->FindLastBinAbove(0));
	  maxEhit[iH]->GetXaxis()->SetTitle("max E_{hit} (MIPS)");
	  maxEhit[iH]->Draw();

	  sprintf(buf,"Maximum = %3.0f MIPs",maxEhit[iH]->GetXaxis()->GetBinCenter(maxEhit[iH]->FindLastBinAbove(0)));
	  lat.SetTextSize(0.05);
	  lat.DrawLatexNDC(0.12,0.85,buf);
	  //sprintf(buf,"%s #mu m Si",pSuffix.c_str());
	  lat.DrawLatexNDC(0.12,0.75,label[iH].c_str());
	  sprintf(buf,"Geant4 Standalone, %3.3f TeV",genEt[iE][iH]*cosh(leta[iH]));
	  lat.DrawLatexNDC(0.1,0.92,buf);

	  mean[iH][iE] = maxEhit[iH]->GetMean();
	  meanerr[iH][iE] = maxEhit[iH]->GetMeanError();
	  rms[iH][iE] = maxEhit[iH]->GetRMS();
	  rmserr[iH][iE] = maxEhit[iH]->GetRMSError();
	  mean3rms[iH][iE] = maxEhit[iH]->GetMean()+3*maxEhit[iH]->GetRMS();
	  mean3rmserr[iH][iE] = maxEhit[iH]->GetMeanError()+3*maxEhit[iH]->GetRMSError();
	  double total = maxEhit[iH]->Integral();
	  for (int iB(0); iB<maxEhit[iH]->GetNbinsX()+2;++iB){
	    double integral = maxEhit[iH]->Integral(0,iB);
	    if (integral/total>0.95) {
	      int95[iH][iE] = maxEhit[iH]->GetXaxis()->GetBinCenter(iB);
	      int95err[iH][iE] = maxEhit[iH]->GetXaxis()->GetBinWidth(iB);
	      break;
	    }
	  }

	  std::cout << label[iH] << " E=" << genEt[iE][iH]*cosh(leta[iH]) 
		    << " & $" << mean[iH][iE] << " #pm " << rms[iH][iE]
		    << "$ & $" << mean3rms[iH][iE] << " #pm " << mean3rmserr[iH][iE]
		    << "$ & $" << int95[iH][iE] << " #pm " << int95err[iH][iE]
		    << "$ \\\\"
		    << std::endl;

	  std::ostringstream lsave;
	  lsave << plotDir.str() << pSuffix << "/BHmaxEhit_" << etastr[iH] << "_" << genEn[iE] << "TeV.pdf";
	  myc->Print(lsave.str().c_str());

	  mycL->cd(nHists*iE+iH+1);
	  maxEhit[iH]->Draw();
	  lat.SetTextSize(0.07);
	  sprintf(buf,"Maximum = %3.0f MIPs",maxEhit[iH]->GetXaxis()->GetBinCenter(maxEhit[iH]->FindLastBinAbove(0)));
	  lat.DrawLatexNDC(0.2,0.8,buf);
	  //sprintf(buf,"%s #mu m Si",pSuffix.c_str());
	  //lat.DrawLatexNDC(0.2,0.7,buf);
	  //sprintf(buf,"Geant4 Standalone, e- %d GeV, pad size = %s",genEn[iE],label[iH].c_str());
	  lat.DrawLatexNDC(0.12,0.75,label[iH].c_str());
	  sprintf(buf,"Geant4 Standalone, %3.3f TeV",genEt[iE][iH]*cosh(leta[iH]));
	  lat.DrawLatexNDC(0.1,0.92,buf);

	}//loop on hists

      }//loop on energies
      std::ostringstream lsave;
      lsave << plotDir.str() << pSuffix<< "/BHmaxEhit_all.pdf";
      mycL->Print(lsave.str().c_str());

      TGraphErrors *gr[nHists][nP];
      TLegend *leg[nP];
      leg[0] = new TLegend(0.6,0.5,0.88,0.78);
      leg[1] = new TLegend(0.15,0.5,0.33,0.78);
      leg[2] = new TLegend(0.15,0.5,0.33,0.78);
      for (unsigned i(0); i<nP;++i){
	leg[i]->SetFillColor(10);
	mycmax[i]->cd();
	//gPad->SetLogx(1);
	//gPad->SetLogy(1);
	gPad->SetGridx(1);
	gPad->SetGridy(1);
	TF1 *fitfunc[nHists];
	for (unsigned iH(0);iH<nHists;++iH){
	  if (i==0) gr[iH][i] = new TGraphErrors(nGenEn,genEnD,mean[iH],genEnErr,rms[iH]);
	  else if (i==1) gr[iH][i] = new TGraphErrors(nGenEn,genEnD,mean3rms[iH],genEnErr,mean3rmserr[iH]); 
	  else gr[iH][i] = new TGraphErrors(nGenEn,genEnD,int95[iH],genEnErr,int95err[iH]); 
	  gr[iH][i]->GetXaxis()->SetTitle("Beam energy (TeV)");
	  if (i==0) gr[iH][i]->GetYaxis()->SetTitle("<maxEhit> (MIPs)");
	  else if (i==1) gr[iH][i]->GetYaxis()->SetTitle("<maxE>+3#times RMS (MIPs)");
	  else gr[iH][i]->GetYaxis()->SetTitle("95% content  (MIPs)");
	  gr[iH][i]->SetTitle("");
	  gr[iH][i]->SetMarkerColor(iH+1);
	  gr[iH][i]->SetLineColor(iH+1);
	  gr[iH][i]->SetMarkerStyle(20+iH);
	  gr[iH][i]->SetMinimum(0);
	  gr[iH][i]->SetMaximum(1000);
	  gr[iH][i]->Draw(iH>0?"P":"AP");
	  leg[i]->AddEntry(gr[iH][i],label[iH].c_str(),"P");
	  
	  gr[iH][i]->Fit("pol1","","same",0.5,3.5);
	  fitfunc[iH] = gr[iH][i]->GetFunction("pol1");
	  if (!fitfunc[iH]) continue;

	  fitfunc[iH]->SetLineColor(iH+1);
	  
	  //if (i==0) sprintf(buf,"<maxE> = (%3.0f #pm %3.0f)  + (%3.1f #pm %3.1f)#times E + (%1.0e#pm%1.0e) #times E^{2}",fitfunc[iH]->GetParameter(0),fitfunc[iH]->GetParError(0),fitfunc[iH]->GetParameter(1),fitfunc[iH]->GetParError(1),fitfunc[iH]->GetParameter(2),fitfunc[iH]->GetParError(2));
	  //else 
	  sprintf(buf,"<maxE> = (%3.1f #pm %3.1f)  + (%3.2f #pm %3.2f)#times E",fitfunc[iH]->GetParameter(0),fitfunc[iH]->GetParError(0),fitfunc[iH]->GetParameter(1),fitfunc[iH]->GetParError(1));
	  lat.SetTextColor(iH+1);
	  lat.SetTextSize(0.04);
	  lat.DrawLatexNDC(0.05,0.95-0.05*iH,buf);

	}
	
	sprintf(buf,"HGCAL Geant4 Standalone Simulation, pi-");
	lat.SetTextColor(1);
	lat.SetTextSize(0.04);
	lat.DrawLatexNDC(0.01,0.01,buf);

	leg[i]->Draw("same");
	lsave.str("");
	lsave << plotDir.str() << pSuffix << "/BHsummary";
	if (i==1) lsave << "_meanplus3rms";
	else if (i==2) lsave << "_95percent";
	lsave << ".pdf";
	mycmax[i]->Update();
	mycmax[i]->Print(lsave.str().c_str());

	//return 1;
      }

    }//loop on scenarios
  }//loop on versions
  

  return 0;
}
