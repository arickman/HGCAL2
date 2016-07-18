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

int plotMaxEhit() {

  std::string pSuffix = "300";

  const unsigned nS = 1;
  std::string scenario[nS] = {
   "e-/"
  };

  const unsigned nV = 1;
  TString version[nV] = {"12"};//,"0"};
 
  const unsigned nHists = 3;
  std::string label[nHists] = {"5x5mm2","10x10mm2","15x15mm2"};

  const unsigned nLayers = 30;
  unsigned genEn[]={5,10,15,20,25,30,40,50,60,100,150,200,300,400,500,1000,2000};
  const unsigned nGenEn=sizeof(genEn)/sizeof(unsigned);
  double genEnErr[nGenEn];
  double genEnD[nGenEn];
  TCanvas *mycL = new TCanvas("mycL","mycL",1500,1000);
  mycL->Divide(nHists,nGenEn);
  TCanvas *myc = new TCanvas("myc","myc",1);

  double mean[nHists][nGenEn];
  double meanerr[nHists][nGenEn];
  double rms[nHists][nGenEn];
  double rmserr[nHists][nGenEn];
  double mean3rms[nHists][nGenEn];
  double mean3rmserr[nHists][nGenEn];

  char buf[500];
  TLatex lat;

  for (unsigned iV(0); iV<nV;++iV){//loop on versions
    for (unsigned iS(0); iS<nS;++iS){//loop on scenarios

      TString plotDir = "../PLOTS/gitV00-02-02/version"+version[iV]+"/"+scenario[iS]+"/";
 
      TH1F *maxEhit[nHists];

      for (unsigned iE(0); iE<nGenEn; ++iE){
	std::cout << "- Processing energy : " << genEn[iE] 
		  << std::endl;
	genEnErr[iE] = 0;
	genEnD[iE] = genEn[iE];

	TFile *inputFile = 0;
	std::ostringstream linputStr;
	linputStr << plotDir << "validation_" << pSuffix << "um_e" << genEn[iE] << ".root";
	inputFile = TFile::Open(linputStr.str().c_str());
	if (!inputFile) {
	  std::cout << " -- Error, input file " << linputStr.str() << " cannot be opened. Exiting..." << std::endl;
	  return 1;
	}
	else std::cout << " -- File " << inputFile->GetName() << " sucessfully opened." << std::endl;
	  
	//maxEhit[0] = (TH1F*)gDirectory->Get("p_maxEhit_2d5");
	maxEhit[0] = (TH1F*)gDirectory->Get("p_maxEhit_5");
	maxEhit[1] = (TH1F*)gDirectory->Get("p_maxEhit_10");
	maxEhit[2] = (TH1F*)gDirectory->Get("p_maxEhit_15");

	for (unsigned iH(0); iH<nHists;++iH){
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
	  sprintf(buf,"%s #mu m Si",pSuffix.c_str());
	  lat.DrawLatexNDC(0.12,0.75,buf);
	  sprintf(buf,"Geant4 Standalone, e- %d GeV, pad size = %s",genEn[iE],label[iH].c_str());
	  lat.DrawLatexNDC(0.1,0.92,buf);

	  mean[iH][iE] = maxEhit[iH]->GetMean();
	  meanerr[iH][iE] = maxEhit[iH]->GetMeanError();
	  rms[iH][iE] = maxEhit[iH]->GetRMS();
	  rmserr[iH][iE] = maxEhit[iH]->GetRMSError();
	  mean3rms[iH][iE] = maxEhit[iH]->GetMean()+3*maxEhit[iH]->GetRMS();
	  mean3rmserr[iH][iE] = maxEhit[iH]->GetMeanError()+3*maxEhit[iH]->GetRMSError();

	  std::ostringstream lsave;
	  lsave << plotDir << pSuffix << "um/maxEhit_" << label[iH] << "_" << genEn[iE] << "GeV.pdf";
	  myc->Print(lsave.str().c_str());

	  mycL->cd(nHists*iE+iH+1);
	  maxEhit[iH]->Draw();
	  lat.SetTextSize(0.07);
	  sprintf(buf,"Maximum = %3.0f MIPs",maxEhit[iH]->GetXaxis()->GetBinCenter(maxEhit[iH]->FindLastBinAbove(0)));
	  lat.DrawLatexNDC(0.2,0.8,buf);
	  sprintf(buf,"%s #mu m Si",pSuffix.c_str());
	  lat.DrawLatexNDC(0.2,0.7,buf);
	  sprintf(buf,"Geant4 Standalone, e- %d GeV, pad size = %s",genEn[iE],label[iH].c_str());
	  lat.DrawLatexNDC(0.1,0.92,buf);

	}//loop on hists

      }//loop on energies
      std::ostringstream lsave;
      lsave << plotDir << pSuffix<< "um/maxEhit_all.pdf";
      mycL->Print(lsave.str().c_str());


      TGraphErrors *gr[nHists][2];
      TLegend *leg[2];
      for (unsigned i(0); i<2;++i){
	leg[i] = new TLegend(0.6,0.2,0.88,0.48);
	leg[i]->SetFillColor(10);
	myc->cd();
	gPad->SetLogx(1);
	gPad->SetLogy(1);
	gPad->SetGridx(1);
	gPad->SetGridy(1);
	TF1 *fitfunc[nHists];
	for (unsigned iH(0);iH<nHists;++iH){
	  if (i==0) gr[iH][i] = new TGraphErrors(nGenEn,genEnD,mean[iH],genEnErr,rms[iH]);
	  else gr[iH][i] = new TGraphErrors(nGenEn,genEnD,mean3rms[iH],genEnErr,mean3rmserr[iH]); 
	  gr[iH][i]->GetXaxis()->SetTitle("Beam energy (GeV)");
	  if (i==0) gr[iH][i]->GetYaxis()->SetTitle("<maxEhit> (MIPs)");
	  else gr[iH][i]->GetYaxis()->SetTitle("<maxE>+3#times RMS (MIPs)");
	  gr[iH][i]->SetTitle("");
	  gr[iH][i]->SetMarkerColor(iH+1);
	  gr[iH][i]->SetLineColor(iH+1);
	  gr[iH][i]->SetMarkerStyle(20+iH);
	  gr[iH][i]->SetMaximum(20000);
	  gr[iH][i]->Draw(iH>0?"P":"AP");
	  leg[i]->AddEntry(gr[iH][i],label[iH].c_str(),"P");
	  
	  if (i==0){
	    gr[iH][i]->Fit("pol2","","same");
	    fitfunc[iH] = gr[iH][i]->GetFunction("pol2");
	  }
	  else {
	    gr[iH][i]->Fit("pol1","","same",20,2000);
	    fitfunc[iH] = gr[iH][i]->GetFunction("pol1");
	  }
	  fitfunc[iH]->SetLineColor(iH+1);
	  
	  if (i==0) sprintf(buf,"<maxE> = (%3.0f #pm %3.0f)  + (%3.1f #pm %3.1f)#times E + (%1.0e#pm%1.0e) #times E^{2}",fitfunc[iH]->GetParameter(0),fitfunc[iH]->GetParError(0),fitfunc[iH]->GetParameter(1),fitfunc[iH]->GetParError(1),fitfunc[iH]->GetParameter(2),fitfunc[iH]->GetParError(2));
	  else sprintf(buf,"<maxE> = (%3.1f #pm %3.1f)  + (%3.2f #pm %3.2f)#times E",fitfunc[iH]->GetParameter(0),fitfunc[iH]->GetParError(0),fitfunc[iH]->GetParameter(1),fitfunc[iH]->GetParError(1));
	  lat.SetTextColor(iH+1);
	  lat.SetTextSize(0.04);
	  //lat.DrawLatexNDC(0.05,0.95-0.05*iH,buf);

	}
	
	sprintf(buf,"HGCAL Geant4 Standalone Simulation, e-");
	lat.SetTextColor(1);
	lat.SetTextSize(0.04);
	lat.DrawLatexNDC(0.01,0.01,buf);

	leg[i]->Draw("same");
	lsave.str("");
	lsave << plotDir << pSuffix << "um/summary";
	if (i==1) lsave << "_meanplus3rms";
	lsave << ".pdf";
	myc->Update();
	myc->Print(lsave.str().c_str());

	//return 1;
      }

    }//loop on scenarios
  }//loop on versions
  

  return 0;
}
