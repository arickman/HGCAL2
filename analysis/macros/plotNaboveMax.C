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

int plotNaboveMax() {

  std::string basename = "p_nAboveMax_";
  //std::string basename = "p_fracEmissed_";

  std::string pSuffix = "300";
  std::string pPrefix = "Thresh3000/";

  bool integrate = true;

  const unsigned nS = 1;
  std::string scenario[nS] = {
   "e-/"
  };

  const unsigned nV = 1;
  TString version[nV] = {"12"};//,"0"};
 
  const unsigned nHists = 3;
  std::string label[nHists] = {"5x5mm2","10x10mm2","15x15mm2"};

  const unsigned nLayers = 30;
  //unsigned genEn[]={1000};
  //unsigned genEn[]={5,10,15,20,25,30,40,50,60,100,150,200,300,400,500,1000,2000};
  unsigned genEn[]={400,500,1000,2000};
  const unsigned nGenEn=sizeof(genEn)/sizeof(unsigned);
  double genEnErr[nGenEn];
  double genEnD[nGenEn];
  TCanvas *mycL = new TCanvas("mycL","mycL",1500,1000);
  mycL->Divide(2,2);
  TCanvas *myc = new TCanvas("myc","myc",1);

  char buf[500];
  TLatex lat;
  

  for (unsigned iV(0); iV<nV;++iV){//loop on versions
    for (unsigned iS(0); iS<nS;++iS){//loop on scenarios

      TString plotDir = "../PLOTS/gitV00-02-02/version"+version[iV]+"/"+scenario[iS]+"/"+pPrefix;
 
      TH1F *maxNhit[nHists];

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

	//maxNhit[0] = (TH1F*)gDirectory->Get(basename+"2d5");
	//if (iE==0) 
	maxNhit[0] = (TH1F*)gDirectory->Get((basename+"5").c_str());
	//else maxNhit[0]->Add((TH1F*)gDirectory->Get((basename+"5").c_str()));
	//if (iE==0) 
	maxNhit[1] = (TH1F*)gDirectory->Get((basename+"10").c_str());
	//else maxNhit[1]->Add((TH1F*)gDirectory->Get((basename+"10").c_str()));
	//if (iE==0) 
	maxNhit[2] = (TH1F*)gDirectory->Get((basename+"15").c_str());
	//else maxNhit[2]->Add((TH1F*)gDirectory->Get((basename+"15").c_str()));
	//	}//loop on energies

	if (integrate){
	  for (unsigned iH(0); iH<nHists;++iH){
	    double inttot = maxNhit[iH]->Integral(0,maxNhit[iH]->GetNbinsX()+1);
	    for (int bin(0); bin<maxNhit[0]->GetNbinsX()+2;++bin){
	      double err = 0;
	      double integral = 
		maxNhit[iH]->IntegralAndError(bin,maxNhit[iH]->GetNbinsX()+1,err);
	      /*for (int bin2(bin); bin2<maxNhit[0]->GetNbinsX()+2;++bin2){
		integral += maxNhit[iH]->GetBinContent(bin2);
		err += pow(maxNhit[iH]->GetBinError(bin2),2);
	      }
	      err = sqrt(err);*/
	      maxNhit[iH]->SetBinContent(bin,integral/inttot);
	      maxNhit[iH]->SetBinError(bin,err/inttot);
	    }
	  }
	}


      mycL->cd();
      TLegend *leg = new TLegend(0.66,0.66,0.89,0.89);
      leg->SetFillColor(10);
      for (unsigned iH(0); iH<nHists;++iH){
	gStyle->SetOptStat(0);
	myc->cd();
	gPad->SetLogy(1);
	gStyle->SetStatX(0.89);
	gStyle->SetStatY(0.89);
	gStyle->SetStatH(0.5);
	gStyle->SetStatW(0.3);
	//maxNhit[iH]->Rebin(genEn[iE]>500? 20*iE : (iE>0?2*iE:1));
	maxNhit[iH]->GetXaxis()->SetRange(maxNhit[iH]->FindFirstBinAbove(0),25);//maxNhit[iH]->FindLastBinAbove(0));
	//maxNhit[iH]->SetTitle(";#hits with E_{hit} > E^{avg}_{max};Events");
	//maxNhit[iH]->SetTitle(";#hits with E_{hit} > E^{avg}_{max};Events");
	maxNhit[iH]->SetLineColor(iH+1);
	maxNhit[iH]->SetMarkerColor(iH+1);
	maxNhit[iH]->SetMarkerStyle(iH+21);

	//maxNhit[iH]->Rebin(2);
	maxNhit[iH]->Draw();
	
	//sprintf(buf,"Maximum = %3.0f MIPs",maxNhit[iH]->GetXaxis()->GetBinCenter(maxNhit[iH]->FindLastBinAbove(0)));
	lat.SetTextSize(0.05);
	//lat.DrawLatexNDC(0.12,0.85,buf);
	sprintf(buf,"%s #mu m Si",pSuffix.c_str());
	lat.DrawLatexNDC(0.12,0.75,buf);
	//sprintf(buf,"Geant4 Standalone, e-, pad size = %s",label[iH].c_str());
	//if (nGenEn==1) 
	sprintf(buf,"Geant4 Standalone, e- %d GeV, pad size = %s",genEn[iE], label[iH].c_str());
	lat.DrawLatexNDC(0.1,0.92,buf);
	
	std::ostringstream lsave;
	lsave << plotDir << pSuffix << "um/" << basename << label[iH];
	//if (nGenEn==1) 
	lsave << "e" << genEn[iE];
	if (integrate) lsave << "_int";
	lsave << ".pdf";
	myc->Print(lsave.str().c_str());

	mycL->cd(iE+1);
	gPad->SetLogy(1);
	//gPad->SetGridx();
	gPad->SetGridy();
	maxNhit[iH]->GetYaxis()->SetTitle("Event fraction");
	maxNhit[iH]->GetXaxis()->SetTitle("N_{hits}(E > E_{sat})");
	//maxNhit[iH]->DrawNormalized(iH==0?"PE":"PEsame");
	maxNhit[iH]->Sumw2();
	if (!integrate && maxNhit[iH]->Integral(0,maxNhit[iH]->GetNbinsX()+1)>0) maxNhit[iH]->Scale(1./maxNhit[iH]->Integral(0,maxNhit[iH]->GetNbinsX()+1));
	maxNhit[iH]->SetMaximum(1);
	maxNhit[iH]->SetMinimum(0.001);
	//maxNhit[iH]->GetXaxis()->SetRangeUser(,20);
	maxNhit[iH]->Draw(iH==0?"PE":"PEsame");
	leg->AddEntry(maxNhit[iH],label[iH].c_str(),"P");
	lat.SetTextSize(0.05);
	//sprintf(buf,"Maximum = %3.0f MIPs",maxNhit[iH]->GetXaxis()->GetBinCenter(maxNhit[iH]->FindLastBinAbove(0)));
	//lat.DrawLatexNDC(0.2,0.8,buf);
	if (iH==0){
	  //sprintf(buf,"%s #mu m Si",pSuffix.c_str());
	  //lat.DrawLatexNDC(0.15,0.3,buf);
	  //sprintf(buf,"Geant4 Standalone, e-, %s #mu m Si",pSuffix.c_str());
	  //if (nGenEn==1) 
	  sprintf(buf,"Geant4 Standalone, e- %d GeV, %s #mu m Si",genEn[iE],pSuffix.c_str());
	  lat.DrawLatexNDC(0.1,0.92,buf);
	}
      }//loop on hists

      mycL->cd(iE+1);
      if (genEn[iE]!=2000) lat.DrawLatexNDC(0.5,0.6,"Saturation: 3000 MIPs");
      else lat.DrawLatexNDC(0.11,0.6,"Saturation: 3000 MIPs");

      leg->Draw("same");
      }//loop on energies

      std::ostringstream lsave;
      lsave << plotDir << pSuffix<< "um/" << basename ;
      //if (nGenEn==1) lsave << "e" << genEn[iE];
      if (integrate) lsave << "int_";
      lsave << "all.pdf";
      mycL->Print(lsave.str().c_str());


      /*      TGraphErrors *gr[nHists][1];
      TLegend *leg[1];
      for (unsigned i(0); i<1;++i){
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
	  //else gr[iH][i] = new TGraphErrors(nGenEn,genEnD,mean3rms[iH],genEnErr,mean3rmserr[iH]); 
	  gr[iH][i]->GetXaxis()->SetTitle("Beam energy (GeV)");
	  if (i==0) gr[iH][i]->GetYaxis()->SetTitle("N(E>maxE)");
	  //else gr[iH][i]->GetYaxis()->SetTitle("<maxN>+3#times RMS (MIPs)");
	  gr[iH][i]->SetTitle("");
	  gr[iH][i]->SetMarkerColor(iH+1);
	  gr[iH][i]->SetLineColor(iH+1);
	  gr[iH][i]->SetMarkerStyle(20+iH);
	  gr[iH][i]->SetMaximum(20000);
	  gr[iH][i]->Draw(iH>0?"P":"AP");
	  leg[i]->AddEntry(gr[iH][i],label[iH].c_str(),"P");
	  /*
	  if (i==0){
	    gr[iH][i]->Fit("pol2","","same");
	    fitfunc[iH] = gr[iH][i]->GetFunction("pol2");
	  }
	  else {
	    gr[iH][i]->Fit("pol1","","same",20,2000);
	    fitfunc[iH] = gr[iH][i]->GetFunction("pol1");
	  }
	  fitfunc[iH]->SetLineColor(iH+1);
	  
	  if (i==0) sprintf(buf,"<maxN> = (%3.0f #pm %3.0f)  + (%3.1f #pm %3.1f)#times E + (%1.0e#pm%1.0e) #times E^{2}",fitfunc[iH]->GetParameter(0),fitfunc[iH]->GetParError(0),fitfunc[iH]->GetParameter(1),fitfunc[iH]->GetParError(1),fitfunc[iH]->GetParameter(2),fitfunc[iH]->GetParError(2));
	  else sprintf(buf,"<maxN> = (%3.1f #pm %3.1f)  + (%3.2f #pm %3.2f)#times E",fitfunc[iH]->GetParameter(0),fitfunc[iH]->GetParError(0),fitfunc[iH]->GetParameter(1),fitfunc[iH]->GetParError(1));
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
	lsave << plotDir << pSuffix << "um/NhitsAboveMax";
	//if (i==1) lsave << "_meanplus3rms";
	lsave << ".pdf";
	myc->Update();
	myc->Print(lsave.str().c_str());

	//return 1;
	}
      */
    }//loop on scenarios
  }//loop on versions
  
  
  return 0;
}
