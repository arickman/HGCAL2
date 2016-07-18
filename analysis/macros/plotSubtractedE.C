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

int plotSubtractedE(){

  SetTdrStyle();

  std::string filename = "PLOTS/PuSubtraction.root";
  TFile *f = TFile::Open(filename.c_str());
  if (!f) {
    std::cout << " -- Error, cannot open input root file " << filename << ". Exiting." << std::endl;
    return 1;
  }

  const unsigned nPu = 1;
  unsigned pu[nPu] = {140};

  const unsigned neta = 7;
  unsigned eta[neta]={17,19,21,23,25,27,29};

  const unsigned nSR = 8;

  const unsigned nCanvas = 3+neta;  
  TCanvas *myc[nCanvas];
  for (unsigned iC(0);iC<nCanvas;++iC){
    std::ostringstream lName;
    lName << "myc" << iC;
    myc[iC] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1500,1000);
    //myc[iC]->Divide(1,2);
    if (iC>0) myc[iC]->Divide(4,2);
  }
  myc[0]->Divide(1,2);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  for (unsigned ipu(0); ipu<nPu; ++ipu){//loop on pu

    TGraphErrors *grmean[neta];
    TGraphErrors *grsigma[neta];

    //TH1F *hist[neta][nSR][nGenEn];

    TLegend *leg = new TLegend(0.15,0.6,0.25,0.94);
    leg->SetFillColor(10);
    char buf[500];
    TLatex lat;

    double maxmean = 0;
    double maxsig = 0;

    for (unsigned ieta(0);ieta<neta;++ieta){//loop on eta
      
      TString srStr = "eta";
      srStr += eta[ieta];
      grmean[ieta] = new TGraphErrors();
      grmean[ieta]->SetName("grmean"+srStr);
      grmean[ieta]->SetTitle(";SR; <puE-E> (GeV)");
      grmean[ieta]->SetMarkerStyle(20+ieta);
      grmean[ieta]->SetMarkerColor(1+ieta);
      grmean[ieta]->SetLineColor(1+ieta);
      grsigma[ieta] = (TGraphErrors *) grmean[ieta]->Clone("grsigma"+srStr);
      grsigma[ieta]->SetTitle(";SR; #sigma(puE-E) (GeV)");
      sprintf(buf,"#eta=%3.1f",eta[ieta]/10.);
      leg->AddEntry(grmean[ieta],buf,"P");
      
      for (unsigned iSR(0); iSR<nSR;++iSR){//loop on signal region
	std::ostringstream ldir;
	ldir << "eta"<< eta[ieta]  << "_pu" << pu[ipu];
	f->cd(ldir.str().c_str());
	
	std::ostringstream label;
	label << "p_sigma_" << iSR << "_" << ldir.str();
	TH1F* hist = (TH1F*)gDirectory->Get(label.str().c_str());
	  
	if (!hist) {
	  std::cout << " -- Error, could not find hist " << label.str() << ". Continue..." << std::endl;
	  //return 1;
	  continue;
	}
	label.str("");
	label << "p_subtrEvspuE_" << iSR << "_" << ldir.str();
	TH2F* hist2d = (TH2F*)gDirectory->Get(label.str().c_str());

	myc[ieta+1]->cd(iSR+1);
	gPad->SetLogy(1);
	gPad->SetGridx(1);
	hist->Rebin(4);

	//TF1 *fitResult = new TF1("fitResult","[0]*TMath::Gaus(x,[1],[2],0)",hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
	//double maxHeight = hist->GetBinContent(hist->GetMaximumBin());
	//double minFit =  hist->GetXaxis()->GetBinCenter(hist->FindFirstBinAbove(maxHeight/50.));
	//double maxFit =  hist->GetXaxis()->GetBinCenter(hist->FindLastBinAbove(maxHeight/50.));
	
	//hist->GetXaxis()->SetRangeUser(minFit,maxFit);

	//double maxbin = hist->GetXaxis()->GetBinCenter(hist->GetMaximumBin());
	//fitResult->SetParameters(maxHeight,
	//			 maxbin,
	//			 (maxFit-minFit)/2.);
	hist->Draw();
	//hist->Fit("fitResult","+","same",0,maxFit);

	sprintf(buf,"#eta=%3.1f, PU %d",eta[ieta]/10.,pu[ipu]);
	lat.DrawLatexNDC(0.4,0.965,buf);
	sprintf(buf,"SR %d",iSR);
	lat.DrawLatexNDC(0.2,0.87,buf);
	sprintf(buf,"Mean = %3.3f +/- %3.3f GeV",hist->GetMean(),hist->GetMeanError());
	lat.DrawLatexNDC(0.2,0.8,buf);
	sprintf(buf,"RMS = %3.3f +/- %3.3f GeV",hist->GetRMS(),hist->GetRMSError());
	lat.DrawLatexNDC(0.2,0.73,buf);

	//sprintf(buf,"<Efit> = %3.3f +/- %3.3f GeV",fitResult->GetParameter(1),fitResult->GetParError(1));
	//lat.DrawLatexNDC(0.2,0.66,buf);
	//sprintf(buf,"RMSfit = %3.3f +/- %3.3f GeV",fitResult->GetParameter(2),fitResult->GetParError(2));
	//lat.DrawLatexNDC(0.2,0.59,buf);
	//sprintf(buf,"#chi^{2}/N = %3.3f/%d",fitResult->GetChisquare(),fitResult->GetNDF());
	//lat.DrawLatexNDC(0.2,0.52,buf);
  
	if (iSR==2){
	  myc[8]->cd(ieta+1);
	  gPad->SetLogy(1);
	  gPad->SetGridx(1);
	  hist->DrawCopy();
	  sprintf(buf,"#eta=%3.1f, PU %d",eta[ieta]/10.,pu[ipu]);
	  lat.DrawLatexNDC(0.4,0.965,buf);
	  sprintf(buf,"SR %d",iSR);
	  lat.DrawLatexNDC(0.2,0.87,buf);
	  sprintf(buf,"Mean = %3.3f +/- %3.3f GeV",hist->GetMean(),hist->GetMeanError());
	  lat.DrawLatexNDC(0.2,0.8,buf);
	  sprintf(buf,"RMS = %3.3f +/- %3.3f GeV",hist->GetRMS(),hist->GetRMSError());
	  lat.DrawLatexNDC(0.2,0.73,buf);

	  myc[9]->cd(ieta+1);
	  gPad->SetLogz(1);
	  gPad->SetGridx(1);
	  gPad->SetGridy(1);
	  hist2d->RebinX(2);
	  hist2d->RebinY(2);
	  hist2d->GetYaxis()->SetRangeUser(0,8);
	  hist2d->GetXaxis()->SetRangeUser(-10,20);
	  hist2d->Draw("colz");
	  sprintf(buf,"#eta=%3.1f, PU %d",eta[ieta]/10.,pu[ipu]);
	  lat.DrawLatexNDC(0.4,0.965,buf);
	  sprintf(buf,"SR %d",iSR);
	  lat.DrawLatexNDC(0.2,0.87,buf);
	  sprintf(buf,"Mean x = %3.2f +/- %3.2f GeV",hist2d->GetMean(1),hist2d->GetMeanError(1));
	  lat.DrawLatexNDC(0.2,0.8,buf);
	  sprintf(buf,"Mean y = %3.2f +/- %3.2f GeV",hist2d->GetMean(2),hist2d->GetMeanError(2));
	  lat.DrawLatexNDC(0.2,0.73,buf);
	  sprintf(buf,"RMS x = %3.2f +/- %3.2f GeV",hist2d->GetRMS(1),hist2d->GetRMSError(1));
	  lat.DrawLatexNDC(0.2,0.66,buf);
	  sprintf(buf,"RMS y = %3.2f +/- %3.2f GeV",hist2d->GetRMS(2),hist2d->GetRMSError(2));
	  lat.DrawLatexNDC(0.2,0.59,buf);
	}

	Int_t np=grmean[ieta]->GetN();
	if (iSR<nSR-1){
	  grmean[ieta]->SetPoint(np,iSR,hist->GetMean());
	  grmean[ieta]->SetPointError(np,0.0,hist->GetMeanError());
	  grsigma[ieta]->SetPoint(np,iSR,hist->GetRMS());
	  grsigma[ieta]->SetPointError(np,0.0,hist->GetRMSError());
	  
	  if (hist->GetMean()>maxmean) maxmean = hist->GetMean();
	  if (hist->GetRMS()>maxsig) maxsig = hist->GetRMS();
	}
      
      }//loop on SR
    }//loop on eta

    //save
    for (unsigned ieta(1);ieta<neta+1;++ieta){//loop on eta
      std::ostringstream lsave;
      lsave << "PLOTS/SubtractedE_pu" << pu[ipu] << "_eta" << eta[ieta-1] << ".pdf";
      myc[ieta]->Print(lsave.str().c_str());
    }
    std::ostringstream lsave;
    lsave << "PLOTS/SubtractedE_pu" << pu[ipu] << "_SR2.pdf";
    myc[8]->Print(lsave.str().c_str());
    lsave.str("");
    lsave << "PLOTS/SubtractedEvsPuE_pu" << pu[ipu] << "_SR2.pdf";
    myc[9]->Print(lsave.str().c_str());

    //draw
    for (unsigned ieta(0);ieta<neta;++ieta){//loop on eta
      
      Int_t np=grmean[ieta]->GetN();
      myc[0]->cd(1);
      grmean[ieta]->SetMinimum(0);
      grmean[ieta]->SetMaximum(maxmean*1.1);
      std::cout << " -- eta=" << eta[ieta] << " npoints = " << np 
		<< " min=" << grmean[ieta]->GetMinimum() << " max=" << grmean[ieta]->GetMaximum()
		<< std::endl;
      grmean[ieta]->Draw((ieta==0)?"APL":"PLsame");
      myc[0]->cd(2);
      //grsigma[ieta]->SetMaximum(grsigma[neta-1]->GetMaximum()*1.1);
      grsigma[ieta]->SetMinimum(0);
      grsigma[ieta]->SetMaximum(maxsig*1.1);
      grsigma[ieta]->Draw((ieta==0)?"APL":"PLsame");
      
    }//loop on eta

    myc[0]->cd(1);
    leg->Draw("same");
    myc[0]->cd(2);
    leg->Draw("same");
    
    lsave.str("");
    lsave << "PLOTS/MeanSigmaSubtr_pu" << pu[ipu] << ".pdf";
    myc[0]->Print(lsave.str().c_str());
    
    
  }//loop on pu

  return 0;
}//main
