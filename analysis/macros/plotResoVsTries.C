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
struct FitResult{
  double chi2;
  unsigned ndf;
  double mean;
  double sigma;
  double meanerr;
  double sigmaerr;
  double reso(){
    return sigma/mean;
  };
};

unsigned fitEnergy(TH1F *hist,
		   TPad *pad,
		   std::string unitStr,
		   FitResult & lres,
		   unsigned isr){
  
  pad->cd();
  //double eMin = hist->GetMean()-5*hist->GetRMS();
  //double eMax = hist->GetMean()+5*hist->GetRMS();
  //hist->GetXaxis()->SetRangeUser(eMin,eMax);
  hist->Draw("PE");


  double nRMSm = isr<1? 1 : 2;
  double nRMSp = 2;
  
  TF1 *fitResult = new TF1("fitResult","[0]*TMath::Gaus(x,[1],[2],0)",hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
  fitResult->SetParameters(hist->Integral(),
			   hist->GetXaxis()->GetBinCenter(hist->GetMaximumBin()),
			   hist->GetRMS());

  //std::cout << " Initial params: "  << fitResult->GetParameter(0) << " "<< fitResult->GetParameter(1) << " " << fitResult->GetParameter(2)
  //<< std::endl;


  int status = hist->Fit("fitResult","L0QEMI","",
			 fitResult->GetParameter(1)-nRMSm*fitResult->GetParameter(2),
			 fitResult->GetParameter(1)+nRMSp*fitResult->GetParameter(2));
  
  
  //std::cout << " First fit: " << status << " " << fitResult->GetParameter(1) << " " << fitResult->GetParameter(2)
  //<< std::endl;

  // if ((status != 0 && status != 4000) || fitResult->GetChisquare()/fitResult->GetNDF()>5){
  //   //std::cout << " -- Bad fit ! Try again..." << std::endl;
  //   status = hist->Fit("fitResult","L0QEMI","",
  // 		       fitResult->GetParameter(1)-nRMSm*fitResult->GetParameter(2),
  // 		       fitResult->GetParameter(1)+nRMSp*fitResult->GetParameter(2));
    
  //   std::cout << " Second fit: " << status << " " << fitResult->GetParameter(1) << " " << fitResult->GetParameter(2)
  // 	      << std::endl;
  // }
  
  // std::cout << " Final fit: " << fitResult->GetParameter(1) << " " << fitResult->GetParameter(2)
  // 	    << std::endl;
  
  fitResult->SetLineColor(2);
  fitResult->Draw("same");
  
  if (status != 0 && status != 4000) {
    std::cout << " Warning! Fit failed with status " << status << "! Please have a look at the verbose output below...." << std::endl;
    hist->Fit("fitResult","L0EMI","",
	      fitResult->GetParameter(1)-nRMSm*fitResult->GetParameter(2),
	      fitResult->GetParameter(1)+nRMSp*fitResult->GetParameter(2));
    //totalE for pu140 is expected to be pathological :/
    //if (isr!=7) return 1;
  }

  // char buf[500];
  // TLatex lat;
  // double latx = hist->GetXaxis()->GetXmin()+(hist->GetXaxis()->GetXmax()-hist->GetXaxis()->GetXmin())/20.;
  // double laty = hist->GetMaximum();
  // sprintf(buf,"<Efit> = %3.3f +/- %3.3f %s",fitResult->GetParameter(1),fitResult->GetParError(1),unitStr.c_str());
  // lat.DrawLatex(latx,laty*0.9,buf);
  // sprintf(buf,"RMSfit = %3.3f +/- %3.3f %s",fitResult->GetParameter(2),fitResult->GetParError(2),unitStr.c_str());
  // lat.DrawLatex(latx,laty*0.8,buf);
  // sprintf(buf,"RMS/meanfit = %3.3f",fitResult->GetParameter(2)/fitResult->GetParameter(1));
  // lat.DrawLatex(latx,laty*0.7,buf);
  
  // sprintf(buf,"#chi^{2}/N = %3.3f/%d = %3.3f",fitResult->GetChisquare(),fitResult->GetNDF(),fitResult->GetChisquare()/fitResult->GetNDF());
  // lat.DrawLatex(latx,laty*0.6,buf);
  
  lres.chi2 = fitResult->GetChisquare();
  lres.ndf = fitResult->GetNDF();
  lres.mean = fitResult->GetParameter(1);
  lres.meanerr = fitResult->GetParError(1);
  lres.sigma = fitResult->GetParameter(2);
  lres.sigmaerr = fitResult->GetParError(2);

  return 0;
};


int plotResoVsTries(){

  SetTdrStyle();


  const unsigned nPts = 1000;
  const unsigned ICval = 10;
  const unsigned nE = 5;

  std::ostringstream label;
  TFile *fcalib[nPts];
  
  TCanvas *mycReso = new TCanvas("mycReso","mycReso",1500,1000);
  mycReso->Divide(2,3);
  TCanvas *mycSigma = new TCanvas("mycSigma","mycSigma",1500,1000);
  mycSigma->Divide(2,3);
  TCanvas *mycMean = new TCanvas("mycMean","mycMean",1500,1000);
  mycMean->Divide(2,3);
  TCanvas *mycC = new TCanvas("mycC","Linearity",1500,1000);
  TCanvas *mycE = new TCanvas("mycE","Efit",1500,1000);

  gStyle->SetOptFit(1111);
  gStyle->SetOptStat("eMR");

  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.5);
  gStyle->SetStatX(0.5);
  gStyle->SetStatY(0.5);

  TLatex lat;
  char buf[500];

  TH1F *mean[nE][2];
  TH1F *sigma[nE][2];
  TH1F *reso[nE][2];
  TH1F *diffMean = new TH1F("diffMean",";#DeltaE/E(400 GeV) - #DeltaE/E(20 GeV)",100,-0.1,0.1);
  TH2F *h400vs20_mean = new TH2F("h400vs20_mean",";#DeltaE/E(20 GeV);#DeltaE/E(400 GeV);trials",100,-0.5,0.5,100,-0.5,0.5);
  TH2F *h400vs20_sigma = new TH2F("h400vs20_sigma",";#Delta#sigma/#sigma(20 GeV);#Delta#sigma/#sigma(400 GeV);trials",100,-3,3,100,-3,3);
  TH2F *h400vs20_reso = new TH2F("h400vs20_reso",";#Delta Reso/Reso(20 GeV);#Delta Reso/Reso(400 GeV);trials",100,-3,3,100,-3,3);

  double genPt[nE] = {5,10,20,50,100};
  double genEn[nE];
   for (unsigned ie(0); ie<nE;++ie){
     genEn[ie] = genPt[ie]*cosh(2.1);
      for (unsigned ip(0); ip<2;++ip){
	label.str("");
	label << "mean_" << ie << "_" << ip;
	mean[ie][ip] = new TH1F(label.str().c_str(),";<E> (GeV);events",100,0.75*genEn[ie],1.25*genEn[ie]);
	label.str("");
	label << "sigma_" << ie << "_" << ip;
	sigma[ie][ip] = new TH1F(label.str().c_str(),";#sigma (GeV);events",100,0,15);
	label.str("");
	label << "reso_" << ie << "_" << ip;
	reso[ie][ip] = new TH1F(label.str().c_str(),";#sigma/<E>;events",100,0,0.06);
      }
   }
  for (unsigned ic(0);ic<nPts;++ic){//loop on intercalib
    label.str("");
    label << "PLOTS/CalibTree";
    label << "_SR7_IC" << ICval;
    label << "_try" << ic;
    label << ".root";
    fcalib[ic] = TFile::Open(label.str().c_str());
    if (!fcalib[ic]) {
      std::cout << " -- failed to open file: " << label.str() << std::endl;
      continue;
    }
    else {
      std::cout << " -- file " << label.str() << " successfully opened." << std::endl;
    }
    fcalib[ic]->cd();
    TH1F *Ereco[nE];
    FitResult lres[nE];
    FitResult lresSmear[nE];

    for (unsigned ie(0); ie<nE;++ie){
      label.str("");     
      label << "TreeEtot_" << genPt[ie];

      TTree *tree = (TTree*)gDirectory->Get(label.str().c_str());
      if (!tree) {
	std::cout << " -- Tree " << label.str() << " not found." << std::endl;
	return 1;
      }

      mycE->cd();
      tree->Draw("(Etot-calibOffset)/calibSlope","","");
      label.str("");
      label << "energy" << genEn[ie];
      Ereco[ie] = (TH1F*)(gPad->GetPrimitive("htemp"))->Clone(label.str().c_str());
      if (!Ereco[ie]){
	std::cout << " -- ERROR, pointer for histogram " << label.str() << " is null." << std::endl;
	return 1;
      }
      Ereco[ie]->SetTitle(";E (GeV);events");
      TPad *lpad = (TPad*)(mycE->cd());
      if (fitEnergy(Ereco[ie],lpad,"GeV",lres[ie],7)!=0) return 1;
      mean[ie][0]->Fill(lres[ie].mean);
      sigma[ie][0]->Fill(lres[ie].sigma);
      reso[ie][0]->Fill(lres[ie].sigma/lres[ie].mean);

      tree->Draw("(Esmear-calibOffset)/calibSlope","","");
      label.str("");
      label << "energySmear" << genEn[ie];
      Ereco[ie] = (TH1F*)(gPad->GetPrimitive("htemp"))->Clone(label.str().c_str());
      if (!Ereco[ie]){
	std::cout << " -- ERROR, pointer for histogram " << label.str() << " is null." << std::endl;
	return 1;
      }
      if (fitEnergy(Ereco[ie],lpad,"GeV",lresSmear[ie],7)!=0) return 1;
      mean[ie][1]->Fill(lresSmear[ie].mean);
      sigma[ie][1]->Fill(lresSmear[ie].sigma);
      reso[ie][1]->Fill(lresSmear[ie].sigma/lresSmear[ie].mean);
      
    }//loop on energies

    diffMean->Fill((lresSmear[nE-1].mean-lres[nE-1].mean)/lres[nE-1].mean - (lresSmear[0].mean-lres[0].mean)/lres[0].mean);
    h400vs20_mean->Fill((lresSmear[0].mean-lres[0].mean)/lres[0].mean,(lresSmear[nE-1].mean-lres[nE-1].mean)/lres[nE-1].mean);
    h400vs20_sigma->Fill((lresSmear[0].sigma-lres[0].sigma)/lres[0].sigma,(lresSmear[nE-1].sigma-lres[nE-1].sigma)/lres[nE-1].sigma);
    h400vs20_reso->Fill((lresSmear[0].reso()-lres[0].reso())/lres[0].reso(),(lresSmear[nE-1].reso()-lres[nE-1].reso())/lres[nE-1].reso());

    fcalib[ic]->Close();
  }//loop on tries

  TLegend *leg = new TLegend(0.75,0.2,0.92,0.4);
  leg->SetFillColor(10);
  gStyle->SetStatX(1.);
  gStyle->SetStatY(1.);

  for (unsigned ie(0); ie<nE;++ie){
    for (int ip(1); ip>=0;--ip){
      std::cout << " Processing ip=" << ip << std::endl;
      mycReso->cd(ie+1);
      reso[ie][ip]->SetLineColor(ip+1);
      //reso[ie][ip]->SetMaximum(1.1*reso[ie][0]->GetMaximum());
      reso[ie][ip]->Draw(ip==1?"":"same");
      if (ie==0) {
	if (ip==0) leg->AddEntry(reso[ie][ip],"No smearing","L");
	if (ip==1) leg->AddEntry(reso[ie][ip],"10% smearing","L");
      }
      if (ip==1) leg->Draw("same");

      lat.SetTextSize(0.1);
      sprintf(buf,"Single #gamma, #eta=2.1, E=%3.1f",genEn[ie]);
      lat.DrawLatexNDC(0.2,0.8,buf);
      sprintf(buf,"ICsmear = %d %%",ICval);
      lat.DrawLatexNDC(0.2,0.7,buf);
      //lat.DrawLatexNDC(0.01,0.01,"HGCAL G4 standalone");
      mycSigma->cd(ie+1);
      sigma[ie][ip]->SetLineColor(ip+1);
      //sigma[ie][ip]->SetMaximum(1.1*sigma[ie][0]->GetMaximum());
      sigma[ie][ip]->Draw(ip==1?"":"same");
      if (ip==1) leg->Draw("same");

      lat.SetTextSize(0.1);
      sprintf(buf,"Single #gamma, #eta=2.1, E=%3.1f",genEn[ie]);
      lat.DrawLatexNDC(0.2,0.8,buf);
      sprintf(buf,"ICsmear = %d %%",ICval);
      lat.DrawLatexNDC(0.2,0.7,buf);

      mycMean->cd(ie+1);
      mean[ie][ip]->SetLineColor(ip+1);
      //mean[ie][ip]->SetMaximum(1.1*mean[ie][0]->GetMaximum());
      mean[ie][ip]->Draw(ip==1?"":"same");
      if (ip==1) leg->Draw("same");
      //mean[ie][ip]->Fit("gaus","","same");

      lat.SetTextSize(0.1);
      sprintf(buf,"Single #gamma, #eta=2.1, E=%3.1f",genEn[ie]);
      lat.DrawLatexNDC(0.2,0.8,buf);
      sprintf(buf,"ICsmear = %d %%",ICval);
      lat.DrawLatexNDC(0.2,0.7,buf);
 
    }
  }

  mycReso->Update();
  mycReso->Print("PLOTS/ResoIC10.pdf");
  mycReso->Print("PLOTS/ResoIC10.C");
  mycMean->Update();
  mycMean->Print("PLOTS/MeanIC10.pdf");
  mycMean->Print("PLOTS/MeanIC10.C");
  mycSigma->Update();
  mycSigma->Print("PLOTS/SigmaIC10.pdf");
  mycSigma->Print("PLOTS/SigmaIC10.C");

  mycC->Divide(2,2);
  mycC->cd(1);
  diffMean->Draw();
  mycC->cd(2);
  h400vs20_mean->Draw("colz");
  mycC->cd(3);
  h400vs20_sigma->Draw("colz");
  mycC->cd(4);
  h400vs20_reso->Draw("colz");
  mycC->Update();
  mycC->Print("PLOTS/E400vs20IC10.pdf");
  mycC->Print("PLOTS/E400vs20IC10.C");

  return 0;

}//main
