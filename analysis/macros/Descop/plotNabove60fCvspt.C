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

#include "../TDRStyle.h"

#define ISPI

int plotNabove60fCvspt(){//main

  bool doresol = false;
  bool docm = false;

  const unsigned nT = 3;
  unsigned thresh[nT] = {30,45,60};
  unsigned threshMip[nT];
  for (unsigned iT(0); iT<nT;++iT) {
    threshMip[iT] = static_cast<unsigned>(26*thresh[iT]/60.+0.5);
  }

  SetTdrStyle();
#ifdef ISPI
  bool isPi = true;
  std::string particle = "pi-";
#else
  bool isPi = false;
  std::string particle = "gamma";
#endif

#ifdef ISPI
  const unsigned npt = 15;
  const double pt[npt] = {1,2,3,4,5,6,7,8,9,10,12,14,16,18,20};
#else
  const unsigned npt = 9;
  const double pt[npt] = {3,5,10,20,30,50,70,100,200};
#endif

#ifdef ISPI
  const unsigned neta = 5;
  const double eta[neta] = {1.75,2.0,2.25,2.5,2.75};
#else
  const unsigned neta = 6;
  const double eta[neta] = {1.7,1.9,2.1,2.3,2.5,2.7};
#endif

  const std::string path = "/afs/cern.ch/work/a/amagnan/PFCalEEAna//HGCalTime/gitV05-02-04/"+particle+"/";

  const unsigned nLim = 6;
  const unsigned limit[nLim] = {1,2,5,10,15,20};

  const unsigned nV = 3;
#ifdef ISPI
  unsigned v[nV] = {33,36,37};
#else
  unsigned v[nV] = {30,34,35};
#endif

  TCanvas *myc2 = new TCanvas("myc2","myc2",1500,1000);
  TCanvas *mycl = new TCanvas("mycl","mycl",1);
  TCanvas *myc[neta];
  TCanvas *myc1[neta];

  for (unsigned ie(0); ie<neta;++ie){//loop on version
    std::ostringstream lname;
    lname << "myc_" << eta[ie]; 
    myc[ie] = new TCanvas(lname.str().c_str(),lname.str().c_str(),1500,1000);
    myc[ie]->Divide(3,2);
    lname.str("");
    lname << "myc1_" << eta[ie]; 
    myc1[ie] = new TCanvas(lname.str().c_str(),lname.str().c_str(),1500,1000);
  }

  for (unsigned iT(0); iT<nT;++iT) {

    std::ostringstream label;
    TGraphErrors *grDummy[nLim];
    for (unsigned il(0);il<nLim;++il){//loop on limits
      grDummy[il] = new TGraphErrors();
      grDummy[il]->SetMinimum(0.001);
      grDummy[il]->SetMaximum(1.2);
      for (unsigned ipt(0);ipt<npt;++ipt){//loop on pt
	grDummy[il]->SetPoint(ipt,pt[ipt],0);
      }
      label.str("");
      label << ";p_{T} (GeV);proba(n_{E>" << thresh[iT] << "fC}#geq" << limit[il] << ")";
      grDummy[il]->SetTitle(label.str().c_str());
      
    }//loop on limits
    TGraphErrors *grDummy1 = new TGraphErrors();
    if (!doresol){
      grDummy1->SetMinimum(isPi ?0 : 0.5);
      grDummy1->SetMaximum(isPi ? 3 : 1000);
  } else if (!docm) {
    grDummy1->SetMinimum(isPi ? 0 : 0.5);
    grDummy1->SetMaximum(isPi ? 3 : 50);
  } else {
    grDummy1->SetMinimum(isPi ? 0 : 0.);
    grDummy1->SetMaximum(isPi ? 3 :1.5);
  }
  for (unsigned ipt(0);ipt<npt;++ipt){//loop on pt
    grDummy1->SetPoint(ipt,pt[ipt],0);
  }
  label.str("");
  if (!doresol) label << ";p_{T} (GeV);<n_{E>" << thresh[iT] << "fC}>";
  else if (!docm) label << ";p_{T} (GeV);<50/#sqrt{n_{E>" << thresh[iT] << "fC}}> (ps)";
  else label << ";p_{T} (GeV);<50#times c/#sqrt{n_{E>" << thresh[iT] << "fC}}> (cm)";
  grDummy1->SetTitle(label.str().c_str());

  TGraphErrors *gr[nV][neta][nLim];
  TGraphErrors *grMean[nV][neta];
  TTree *tree = 0;
  TFile *ftmp = 0;

  TLatex lat;
  char buf[200];

  TLegend *leg;
  if (!doresol) leg = new TLegend(0.68,0.15,0.92,0.4);
  else leg = new TLegend(0.68,0.7,0.92,0.94);
  leg->SetFillColor(10);
  for (unsigned ieta(0); ieta<neta;++ieta){//loop on eta
    for (unsigned il(0);il<nLim;++il){//loop on limits
      myc[ieta]->cd(il+1);
      grDummy[il]->Draw("AL");
    }
    myc1[ieta]->cd();
    if (!isPi) gPad->SetLogy(1);
    gPad->SetGridy(1);
    grDummy1->Draw("AL");
    grDummy1->GetYaxis()->SetNdivisions(121);
    
    for (unsigned iv(0); iv<nV;++iv){//loop on version
      for (unsigned il(0);il<nLim;++il){//loop on limits
	std::ostringstream lname;
	lname << "grv" << v[iv]<< "_eta" << eta[ieta] << "_lim" << limit[il];
	gr[iv][ieta][il] = new TGraphErrors();
	gr[iv][ieta][il]->SetName(lname.str().c_str());
      }
      grMean[iv][ieta] = new TGraphErrors();
      std::ostringstream lname;
      lname << "grv" << v[iv]<< "_eta" << eta[ieta] << "_mean";
      grMean[iv][ieta]->SetName(lname.str().c_str());
      //set dummy values to have full range for sure...
      for (unsigned ipt(0);ipt<npt;++ipt){//loop on pt
	std::ostringstream fname;
	fname << path << "v" << v[iv] << "_et" << pt[ipt] << "_eta" << eta[ieta];
	if (eta[ieta]==2) fname << ".0";
	if (thresh[iT]!=60) fname << "_thresh" << threshMip[iT];
	if (ieta==4 && (pt[ipt]==50 || pt[ipt]==70)) fname << "/NAbove" << threshMip[iT] << "_old.root";
	else fname  << "/NAbove" << threshMip[iT] << ".root";
	ftmp = TFile::Open(fname.str().c_str());
	if (!ftmp) {
	  continue;
	}
	ftmp->cd();
	tree = (TTree*)gDirectory->Get("outtree");
	if (!tree) {
	  continue;
	}
	const unsigned nbins = 2000;
	fname.str("");
	fname << ";n(E>" << threshMip[iT] << " mips);proba";
	TH1F *hprob = new TH1F("hprob",fname.str().c_str(),nbins,0,nbins);
	//hprob->StatOverflows(1);
	hprob->Sumw2();
	myc2->cd();
	if (!doresol) tree->Draw("nover>>hprob");
	else if (!docm) tree->Draw("50./sqrt(nover)>>hprob","nover>=1");
	else tree->Draw("50.*0.03/sqrt(nover)>>hprob","nover>=1");
	hprob->GetXaxis()->SetRangeUser(0,1000);
	unsigned nentries = hprob->GetEntries();
	std::cout << "v " << v[iv] 
		  << " et " << pt[ipt] 
		  << " eta " << eta[ieta];
	if (nentries>0) {
	  std::cout << " nentries " << nentries << " mean " << hprob->GetMean() << " " << hprob->GetMeanError() ;

	  //if (!doresol) {
	  double xval = pt[ipt];
	  double yval = 0;
	  if (iv==0){
	    if (grMean[0][ieta]->GetN()>0) grMean[0][ieta]->GetPoint(grMean[iv][ieta]->GetN(),xval,yval);
	    else yval = hprob->GetMean();
	  }
	  grMean[iv][ieta]->SetPoint(grMean[iv][ieta]->GetN(),pt[ipt],hprob->GetMean()/yval);
	  grMean[iv][ieta]->SetPointError(grMean[iv][ieta]->GetN()-1,0,hprob->GetMeanError()/yval);
	    //}
	    //else if (!docm) {
	    //grMean[iv][ieta]->SetPoint(grMean[iv][ieta]->GetN(),pt[ipt],50./sqrt(hprob->GetMean()));
	    //grMean[iv][ieta]->SetPointError(grMean[iv][ieta]->GetN()-1,0,0.5*50*hprob->GetMeanError()*pow(hprob->GetMean(),-3./2));
	    //}	  
	    //else {
	    //grMean[iv][ieta]->SetPoint(grMean[iv][ieta]->GetN(),pt[ipt],50.*0.03/sqrt(hprob->GetMean()));
	    //grMean[iv][ieta]->SetPointError(grMean[iv][ieta]->GetN()-1,0,0.5*50*0.03*hprob->GetMeanError()*pow(hprob->GetMean(),-3./2));
	    //}
	  hprob->Scale(1./nentries);
	}
	else {
	  continue;
	}
	for (unsigned il(0);il<nLim;++il){//loop on limits
	  double error = 0;
	  double integral = hprob->IntegralAndError(limit[il]+1,nbins+1,error);
 
	  std::cout << " entries=" << hprob->GetEntries() 
		    << " integral(" << limit[il]+1 << ",N) " << integral << "+/-" << error
		    << std::endl;
	  gr[iv][ieta][il]->SetPoint(gr[iv][ieta][il]->GetN(),pt[ipt],integral);
	  gr[iv][ieta][il]->SetPointError(gr[iv][ieta][il]->GetN()-1,0,error);
	}
	hprob->Delete();
	ftmp->Close();
      }//loop on pt
      for (unsigned il(0);il<nLim;++il){//loop on limits
	myc[ieta]->cd(il+1);
	std::cout << "v " << v[iv] 
		  << " eta " << eta[ieta]
		  << " limit " << limit[il]
		  << " npoints = " << gr[iv][ieta][il]->GetN()
		  << " nmean = " << grMean[iv][ieta]->GetN()
		  << std::endl;
	
	gr[iv][ieta][il]->SetMarkerColor(iv+1);
	gr[iv][ieta][il]->SetLineColor(iv+1);
	gr[iv][ieta][il]->SetMarkerStyle(20+iv);
	gr[iv][ieta][il]->SetMinimum(0.0);
	gr[iv][ieta][il]->SetMaximum(1.1);
	gr[iv][ieta][il]->Draw("PE");
      }//loop on limits

      myc1[ieta]->cd();
      grMean[iv][ieta]->SetMarkerColor(iv+1);
      grMean[iv][ieta]->SetLineColor(iv+1);
      grMean[iv][ieta]->SetMarkerStyle(20+iv);
      grMean[iv][ieta]->SetMinimum(0.01);
      grMean[iv][ieta]->SetMaximum(25);
      grMean[iv][ieta]->Draw("PEL");

    }//loop on versions

    std::cout << " End loopon versions" << std::endl;

    for (unsigned il(0);il<nLim;++il){//loop on limits
 
      myc[ieta]->cd(il+1);
      if (il>2) gPad->SetLogy(1);
      if (ieta==0 && il==0) {
	leg->AddEntry(gr[0][ieta][il],"TP-28-12","P");
	if (nV>1) leg->AddEntry(gr[1][ieta][il],"TP-24-11","P");
	if (nV>2) leg->AddEntry(gr[2][ieta][il],"TP-18-09","P");
      }
      
      if (il==0) leg->Draw("same");
      if (isPi) sprintf(buf,"#pi^{-}, #eta = %3.2f",eta[ieta]);
      else sprintf(buf,"#gamma, #eta = %3.2f",eta[ieta]);
      lat.SetTextSize(0.06);
      if (il==0) lat.DrawLatexNDC(0.5,0.85,buf);
      lat.SetTextSize(0.03);
      if (il==3) lat.DrawLatexNDC(0.02,0.02,"HGCAL Standalone Simulation");

      if (!doresol && !docm) {
	mycl->cd();
	grDummy[il]->Draw("AL");
	for (unsigned iv(0); iv<nV;++iv){//loop on version
	  gr[iv][ieta][il]->Draw("PE");
	}
	leg->Draw("same");
	lat.SetTextSize(0.06);
	lat.DrawLatexNDC(0.5,0.85,buf);
	lat.SetTextSize(0.03);
	lat.DrawLatexNDC(0.02,0.02,"HGCAL Standalone Simulation");
	std::ostringstream lsave;
	lsave << "SummaryNabove" << thresh[iT] << "fC_" << particle << "_vspt_eta" << eta[ieta] << "_limit" << limit[il];
	mycl->Update();
	mycl->Print((lsave.str()+".pdf").c_str());
	mycl->Print((lsave.str()+".C").c_str());
      }


    }
    myc[ieta]->Update();
    std::ostringstream lsave;
    lsave << "SummaryNabove"<< thresh[iT] << "fC_" << particle << "_vspt_eta" << eta[ieta];
    if (!doresol && !docm) {
      myc[ieta]->Print((lsave.str()+".pdf").c_str());
      myc[ieta]->Print((lsave.str()+".C").c_str());
    }

    myc1[ieta]->cd();
    leg->Draw("same");
    lat.SetTextSize(0.06);
    lat.DrawLatexNDC(0.2,0.85,buf);
    lat.SetTextSize(0.03);
    lat.DrawLatexNDC(0.02,0.02,"HGCAL Standalone Simulation");
    myc1[ieta]->Update();
    lsave.str("");
    if (!doresol) lsave << "SummaryNMeanabove" << thresh[iT] << "fC_" << particle << "_vspt_eta" << eta[ieta];
    else if (!docm) lsave << "SummaryResolabove" << thresh[iT] << "fC_" << particle << "_vspt_eta" << eta[ieta];
    else lsave << "SummaryResolAtVertexabove" << thresh[iT] << "fC_" << particle << "_vspt_eta" << eta[ieta] ;
    myc1[ieta]->Print((lsave.str()+".pdf").c_str());
    myc1[ieta]->Print((lsave.str()+".C").c_str());


  }//loop on eta

  }//loop on thresholds
  return 0;
}//main
