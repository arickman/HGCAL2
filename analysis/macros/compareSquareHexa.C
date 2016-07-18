#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>

#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH2F.h"
#include "TProfile.h"
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

int compareSquareHexa(){//main

  const unsigned nLayers = 28;
  const unsigned nF = 2;
  TFile *fin[nF];
  fin[0] = TFile::Open("/afs/cern.ch/work/a/amagnan/PFCalEEAna/HGCalJack/gitV06-03-04/version100/model3/gamma/eta20_et20_pu0_IC3_Si2.root");
  fin[1] = TFile::Open("/afs/cern.ch/work/a/amagnan/PFCalEEAna/HGCalHexa/gittestHexa/version100/model3/gamma/eta20_et20_pu0_IC3_Si2.root");

  std::string label[nF] = {"Square","Hexagons"};
  
  TCanvas *myc = new TCanvas("myc","myc",1);
  TCanvas *mycSR2 = new TCanvas("mycSR2","mycSR2",1);
  TCanvas *myc1 = new TCanvas("myc1","myc1",1500,1000);
  TCanvas *myc17 = new TCanvas("myc17","myc17",1500,1000);
  TCanvas *myc7all = new TCanvas("myc7all","myc7all",1500,1000);
  TCanvas *mycmaxout7 = new TCanvas("mycmaxout7","mycmaxout7",1500,1000);
  TCanvas *myc2D = new TCanvas("myc2D","myc2D",1500,1000);
  TCanvas *mycprof = new TCanvas("mycprof","mycprof",1500,1000);
  TCanvas *mycCor = new TCanvas("mycCor","mycCor",1);

  const unsigned iSR = 2;

  TH1F *Etot[nF];
  TH1F *ESR2[nF];
  TTree *tree[nF];
  TH1F *E1[nF][nLayers];
  TH1F *E1over7[nF][nLayers];
  TH1F *E7overtot[nF][nLayers];
  TH1F *EmaxOut7[nF][nLayers];
  TH2F *transverseCor[nF][nLayers];
  TProfile *prof[nF][nLayers];
  TH1F *ECor[nF];

  std::ostringstream lname;
  std::ostringstream lvar;
  std::ostringstream lcut;
  myc->cd();

  for (unsigned iF(0); iF<nF;++iF){//loop on files
    fin[iF]->cd("Energies");
    tree[iF] = (TTree*)gDirectory->Get("Ereso");
    if (!tree[iF]) {
      std::cout << " Tree not found for file " << fin[iF]->GetName() << std::endl;
      return 1;
    }
    lname.str("");
    lname << "Etot" << label[iF];
    Etot[iF] = new TH1F(lname.str().c_str(),";Etot (mips);showers",100,6000,9500);
    Etot[iF]->SetLineColor(1+iF);
    Etot[iF]->SetMarkerColor(1+iF);
    Etot[iF]->SetMarkerStyle(21+iF);
    Etot[iF]->Sumw2();
    lvar.str("");
    lvar << "wgtEtotal>>" << lname.str();
    tree[iF]->Draw(lvar.str().c_str());

    lname.str("");
    lname << "ESR2" << label[iF];
    ESR2[iF] = new TH1F(lname.str().c_str(),";ESR2 (mips);showers",100,4700,7100);
    ESR2[iF]->SetLineColor(1+iF);
    ESR2[iF]->SetMarkerColor(1+iF);
    ESR2[iF]->SetMarkerStyle(21+iF);
    ESR2[iF]->Sumw2();
    lvar.str("");
    //lname << "(";
    for (unsigned iL(0);iL<nLayers;++iL){	      
      if (iL==0) lvar << "absweight_" << iL  << "*energy_" << iL << "_SR" << iSR ;
      else lvar << "+" << "absweight_" << iL  << "*energy_" << iL << "_SR" << iSR ;
    }
    lvar << ">>" << lname.str();
    tree[iF]->Draw(lvar.str().c_str());

    for (unsigned iL(0);iL<nLayers;++iL){	      
      lname.str("");
      lname << "E1" << label[iF] << "_" << iL;
      E1[iF][iL] = new TH1F(lname.str().c_str(),";ESR0 (mips);showers",200,0,500);
      E1[iF][iL]->SetLineColor(1+iF);
      E1[iF][iL]->SetMarkerColor(1+iF);
      E1[iF][iL]->SetMarkerStyle(21+iF);
      E1[iF][iL]->Sumw2();
      lvar.str("");
      lvar << "energy_" << iL << "_SR0>>" << lname.str();
      tree[iF]->Draw(lvar.str().c_str());

      lname.str("");
      lname << "E1over7" << label[iF] << "_" << iL;
      E1over7[iF][iL] = new TH1F(lname.str().c_str(),";ESR0/ESR2;showers",50,0,1.);
      E1over7[iF][iL]->SetLineColor(1+iF);
      E1over7[iF][iL]->SetMarkerColor(1+iF);
      E1over7[iF][iL]->SetMarkerStyle(21+iF);
      E1over7[iF][iL]->Sumw2();
      for (unsigned iS(iL);iS<nLayers;++iS){
	lvar.str("");
	lvar << "energy_" << iS << "_SR0/energy_" << iS << "_SR2>>+" << lname.str();
	lcut.str("");
	lcut << "energy_" << iS << "_SR0 > 5";
	if (iS>0) lcut << " && energy_" << iS-1 << "_SR0 < 5";
	tree[iF]->Draw(lvar.str().c_str(),lcut.str().c_str());
      }

      lname.str("");
      lname << "E7overtot" << label[iF] << "_" << iL;
      E7overtot[iF][iL] = new TH1F(lname.str().c_str(),";ESR2/Etot;showers",50,0,1.);
      E7overtot[iF][iL]->SetLineColor(1+iF);
      E7overtot[iF][iL]->SetMarkerColor(1+iF);
      E7overtot[iF][iL]->SetMarkerStyle(21+iF);
      E7overtot[iF][iL]->Sumw2();
      for (unsigned iS(iL);iS<nLayers;++iS){
	lvar.str("");
	lvar << "energy_" << iS << "_SR2/energy_" << iS << "_SR5>>+" << lname.str();
	lcut.str("");
	lcut << "energy_" << iS << "_SR0 > 5";
	if (iS>0) lcut << " && energy_" << iS-1 << "_SR0 < 5";
	tree[iF]->Draw(lvar.str().c_str(),lcut.str().c_str());
      }

      lname.str("");
      lname << "EmaxOut7" << label[iF] << "_" << iL;
      EmaxOut7[iF][iL] = new TH1F(lname.str().c_str(),";max Ehit outside SR2 (mips);showers",100,0,150.);
      EmaxOut7[iF][iL]->SetLineColor(1+iF);
      EmaxOut7[iF][iL]->SetMarkerColor(1+iF);
      EmaxOut7[iF][iL]->SetMarkerStyle(21+iF);
      EmaxOut7[iF][iL]->Sumw2();
      lvar.str("");
      lvar << "maxhitEoutside_" << iL << "_SR2>>+" << lname.str();
      lcut.str("");
      tree[iF]->Draw(lvar.str().c_str(),lcut.str().c_str());
      
      lname.str("");
      lname << "transverseCor" << label[iF] << "_" << iL;
      transverseCor[iF][iL] = new TH2F(lname.str().c_str(),";ESR0/ESR2;ESR2/Etot;showers",50,0.,1.,50,iL<14? 0.5 : 0,1.);
      transverseCor[iF][iL]->SetLineColor(1+iF);
      transverseCor[iF][iL]->SetMarkerColor(1+iF);
      transverseCor[iF][iL]->SetMarkerStyle(21+iF);
      transverseCor[iF][iL]->Sumw2();
      for (unsigned iS(iL);iS<nLayers;++iS){
	lvar.str("");
	lvar << "energy_" << iS << "_SR2/energy_" << iS << "_SR5:energy_" << iS << "_SR0/energy_" << iS << "_SR2>>+" << lname.str();
	lcut.str("");
	lcut << "energy_" << iS << "_SR0 > 5";
	if (iS>0) lcut << " && energy_" << iS-1 << "_SR0 < 5";
	tree[iF]->Draw(lvar.str().c_str(),lcut.str().c_str());
      }
    }//loop on layers

  }//loop on files


  TLegend *leg = new TLegend(0.63,0.8,0.99,0.99);
  leg->SetFillColor(10);

  gStyle->SetOptStat(0);

  TLatex lat;
  char buf[500];
  myc->cd();
  for (unsigned iF(0); iF<nF;++iF){//loop on files
    Etot[iF]->Scale(Etot[1]->GetEntries()/Etot[iF]->GetEntries());
    Etot[iF]->Draw(iF==0?"PE":"PEsame");
    leg->AddEntry(Etot[iF],label[iF].c_str(),"P");

    Etot[iF]->Fit("gaus","","same");
    TF1 *fit = (TF1*)Etot[iF]->GetFunction("gaus");
    if (!fit) continue;
    fit->SetLineColor(1+iF);
    fit->Draw("same");
    double reso = fit->GetParameter(2)/fit->GetParameter(1);
    //double err = reso*sqrt(pow(fit->GetParError(2)/fit->GetParameter(2),2)+pow(fit->GetParError(1)/fit->GetParameter(1),2));
    double err = fit->GetParError(2)/fit->GetParameter(1);
    sprintf(buf,"#sigma/E=%3.4f #pm %3.4f",reso,err);
    lat.SetTextColor(1+iF);
    lat.DrawLatexNDC(0.15,0.7-0.1*iF,buf);

  }//loop on files
  leg->Draw("same");
  lat.SetTextColor(1);
  lat.DrawLatexNDC(0.15,0.85,"#gamma p_{T}=20 GeV, #eta=2.0");
  lat.DrawLatexNDC(0.01,0.01,"HGCAL Geant4 Standalone Simulation");

  myc->Update();
  myc->Print("ComparisonSquareHexa_Etot_gamma_pt20_eta2.pdf");

  mycSR2->cd();
  for (unsigned iF(0); iF<nF;++iF){//loop on files
    ESR2[iF]->Scale(ESR2[1]->GetEntries()/ESR2[iF]->GetEntries());
    ESR2[iF]->Draw(iF==0?"PE":"PEsame");

    ESR2[iF]->Fit("gaus","","same");
    TF1 *fit = (TF1*)ESR2[iF]->GetFunction("gaus");
    if (!fit) continue;
    fit->SetLineColor(1+iF);
    fit->Draw("same");
    double reso = fit->GetParameter(2)/fit->GetParameter(1);
    //double err = reso*sqrt(pow(fit->GetParError(2)/fit->GetParameter(2),2)+pow(fit->GetParError(1)/fit->GetParameter(1),2));
    double err = fit->GetParError(2)/fit->GetParameter(1);
    sprintf(buf,"#sigma/E=%3.4f #pm %3.4f",reso,err);
    lat.SetTextColor(1+iF);
    lat.DrawLatexNDC(0.15,0.7-0.1*iF,buf);
  }//loop on files
  leg->Draw("same");
  lat.SetTextColor(1);
  lat.DrawLatexNDC(0.15,0.85,"#gamma p_{T}=20 GeV, #eta=2.0");
  lat.DrawLatexNDC(0.01,0.01,"HGCAL Geant4 Standalone Simulation");

  mycSR2->Update();
  mycSR2->Print("ComparisonSquareHexa_ESR2_gamma_pt20_eta2.pdf");

  myc1->Divide(7,4);
  for (unsigned iL(0);iL<nLayers;++iL){	      
    myc1->cd(iL+1);
    gPad->SetLogx(1);
    gPad->SetGridx(1);
   for (unsigned iF(0); iF<nF;++iF){//loop on files
      E1[iF][iL]->Scale(E1[1][iL]->GetEntries()/E1[iF][iL]->GetEntries());
      E1[iF][iL]->Draw(iF==0?"PE":"PEsame");
    }//loop on files
    leg->Draw("same");
    sprintf(buf,"Layer %d",iL);
    lat.DrawLatexNDC(0.15,0.85,buf);
  }
  myc1->Update();
  myc1->Print("ComparisonSquareHexa_ESR0_gamma_pt20_eta2.pdf");

  myc17->Divide(7,4);
  for (unsigned iL(0);iL<nLayers;++iL){	      
    myc17->cd(iL+1);
    for (unsigned iF(0); iF<nF;++iF){//loop on files
      E1over7[iF][iL]->Scale(E1over7[1][iL]->GetEntries()/E1over7[iF][iL]->GetEntries());
      E1over7[iF][iL]->Draw(iF==0?"PE":"PEsame");
    }//loop on files
    leg->Draw("same");
    sprintf(buf,"Layer %d",iL);
    lat.DrawLatexNDC(0.15,0.85,buf);
  }
  myc17->Update();
  myc17->Print("ComparisonSquareHexa_ESR0overSR2_gamma_pt20_eta2.pdf");

  myc7all->Divide(7,4);
  for (unsigned iL(0);iL<nLayers;++iL){	      
    myc7all->cd(iL+1);
    for (unsigned iF(0); iF<nF;++iF){//loop on files
      E7overtot[iF][iL]->Scale(E7overtot[1][iL]->GetEntries()/E7overtot[iF][iL]->GetEntries());
      E7overtot[iF][iL]->Draw(iF==0?"PE":"PEsame");
    }//loop on files
    leg->Draw("same");
    sprintf(buf,"Layer %d",iL);
    lat.DrawLatexNDC(0.15,0.85,buf);
  }
  myc7all->Update();
  myc7all->Print("ComparisonSquareHexa_ESR2overTotal_gamma_pt20_eta2.pdf");

  mycmaxout7->Divide(7,4);
  for (unsigned iL(0);iL<nLayers;++iL){	      
    mycmaxout7->cd(iL+1);
    gPad->SetLogx(1);
    gPad->SetGridx(1);
    for (unsigned iF(0); iF<nF;++iF){//loop on files
      EmaxOut7[iF][iL]->Scale(EmaxOut7[1][iL]->GetEntries()/EmaxOut7[iF][iL]->GetEntries());
      EmaxOut7[iF][iL]->Draw(iF==0?"PE":"PEsame");
    }//loop on files
    leg->Draw("same");
    sprintf(buf,"Layer %d",iL);
    lat.DrawLatexNDC(0.15,0.85,buf);
  }
  mycmaxout7->Update();
  mycmaxout7->Print("ComparisonSquareHexa_EmaxHitOutOfSR2_gamma_pt20_eta2.pdf");

  myc2D->Divide(7,4);
  for (unsigned iL(0);iL<nLayers;++iL){	      
    myc2D->cd(iL+1);
    for (unsigned iF(0); iF<nF;++iF){//loop on files
      transverseCor[iF][iL]->Scale(transverseCor[1][iL]->GetEntries()/transverseCor[iF][iL]->GetEntries());
      lname.str("");
      lname << transverseCor[iF][iL]->GetName() << "_prof";
      prof[iF][iL] = transverseCor[iF][iL]->ProfileX(lname.str().c_str());
      prof[iF][iL]->SetLineColor(1+iF);
      prof[iF][iL]->SetMarkerColor(1+iF);
      prof[iF][iL]->SetMarkerStyle(21+iF);

      transverseCor[iF][iL]->Draw(iF==0?"box":"boxsame");
    }//loop on files
    //leg->Draw("same");
    sprintf(buf,"Layer %d",iL);
    lat.DrawLatexNDC(0.15,0.85,buf);
  }

  myc2D->Update();
  myc2D->Print("ComparisonSquareHexa_E7totvsE17_gamma_pt20_eta2.pdf");

  const unsigned nPars = 3;
  double par[nF][nLayers][nPars];

  mycprof->Divide(7,4);
  for (unsigned iL(0);iL<nLayers;++iL){	      
    mycprof->cd(iL+1);
    for (unsigned iF(0); iF<nF;++iF){//loop on files
      prof[iF][iL]->Rebin(2);
      if (iL<16) prof[iF][iL]->GetYaxis()->SetRangeUser(0.5,1.);
      else prof[iF][iL]->GetYaxis()->SetRangeUser(0.1,1.);
      prof[iF][iL]->Draw(iF==0?"PE":"PEsame");
      std::ostringstream funcName;
      if (iL<8 || iL>20) funcName << "pol0";
      else funcName << "pol" << nPars-1;
      prof[iF][iL]->Fit(funcName.str().c_str(),"R","same",0.2,0.8);
      TF1 *fit = (TF1*)prof[iF][iL]->GetFunction(funcName.str().c_str());
      if (!fit) continue;
      fit->SetLineColor(7-iF);
      fit->Draw("same");
      for (int iP(0);iP<fit->GetNpar();++iP){
	par[iF][iL][iP] = fit->GetParameter(iP);
      }
    }//loop on files
    leg->Draw("same");
    sprintf(buf,"Layer %d",iL);
    lat.DrawLatexNDC(0.15,0.85,buf);
  }

  mycprof->Update();
  mycprof->Print("ComparisonSquareHexa_E7totvsE17prof_gamma_pt20_eta2.pdf");

  mycCor->cd();
  for (unsigned iF(0); iF<nF;++iF){//loop on files
    fin[iF]->cd("Energies");
    lname.str("");
    lname << "ECor" << label[iF];
    ECor[iF] = new TH1F(lname.str().c_str(),";Ecor (mips);showers",100,5000,15000);
    ECor[iF]->SetLineColor(1+iF);
    ECor[iF]->SetMarkerColor(1+iF);
    ECor[iF]->SetMarkerStyle(21+iF);
    ECor[iF]->Sumw2();

    lvar.str("");
    //lname << "(";
    for (unsigned iL(0);iL<nLayers;++iL){	      
      if (iL!=0) lvar << "+";
      lvar << "absweight_" << iL  << "*energy_" << iL << "_SR" << iSR
	   << "/(" << par[iF][iL][0];
      if (iL>7 && iL<21) lvar << "+" << par[iF][iL][1] << "*energy_" << iL << "_SR0/energy_" << iL << "_SR2";
      if (nPars>2) lvar << "+" << par[iF][iL][2] << "*pow(energy_" << iL << "_SR0/energy_" << iL << "_SR2,2)" ;
      if (nPars>3) lvar << "+" << par[iF][iL][3] << "*pow(energy_" << iL << "_SR0/energy_" << iL << "_SR2,3)";
      lvar << ")" ;
    }
    lvar << ">>" << lname.str();
    tree[iF]->Draw(lvar.str().c_str());
  }


  mycCor->cd();
  for (unsigned iF(0); iF<nF;++iF){//loop on files
    ECor[iF]->Scale(ECor[1]->GetEntries()/ECor[iF]->GetEntries());
    ECor[iF]->Draw(iF==0?"PE":"PEsame");

    ECor[iF]->Fit("gaus","","same");
    TF1 *fit = (TF1*)ECor[iF]->GetFunction("gaus");
    if (!fit) continue;
    fit->SetLineColor(1+iF);
    fit->Draw("same");
    double reso = fit->GetParameter(2)/fit->GetParameter(1);
    //double err = reso*sqrt(pow(fit->GetParError(2)/fit->GetParameter(2),2)+pow(fit->GetParError(1)/fit->GetParameter(1),2));
    double err = fit->GetParError(2)/fit->GetParameter(1);
    sprintf(buf,"#sigma/E=%3.4f #pm %3.4f",reso,err);
    lat.SetTextColor(1+iF);
    lat.DrawLatexNDC(0.15,0.7-0.1*iF,buf);
  }//loop on files
  leg->Draw("same");
  lat.SetTextColor(1);
  lat.DrawLatexNDC(0.15,0.85,"#gamma p_{T}=20 GeV, #eta=2.0");
  lat.DrawLatexNDC(0.01,0.01,"HGCAL Geant4 Standalone Simulation");

  mycCor->Update();
  mycCor->Print("ComparisonSquareHexa_ECor_gamma_pt20_eta2.pdf");




  return 0;

}//main
