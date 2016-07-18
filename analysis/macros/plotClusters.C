#include "TH1F.h"
#include <string>
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TDRStyle.h"

int plotClusters(){//main

  SetTdrStyle();

  const unsigned nPu = 2;
  TFile *f[nPu];
  f[0] = TFile::Open("200um/eta17_et50_pu0.root");
  f[1] = TFile::Open("200um/eta17_et50_pu140.root");

  TCanvas *myc0 = new TCanvas("myc0","myc0",1500,1000);
  myc0->Divide(3,2);
  TCanvas *myc1 = new TCanvas("myc1","myc1",1500,1000);
  myc1->Divide(3,2);
  TCanvas *myc2 = new TCanvas("myc2","myc2",1500,1000);
  myc2->Divide(3,2);
 
  gStyle->SetOptStat(0);

  TH1F *p_nClusters[2];

  TH1F *p_clusnHits_all[2];
  TH1F *p_seedEoverE_all[2];
  TH1F *p_clusLayer_all[2];
  TH1F *p_clusWidth_all[2];
  TH1F *p_seeddeta_all[2];
  TH1F *p_seeddphi_all[2];

  TH1F *p_clusnHits_sel[2];
  TH1F *p_seedEoverE_sel[2];
  TH1F *p_clusLayer_sel[2];
  TH1F *p_clusWidth_sel[2];
  TH1F *p_seeddeta_sel[2];
  TH1F *p_seeddphi_sel[2];

  TString puStr[nPu] = {"pu0","pu140"};

  for (unsigned ipu(0); ipu<nPu;++ipu){//loop on pu
    f[ipu]->cd("PositionFit");
    p_nClusters[ipu] = (TH1F*)((TH1F*)(gDirectory->Get("p_nClusters"))->Clone("p_nClusters"+puStr[ipu]));
    
    p_clusnHits_all[ipu] = (TH1F*)gDirectory->Get("p_clusnHits_all");
    p_clusnHits_all[ipu]->SetTitle(";n_{hits};n_{clusters}");
    //if (ipu==1) 
    p_clusnHits_all[ipu]->SetDirectory(0);
    p_seedEoverE_all[ipu] = (TH1F*)gDirectory->Get("p_seedEoverE_all")->Clone("p_seedEoverE_all"+puStr[ipu]);
    p_clusLayer_all[ipu] = (TH1F*)gDirectory->Get("p_clusLayer_all")->Clone("p_clusLayer_all"+puStr[ipu]);
    p_clusWidth_all[ipu] = (TH1F*)gDirectory->Get("p_clusWidth_all")->Clone("p_clusWidth_all"+puStr[ipu]);
    p_seeddeta_all[ipu] = (TH1F*)gDirectory->Get("p_seeddeta_all")->Clone("p_seeddeta_all"+puStr[ipu]);
    p_seeddphi_all[ipu] = (TH1F*)gDirectory->Get("p_seeddphi_all")->Clone("p_seeddphi_all"+puStr[ipu]);
    
    p_clusnHits_sel[ipu] = (TH1F*)gDirectory->Get("p_clusnHits_sel")->Clone("p_clusnHits_sel"+puStr[ipu]);
    p_clusnHits_sel[ipu]->SetTitle(";n_{hits};n_{events}");
    p_seedEoverE_sel[ipu] = (TH1F*)gDirectory->Get("p_seedEoverE_sel")->Clone("p_seedEoverE_sel"+puStr[ipu]);
    p_clusLayer_sel[ipu] = (TH1F*)gDirectory->Get("p_clusLayer_sel")->Clone("p_clusLayer_sel"+puStr[ipu]);
    p_clusWidth_sel[ipu] = (TH1F*)gDirectory->Get("p_clusWidth_sel")->Clone("p_clusWidth_sel"+puStr[ipu]);
    p_seeddeta_sel[ipu] = (TH1F*)gDirectory->Get("p_seeddeta_sel")->Clone("p_seeddeta_sel"+puStr[ipu]);
    p_seeddphi_sel[ipu] = (TH1F*)gDirectory->Get("p_seeddphi_sel")->Clone("p_seeddphi_sel"+puStr[ipu]);
  }    
    
  //all vs sel for 0 PU

  TLegend *leg1 = new TLegend(0.7,0.7,0.94,0.94);
  leg1->SetFillColor(10);

  myc0->cd(1);
  gPad->SetLogy(1);
  p_clusnHits_all[0]->SetMinimum(0.1);
  p_clusnHits_all[0]->Draw();
  p_clusnHits_sel[0]->SetLineColor(2);
  p_clusnHits_sel[0]->Draw("same");
  leg1->AddEntry(p_clusnHits_all[0],"All","L");
  leg1->AddEntry(p_clusnHits_sel[0],"Sel","L");
  leg1->Draw("same");

  myc0->cd(2);
  gPad->SetLogy(1);
  p_seedEoverE_all[0]->SetMinimum(0.1);
  p_seedEoverE_all[0]->Draw();
  p_seedEoverE_sel[0]->SetLineColor(2);
  p_seedEoverE_sel[0]->Draw("same");
  leg1->Draw("same");

  myc0->cd(3);
  gPad->SetLogy(1);
  p_clusLayer_all[0]->SetMinimum(0.1);
  p_clusLayer_all[0]->Draw();
  p_clusLayer_sel[0]->SetLineColor(2);
  p_clusLayer_sel[0]->Draw("same");
  leg1->Draw("same");

  myc0->cd(4);
  gPad->SetLogy(1);
  p_clusWidth_all[0]->SetMinimum(0.1);
  p_clusWidth_all[0]->Draw();
  p_clusWidth_sel[0]->SetLineColor(2);
  p_clusWidth_sel[0]->Draw("same");
  leg1->Draw("same");

  myc0->cd(5);
  gPad->SetLogy(1);
  p_seeddeta_all[0]->SetMinimum(0.1);
  p_seeddeta_all[0]->Draw();
  p_seeddeta_sel[0]->SetLineColor(2);
  p_seeddeta_sel[0]->Draw("same");
  leg1->Draw("same");


  myc0->cd(6);
  gPad->SetLogy(1);
  p_seeddphi_all[0]->SetMinimum(0.1);
  p_seeddphi_all[0]->Draw();
  p_seeddphi_sel[0]->SetLineColor(2);
  p_seeddphi_sel[0]->Draw("same");
  leg1->Draw("same");


  myc0->Update();
  myc0->Print("clusterSigvsBkg_0pu.pdf");

  //all vs sel for 140 PU
 
  myc1->cd(1);
  gPad->SetLogy(1);
  p_clusnHits_all[1]->SetMinimum(0.1);
  p_clusnHits_all[1]->Draw();
  p_clusnHits_sel[1]->SetLineColor(2);
  p_clusnHits_sel[1]->Draw("same");
  leg1->Draw("same");

  myc1->cd(2);
  gPad->SetLogy(1);
  p_seedEoverE_all[1]->SetMinimum(0.1);
  p_seedEoverE_all[1]->Draw();
  p_seedEoverE_sel[1]->SetLineColor(2);
  p_seedEoverE_sel[1]->Draw("same");
  leg1->Draw("same");

  myc1->cd(3);
  gPad->SetLogy(1);
  p_clusLayer_all[1]->SetMinimum(0.1);
  p_clusLayer_all[1]->Draw();
  p_clusLayer_sel[1]->SetLineColor(2);
  p_clusLayer_sel[1]->Draw("same");
  leg1->Draw("same");

  myc1->cd(4);
  gPad->SetLogy(1);
  p_clusWidth_all[1]->SetMinimum(0.1);
  p_clusWidth_all[1]->Draw();
  p_clusWidth_sel[1]->SetLineColor(2);
  p_clusWidth_sel[1]->Draw("same");
  leg1->Draw("same");

  myc1->cd(5);
  gPad->SetLogy(1);
  p_seeddeta_all[1]->SetMinimum(0.1);
  p_seeddeta_all[1]->Draw();
  p_seeddeta_sel[1]->SetLineColor(2);
  p_seeddeta_sel[1]->Draw("same");
  leg1->Draw("same");


  myc1->cd(6);
  gPad->SetLogy(1);
  p_seeddphi_all[1]->SetMinimum(0.1);
  p_seeddphi_all[1]->Draw();
  p_seeddphi_sel[1]->SetLineColor(2);
  p_seeddphi_sel[1]->Draw("same");
  leg1->Draw("same");


  myc1->Update();
  myc1->Print("clusterSigvsBkg_140pu.pdf");

  TLegend *leg2 = new TLegend(0.7,0.7,0.94,0.94);
  leg2->SetFillColor(10);

  //sel 0 vs sel 140pu
  for (unsigned ipu(0); ipu<nPu;++ipu){//loop on pu
    myc2->cd(1);
    p_clusnHits_sel[ipu]->SetLineColor(1+ipu);
    p_clusnHits_sel[ipu]->Draw(ipu==0?"":"same");
    leg2->AddEntry(p_clusnHits_sel[ipu],puStr[ipu],"L");
    if (ipu==1) leg2->Draw("same");

    myc2->cd(2);
    p_seedEoverE_sel[ipu]->SetLineColor(1+ipu);
    p_seedEoverE_sel[ipu]->Draw(ipu==0?"":"same");
    if (ipu==1) leg2->Draw("same");

    myc2->cd(3);
    p_clusLayer_sel[ipu]->SetLineColor(1+ipu);
    p_clusLayer_sel[ipu]->Draw(ipu==0?"":"same");
    if (ipu==1) leg2->Draw("same");

    myc2->cd(4);
    p_clusWidth_sel[ipu]->SetLineColor(1+ipu);
    p_clusWidth_sel[ipu]->Draw(ipu==0?"":"same");
    if (ipu==1) leg2->Draw("same");

    myc2->cd(5);
    p_seeddeta_sel[ipu]->SetLineColor(1+ipu);
    p_seeddeta_sel[ipu]->Draw(ipu==0?"":"same");
    if (ipu==1) leg2->Draw("same");

    myc2->cd(6);
    p_seeddphi_sel[ipu]->SetLineColor(1+ipu);
    p_seeddphi_sel[ipu]->Draw(ipu==0?"":"same");
    if (ipu==1) leg2->Draw("same");

  }//loop on pu

  myc2->Update();
  myc2->Print("cluster0puvs140pu.pdf");







}//main
