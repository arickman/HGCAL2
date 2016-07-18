#include "TFile.h"
#include "TString.h"
#include "TNtuple.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TMath.h"
#include "TSpline.h"

#include "Math/WrappedTF1.h"
#include "Math/BrentMinimizer1D.h"

#include <map>

#include "HGCSSCaloProperties.hh"

//
int main()
{
  setStyle();

  enum DetectorVersion { v_CALICE=0,
                         v_HGCALEE_Si80=1,
                         v_HGCALEE_Si120=2,
                         v_HGCALEE_Si200=3,
                         v_HGCALEE_Si500=4,
                         v_HGCALEE_gap1=5,
                         v_HGCALEE_CALICE=6,
                         v_HGCALEE_inverted=7,
                         v_HGCALEE_concept=8,
			 v_HGCALEE_W=9,
			 v_HGCALEE_gap4=10 ,
			 v_HGCALEE_fineSampling=11,
			 v_HGCALEE_samplingOptim=12
  };

  /*
  size_t versions[]={
    v_CALICE,
    v_HGCALEE_CALICE,
    v_HGCALEE_Si80,
    v_HGCALEE_Si120,
    v_HGCALEE_Si200,
    v_HGCALEE_Si500,
    v_HGCALEE_W,
    //    v_HGCALEE_fineSampling,
    v_HGCALEE_inverted,
    v_HGCALEE_concept,
    v_HGCALEE_gap1,
    v_HGCALEE_gap4,
    v_HGCALEE_gap4
  };
  */
  size_t versions[]={v_HGCALEE_concept};
  //size_t versions[]={v_HGCALEE_fineSampling}; 
  //  size_t versions[]={v_HGCALEE_Si200};
  //size_t versions[]={v_CALICE,v_HGCALEE_Si200, v_HGCALEE_Si500,v_HGCALEE_CALICE, v_HGCALEE_W};
  //size_t versions[]={v_HGCALEE_Si80,v_HGCALEE_Si120,v_HGCALEE_Si200,v_HGCALEE_Si500};
  //size_t versions[]={v_HGCALEE_Si200,v_HGCALEE_gap1,v_HGCALEE_concept,v_HGCALEE_inverted};
  size_t nversions( sizeof(versions)/sizeof(size_t) );

  std::vector<TH1F *> stochGr, constGr;
  for(size_t i=0; i<4; i++)
    {
      TString pf(""); pf+=i;
      TH1F *h=new TH1F("stochgr"+pf,";Detector version;(#sigma/E)_{stochastic}",nversions,0,nversions); h->SetMarkerStyle(20+i); h->SetDirectory(0); stochGr.push_back(h);
      h=new TH1F("constgr"+pf,";Detector version;(#sigma/E)_{cte}",nversions,0,nversions);              h->SetMarkerStyle(20+i); h->SetDirectory(0); constGr.push_back(h);
    }  

  for(size_t iv=0; iv<nversions; iv++){
    size_t i=versions[iv];
   
    TString binLabel("CALICE");
    if(i==v_HGCALEE_Si80)     binLabel="Si 80 #mum";
    if(i==v_HGCALEE_Si120)    binLabel="Si 120 #mum";
    if(i==v_HGCALEE_Si200)    binLabel="Si 200 #mum";
    if(i==v_HGCALEE_Si500)    binLabel="Si 500 #mum";
    if(i==v_HGCALEE_gap1)     binLabel="Air gap 1mm";
    if(i==v_HGCALEE_CALICE)   binLabel="#splitline{CALICE-like}{sampling}";
    if(i==v_HGCALEE_fineSampling)   binLabel="Fine sampling";
    if(i==v_HGCALEE_inverted) binLabel="Inverted";
    if(i==v_HGCALEE_concept)  binLabel="Concept";
    if(i==v_HGCALEE_W)        binLabel="W absorber";
    if(i==v_HGCALEE_gap4)     binLabel="Air gap 4mm";

    //TString ver("version_"); ver+=i;
    TString ver("version"); ver+=i;

    CaloProperties props(ver);
    props.characterizeCalo();

    for(size_t ialgo=0; ialgo<3; ialgo++)
      {
	stochGr[ialgo]->SetTitle(props.resCurve_[ialgo]->GetTitle());
	stochGr[ialgo]->GetXaxis()->SetBinLabel(iv+1,binLabel);
	stochGr[ialgo]->SetBinContent          (iv+1,props.stochTerms_[ialgo].first);
	stochGr[ialgo]->SetBinError            (iv+1,props.stochTerms_[ialgo].second);
	constGr[ialgo]->SetTitle(props.resCurve_[ialgo]->GetTitle());
	constGr[ialgo]->GetXaxis()->SetBinLabel(iv+1,binLabel);
	constGr[ialgo]->SetBinContent          (iv+1,props.constTerms_[ialgo].first);
	constGr[ialgo]->SetBinError            (iv+1,props.constTerms_[ialgo].second);
      }
  }

  TCanvas *csum=new TCanvas("csum","csum",1200,600);
  csum->Divide(2,1);
  csum->cd(1);
  TLegend *leg=new TLegend(0.15,0.85,0.9,0.95);
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->SetNColumns(3);
  for(Int_t i=0; i<2; i++)
    {
      TPad *p=(TPad *)csum->cd(i+1);
      p->SetTopMargin(0.05);
      p->SetBottomMargin(0.15);

      for(size_t ialgo=0; ialgo<3; ialgo++)
	{	  
	  TH1F *gr=(i==0 ? constGr[ialgo] : stochGr[ialgo]);
	  gr->Draw(ialgo==0 ? "e1" : "e1same");
	  gr->GetYaxis()->SetTitleOffset(1.4);
	  gr->GetYaxis()->SetTitleSize(0.05);
	  gr->GetYaxis()->SetLabelSize(0.04);
	  gr->GetYaxis()->SetRangeUser(0,i==0 ? 0.05 : 0.5);
	  gr->GetYaxis()->SetNdivisions(10);
	  gr->GetXaxis()->SetTitleSize(0.05);
	  gr->GetXaxis()->SetLabelSize(0.04);
	  gr->GetXaxis()->SetTitleOffset(1.1);
	  if(i==0) leg->AddEntry(gr,gr->GetTitle(),"p");
	  if(i==0 && ialgo==0) drawHeader();
	}
    }
  csum->cd(1);
  leg->Draw();
  csum->cd();
  csum->Modified();
  csum->Update();
  csum->SaveAs("PLOTS/CaloPerformanceSummary.png");
  csum->SaveAs("PLOTS/CaloPerformanceSummary.root");

  return 0;

}//main
