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
#include "TColor.h"
#include "TEveManager.h"
#include "TGLViewer.h"
#include "TF2.h"
#include "TEvePlot3D.h"
#include "TEveTrans.h"

// void glplot(TH3F *h31)
// {
//   TEveManager::Create();
//   gEve->GetDefaultGLViewer()->SetCurrentCamera(TGLViewer::kCameraPerspXOY);
  
//   h31->SetFillColor(2);
//   TEvePlot3D* x = new TEvePlot3D("EvePlot - TH3F");
//   x->SetPlot(h31, "glbox");
//   x->RefMainTrans().MoveLF(2, -1);
//   gEve->AddElement(x);

//   gEve->Redraw3D(kTRUE);
// }

void set_plot_style()
{
    const Int_t NRGBs = 7;
    const Int_t NCont = 255;

    //Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    //Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    //Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    //Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    //Double_t stops[NRGBs] = { 0.00, 0.20, 0.30, 0.45, 0.61, 0.84, 1.00 };
    Double_t stops[NRGBs] = { 0.00, 0.10, 0.20, 0.35, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 1.00, 1.00, 0.00, 0.00, 0.50, 1.00, 0.50 };
    Double_t green[NRGBs] = { 1.00, 1.00, 1.00, 0.81, 0.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.50, 0.00, 0.10, 1.00, 0.50, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}

int main() {//EvtDisplay(){//main  

  //const unsigned nS = 7;
  //TString scenario[nS] = {"0","1","2","3","4","5","6"};
  TString scenario = "VBFH/concept/";
  bool doAll = true;
  TString version = "20";
  unsigned mipThresh = 1;

  //double minX=-150,maxX=150;
  double minX=-400,maxX=50;
  //double minY=420,maxY=720;
  //double minY=1150,maxY=1450;
  double minY=320,maxY=900;

  std::ostringstream ltitleBase;
  //ltitle << "#gamma 200 GeV"
  ltitleBase << "VBF jet"
    //<< " Event #" << event[ievt]
	     << ", E_{1#times1 cm^{2}} > "
	     << mipThresh << " Mips";


  const unsigned nLayers = 64;
  const unsigned nEcalLayers = 31;

  const unsigned nEvts = 10;//1000;
  unsigned event[nEvts];//={1,6,12};//103,659,875};
  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on events
    event[ievt] = ievt;
  }

  double minZ=3170,maxZ=5000;

  if (doAll){
    minX=-1700,maxX=1700;
    minY=-1700,maxY=1700;
  }

  std::ostringstream lName;

  const unsigned nCanvas = nEvts;  
  TCanvas *myc[nCanvas];
  for (unsigned iC(0);iC<nCanvas;++iC){
    lName.str("");
    lName << "myc" << iC;
    myc[iC] = new TCanvas(lName.str().c_str(),lName.str().c_str(),800,800);
  }

  std::ostringstream saveName;
  
  TString inputDir = "PLOTS/version_"+version+"/"+scenario+"/";

  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on events
    lName.str("");
    lName << inputDir << "CalibHistos_E200_evt" << event[ievt] << ".root";
    TFile *inputFile = TFile::Open(lName.str().c_str());
    if (!inputFile) {
      std::cout << " -- Error, input file " << lName.str() << " cannot be opened. Going to next..." << std::endl;
      continue;
      //return 1;
    }


    TString plotDir = inputDir;
    if (doAll) plotDir += "Overview/";

    inputFile->cd();
    TH3F *p_xyz = 0;
    p_xyz = (TH3F*)gDirectory->Get("p_xyz")->Clone();
      
    if (!p_xyz) {
      std::cout << " -- ERROR, pointer for XYZ histogram is null. Exiting..." << std::endl;
      return 1;
    }
    p_xyz->Sumw2();

    std::cout << " --min and max of hist = " 
	      <<  p_xyz->GetMinimum() << " " 
	      << p_xyz->GetMaximum() 
	      << std::endl;

    set_plot_style();

    std::ostringstream ltitle;
    ltitle << ltitleBase.str();


    const unsigned nSlices = 6;
    int lColor[nSlices] = {1,4,7,3,5,2};
    TH3F *p_slice[nSlices];
    double val[nSlices+1];
    for (unsigned iS(0); iS<nSlices; ++iS){
      lName.str("");
      lName << "slice_" << iS ; 
      p_slice[iS] = new TH3F(lName.str().c_str(),";z(mm);x(mm);y(mm)",
			     p_xyz->GetNbinsX(),p_xyz->GetXaxis()->GetBinLowEdge(1),p_xyz->GetXaxis()->GetBinLowEdge(p_xyz->GetNbinsX()+1),
			     p_xyz->GetNbinsY(),p_xyz->GetYaxis()->GetBinLowEdge(1),p_xyz->GetYaxis()->GetBinLowEdge(p_xyz->GetNbinsY()+1),
			     p_xyz->GetNbinsZ(),p_xyz->GetZaxis()->GetBinLowEdge(1),p_xyz->GetZaxis()->GetBinLowEdge(p_xyz->GetNbinsZ()+1)
			     );
      p_slice[iS]->SetMarkerColor(lColor[iS]);
      //p_slice[iS]->SetLineColor(iS+1);
      //p_slice[iS]->SetFillColor(iS+1);
      //p_slice[iS]->SetMarkerStyle(7);
      p_slice[iS]->GetXaxis()->SetRangeUser(minZ,maxZ);
      p_slice[iS]->GetYaxis()->SetRangeUser(minX,maxX);
      p_slice[iS]->GetZaxis()->SetRangeUser(minY,maxY);
    
      p_slice[iS]->GetXaxis()->SetLabelSize(0.05);
      p_slice[iS]->GetYaxis()->SetLabelSize(0.05);
      p_slice[iS]->GetZaxis()->SetLabelSize(0.05);
      p_slice[iS]->GetXaxis()->SetTitleSize(0.05);
      p_slice[iS]->GetYaxis()->SetTitleSize(0.05);
      p_slice[iS]->GetZaxis()->SetTitleSize(0.05);
      p_slice[iS]->SetTitle(ltitle.str().c_str());

      val[iS]= exp(iS*log(p_xyz->GetMaximum())/nSlices);
    }
    val[nSlices] = p_xyz->GetMaximum();

    double lmin = 10;
    double lmax = p_xyz->GetBinContent(1,1,1);

    for (int xb(1); xb<p_xyz->GetNbinsX()+1;++xb){
      std::cout << "--- xb=" << xb << "" << std::endl;
      for (int yb(1); yb<p_xyz->GetNbinsY()+1;++yb){
	for (int zb(1); zb<p_xyz->GetNbinsZ()+1;++zb){
	  double ltmp = p_xyz->GetBinContent(xb,yb,zb);
	  if (ltmp<lmin && ltmp>0) lmin=ltmp;
	  if (ltmp>lmax) lmax=ltmp;

	  if (ltmp<1) continue;
	  //std::cout << xb << " " << yb << " " << zb << " " << p_xyz->GetBinContent(xb,yb,zb) << std::endl;
	  if (ltmp < mipThresh) {
		  p_xyz->SetBinContent(xb,yb,zb,0);
	  }
	  else {
	    p_xyz->SetBinContent(xb,yb,zb,ltmp);
	    for (unsigned iS(0); iS<nSlices; ++iS){
	      if (ltmp<=val[iS+1] && ltmp>val[iS])
		p_slice[iS]->SetBinContent(xb,yb,zb,ltmp);
	    }
	  }
	}
      }
    }

    std::cout << " --min and max by hand = " 
	      << lmin << " " 
	      << lmax 
	      << std::endl;

    myc[ievt]->cd();
    TLatex lat;
    for (unsigned iS(0); iS<nSlices; ++iS){
      std::cout << " -- slice " << iS << std::endl;
      //set artificially the boundaries: min and max in dummy bins
      p_slice[iS]->SetBinContent(1,1,1,1);
      p_slice[iS]->SetBinContent(2,1,1,lmax);

      p_slice[iS]->Draw(iS==0? "" : "same");
      char buf[500];
      sprintf(buf,"%3.1f < Mips < %3.1f",val[iS],val[iS+1]);
      lat.SetTextColor(lColor[iS]);
      lat.DrawLatexNDC(0.12,(nSlices-iS)/nSlices,buf);
    }


    myc[ievt]->Update();
    saveName.str("");
    saveName << plotDir << "/evtDisplay_evt" << event[ievt] << "_mipThresh" << mipThresh;
    myc[ievt]->Print((saveName.str()+".png").c_str());
    myc[ievt]->Print((saveName.str()+".pdf").c_str());
    
  }//loop on events

  return 0;
  
  
}//main

