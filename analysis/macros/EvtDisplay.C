#include<string>
#include<iostream>
#include<fstream>
#include<sstream>

#include "TFile.h"
#include "TTree.h"
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

double absWeight(const unsigned layer){
if (layer==0) return 1;
if (layer==1) return 1.00258;
if (layer==2) return 0.984423;
if (layer==3) return 1.00258;
if (layer==4) return 0.984423;
if (layer==5) return 1.00258;
if (layer==6) return 0.984423;
if (layer==7) return 1.00258;
if (layer==8) return 0.984423;
if (layer==9) return 1.00258;
if (layer==10) return 1.33536;
if (layer==11) return 1.3627;
if (layer==12) return 1.33536;
if (layer==13) return 1.3627;
if (layer==14) return 1.33536;
if (layer==15) return 1.3627;
if (layer==16) return 1.33536;
if (layer==17) return 1.3627;
if (layer==18) return 1.33536;
if (layer==19) return 1.3627;
if (layer==20) return 1.9495;
if (layer==21) return 1.9629;
if (layer==22) return 1.9495;
if (layer==23) return 1.9629;
if (layer==24) return 1.9495;
if (layer==25) return 1.9629;
if (layer==26) return 1.9495;
if (layer==27) return 2.01643;
if (layer==28) return 6.00121;
if (layer==29) return 5.31468;
if (layer==30) return 5.31468;
if (layer==31) return 5.31468;
if (layer==32) return 5.31468;
if (layer==33) return 5.31468;
if (layer==34) return 5.31468;
if (layer==35) return 5.31468;
if (layer==36) return 5.31468;
if (layer==37) return 5.31468;
if (layer==38) return 5.31468;
if (layer==39) return 5.31468;
if (layer==40) return 8.71728;
if (layer==41) return 8.00569;
if (layer==42) return 8.00569;
if (layer==43) return 8.00569;
if (layer==44) return 8.00569;
if (layer==45) return 8.00569;
if (layer==46) return 8.00569;
if (layer==47) return 8.00569;
if (layer==48) return 8.00569;
if (layer==49) return 8.00569;
if (layer==50) return 8.00569;
if (layer==51) return 8.00569;
 return 1;
}

int EvtDisplay(){//main  

  //const unsigned nS = 7;
  //TString scenario[nS] = {"0","1","2","3","4","5","6"};
  TString scenario = "VBFHgg/";
  //TString scenario = "ZHtautau/concept/";
  //TString scenario = "ZH125/concept/";
  //TString scenario = "GluGlu/concept/";
  bool doAll = true;
  unsigned versionNumber = 33;

  TString version = "";
  version += versionNumber;

  unsigned mipThresh = 1;


  double minRadius = 316;//mm

  ////double minX=-150,maxX=150;
  //double minX=-700,maxX=700;
  ////double minY=420,maxY=720;
  ////double minY=1150,maxY=1450;
  //double minY=-1100,maxY=500;

  std::ostringstream ltitleBase;
  //ltitle << "#gamma 200 GeV"
  ltitleBase << "VBF H,H#rightarrow#gamma#gamma";
  //ltitleBase << "ZH,H#rightarrow#tau#tau";
  //ltitleBase << "pp #rightarrow gg";
    //<< " Event #" << event[ievt]
    //	     << ", E_{1#times1 cm^{2}} > "
    //	     << mipThresh << " Mips";


  const unsigned nLayers = 52;
  const unsigned nEcalLayers = 28;

  //VBFH
  const unsigned nEvts = 7;//1000;
  unsigned event[nEvts]={4,5,6,7,9,11,12};//,6,12};//103,659,875};
  //Htautau 1000
  //const unsigned nEvts = 7;//1000;
  //unsigned event[nEvts]={1,2,5,8,10,11,12};//,6,12};//103,659,875};
  //Htautau 125
  //const unsigned nEvts = 5;//1000;
  //unsigned event[nEvts]={0,1,2,10,14};//,6,12};//103,659,875};
  //gluons
  //const unsigned nEvts = 3;//1000;
  //unsigned event[nEvts]={7,16,22};//,2,3,4};//,6,12};//103,659,875};

  //  const unsigned nEvts = 13;//1000;
  //unsigned event[nEvts];//={1,6,12};//103,659,875};
  //for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on events
  //event[ievt] = ievt;
  //}

  double minZ=3170,maxZ=5000;


  //VBFH
  double minX[nEvts] = {-700,-600,-400,-400,0,0,0};
  double maxX[nEvts] = {-100,0,200,100,700,700,500};
  double minY[nEvts] = {100,-600,-1100,-800,-100,-200,-700};
  double maxY[nEvts] = {480,0,0,0,400,300,0};

  //Htautau 1000
  /*double minX[nEvts] = {0,-500,0,-400,0,-500,-700};
  double maxX[nEvts] = {600,0,700,500,700,0,600};
  double minY[nEvts] = {0,-700,-100,-700,-400,0,-400};
  double maxY[nEvts] = {500,0,400,900,200,700,700};*/

  //Htautau 125
  /*double minX[nEvts] = {-300,-850,-200,-450,-550};
  double maxX[nEvts] = {100,-450,300,0,-150};
  double minY[nEvts] = {-900,-280,200,-950,-800};
  double maxY[nEvts] = {-500,150,650,-550,-400};*/

  //pp to gg
  //double minX[nEvts] = {100,100,-200};
  //double maxX[nEvts] = {700,700,500};
  //double minY[nEvts] = {100,100,200};
  //double maxY[nEvts] = {700,700,900};

  if (doAll){
    for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on events
      minX[ievt]=-1700;maxX[ievt]=1700;
      minY[ievt]=-1700;maxY[ievt]=1700;
    }
  }

  //if (nLayers == nEcalLayers) {
  //maxZ = 3500;
  //}
  unsigned nZ=(maxZ-minZ)/10.;


  //if (doAll){
  //minX=-1700,maxX=1700;
  //minY=-1700,maxY=1700;
  //}

  std::ostringstream lName;

  const unsigned nCanvas = nEvts;  
  TCanvas *myc[nCanvas];
  //TPad *pad1 = 0;

  for (unsigned iC(0);iC<nCanvas;++iC){
    lName.str("");
    lName << "myc" << iC;
    myc[iC] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1400,1000);
    //myc[iC]->Divide(2,1);
  }

  std::ostringstream saveName;
  
  std::ostringstream inputDir;
  inputDir << "/afs/cern.ch/work/a/amagnan/public/HGCalGeant4/git_V00-03-06/version_33/model_2/";
  inputDir << scenario ;
  inputDir << "/BON/DigiPFcal.root" ;

  TString outputDir = "../PLOTS/version_"+version+"/"+scenario+"/";
  if (doAll) outputDir += "Overview/";

  TFile *inputFile = TFile::Open(inputDir.str().c_str());
  if (!inputFile) {
    std::cout << " -- Error, input file " << inputDir.str() << " cannot be opened." << std::endl;
    return 1;
  }
  inputFile->cd();
  TTree *tree = (TTree*)gDirectory->Get("RecoTree");
  TH3F *p_xyz[nEvts];

  //for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on events
  for (unsigned ievt(0); ievt<2; ++ievt){//loop on events


    TString plotDir = outputDir;
    //if (doAll) plotDir += "Overview/";

    unsigned nX=(maxX[ievt]-minX[ievt])/10.,nY=(maxY[ievt]-minY[ievt])/10.;
    lName.str("");
    lName << "xyz_" << ievt ; 

    p_xyz[ievt] = new TH3F(lName.str().c_str(),";z(mm);x(mm);y(mm)",
			   nZ,minZ,maxZ,
			   nX,minX[ievt],maxX[ievt],
			   nY,minY[ievt],maxY[ievt]);
      
    p_xyz[ievt]->Sumw2();
    //p_xyz->RebinY(2);
    //p_xyz->RebinZ(2);


    set_plot_style();
    gStyle->SetOptStat(1111111);

    myc[ievt]->cd();
    std::ostringstream lvar;
    lvar << "HGCSSRecoHitVec.ypos_:HGCSSRecoHitVec.xpos_:HGCSSRecoHitVec.zpos_>>" << lName.str();
    std::ostringstream lcut;
    lcut << "event_==" << event[ievt];
    tree->Draw(lvar.str().c_str(),lcut.str().c_str());

    std::cout << " --bins, min and max of hist = " 
	      << p_xyz[ievt]->GetNbinsY() << " " 
	      << p_xyz[ievt]->GetMinimum() << " " 
	      << p_xyz[ievt]->GetMaximum() 
	      << std::endl;

    //continue;

    std::ostringstream ltitle;
    ltitle << ltitleBase.str();


    const unsigned nSlices = 6;
    int lColor[nSlices] = {1,4,7,3,5,2};
    TH3F *p_slice[nSlices];
    double val[nSlices+1] = {1,5,20,100,300,500,10000};
    for (unsigned iS(0); iS<nSlices; ++iS){
      lName.str("");
      lName << "slice_" << iS ; 
      p_slice[iS] = new TH3F(lName.str().c_str(),";z(mm);x(mm);y(mm)",
			     p_xyz[ievt]->GetNbinsX(),p_xyz[ievt]->GetXaxis()->GetBinLowEdge(1),p_xyz[ievt]->GetXaxis()->GetBinLowEdge(p_xyz[ievt]->GetNbinsX()+1),
			     p_xyz[ievt]->GetNbinsY(),p_xyz[ievt]->GetYaxis()->GetBinLowEdge(1),p_xyz[ievt]->GetYaxis()->GetBinLowEdge(p_xyz[ievt]->GetNbinsY()+1),
			     p_xyz[ievt]->GetNbinsZ(),p_xyz[ievt]->GetZaxis()->GetBinLowEdge(1),p_xyz[ievt]->GetZaxis()->GetBinLowEdge(p_xyz[ievt]->GetNbinsZ()+1)
			     );
      p_slice[iS]->SetMarkerColor(lColor[iS]);
      //p_slice[iS]->SetLineColor(iS+1);
      //p_slice[iS]->SetFillColor(iS+1);
      //p_slice[iS]->SetMarkerStyle(7);
      p_slice[iS]->GetXaxis()->SetRangeUser(minZ,maxZ);
      p_slice[iS]->GetYaxis()->SetRangeUser(minX[ievt],maxX[ievt]);
      p_slice[iS]->GetZaxis()->SetRangeUser(minY[ievt],maxY[ievt]);
    
      p_slice[iS]->GetXaxis()->SetLabelSize(0.05);
      p_slice[iS]->GetYaxis()->SetLabelSize(0.05);
      p_slice[iS]->GetZaxis()->SetLabelSize(0.05);
      p_slice[iS]->GetXaxis()->SetTitleSize(0.05);
      p_slice[iS]->GetYaxis()->SetTitleSize(0.05);
      p_slice[iS]->GetZaxis()->SetTitleSize(0.05);
      p_slice[iS]->SetTitle(ltitle.str().c_str());

      //val[iS]= exp(iS*log(p_xyz[ievt]->GetMaximum())/nSlices);
    }

    //val[nSlices] = p_xyz[ievt]->GetMaximum();

    double lmin = 10;
    double lmax = p_xyz[ievt]->GetBinContent(1,1,1);

    for (int xb(1); xb<p_xyz[ievt]->GetNbinsX()+1;++xb){
      //std::cout << "--- xb=" << xb << "" << std::endl;
      for (int yb(1); yb<p_xyz[ievt]->GetNbinsY()+1;++yb){
	double xTmp = p_xyz[ievt]->GetYaxis()->GetBinCenter(yb);
	for (int zb(1); zb<p_xyz[ievt]->GetNbinsZ()+1;++zb){
	  double yTmp = p_xyz[ievt]->GetZaxis()->GetBinCenter(zb);
	  double radius = sqrt(yTmp*yTmp+xTmp*xTmp);
	  if (radius < minRadius){
	    p_xyz[ievt]->SetBinContent(xb,yb,zb,0);
	    continue;
	  }
	  double ltmp = p_xyz[ievt]->GetBinContent(xb,yb,zb);
	  if (ltmp<lmin && ltmp>0) lmin=ltmp;
	  if (ltmp>lmax) lmax=ltmp;

	  if (ltmp<1) continue;
	  //std::cout << xb << " " << yb << " " << zb << " " << p_xyz[ievt]->GetBinContent(xb,yb,zb) << std::endl;
	  if (ltmp < mipThresh) {
		  p_xyz[ievt]->SetBinContent(xb,yb,zb,0);
	  }
	  else {
	    p_xyz[ievt]->SetBinContent(xb,yb,zb,ltmp);
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
    //TPad *pad1 = new TPad("pad1","This is pad1",0.01,0.01,0.83,0.99);
    //pad1->Draw();
    //pad1->cd();    
    //myc[ievt]->cd(1);
    for (unsigned iS(0); iS<nSlices; ++iS){
      std::cout << " -- slice " << iS << std::endl;
      //set artificially the boundaries: min and max in dummy bins
      //p_slice[iS]->SetBinContent(1,1,1,1);
      //p_slice[iS]->SetBinContent(2,1,1,lmax);
      p_slice[iS]->Draw(iS==0? "" : "same");
    }
    //myc[ievt]->cd(2);
    //TPad *pad2 = new TPad("pad2","This is pad2",0.85,0.01,0.99,0.99);
    //pad2->Draw();
    //pad2->cd();
    TLatex lat;
    lat.SetTextAlign(31);
    //lat.SetTextSize(0.2);
    for (unsigned iS(0); iS<nSlices; iS++){
      //std::cout << iS << std::endl;
      char buf[500];
      sprintf(buf,"%3.1f < Mips < %3.1f",val[iS],val[iS+1]);
      lat.SetTextColor(lColor[iS]);
      lat.DrawLatexNDC(1.,(iS*1.0/nSlices)/2.+0.5,buf);
    }

    myc[ievt]->Update();
    saveName.str("");
    saveName << plotDir << "/evtDisplay_evt" << event[ievt] << "_mipThresh" << mipThresh;
    myc[ievt]->Print((saveName.str()+".png").c_str());
    myc[ievt]->Print((saveName.str()+".pdf").c_str());
    
  }//loop on events

  return 0;
  
  
}//main

