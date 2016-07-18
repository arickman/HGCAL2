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

double minR(const unsigned layer){
  if (layer<28) return 750;
  else if (layer<40) return 600;
  else return 0;
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

void set_plot_style()
{
    const Int_t NRGBs = 7;
    const Int_t NCont = 64;

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

int plotXYZ(){//main  

  //const unsigned nS = 7;
  //TString scenario[nS] = {"0","1","2","3","4","5","6"};
  const unsigned nS = 1;
  std::string scenario[nS] = {
    //"quark_u/eta30/",
    //"quark_u/PU/eta30/"
    //"gamma_eta25/GeVCal/",
    //"gamma_eta25/PU/GeVCal/"
    "VBFHgg/",
    //"VBFH/PU/concept/"
    //"ZH125/concept/",
    //"ZH125/PU/concept/"
    //"GluGlu/concept/",
    //"GluGlu/PU/concept/"
  };

  bool doVBF = true;
  bool doGlu = false;

  std::ostringstream ltitleBase;
  //ltitle << "#gamma 200 GeV"
  if (doVBF) ltitleBase << "VBF jet";
  else if (doGlu) ltitleBase << "Gluon jet";
  else ltitleBase << "ZH,H#rightarrow#tau#tau";
  //<< " Event #" << event[ievt]
  //ltitleBase << ", E_{1#times1 cm^{2}} > ";

  const unsigned nV = 1;
  TString version[nV] = {"33"};
  //double Emip = 0.0548;//in MeV
  // double Emip[nLayers];
  // for (unsigned iL(0);iL<nLayers;++iL){
  //   if (nLayers < nEcalLayers) Emip[iL] = 0.0548;//in MeV
  //   else if (nLayers < nEcalLayers+24)  Emip[iL] = 0.0849;
  //   else Emip[iL] = 1.49;
  // }

  const unsigned nmips = 3;
  unsigned mipThresh = 1;//Base[nmips] = {1,5,10};
  bool doThresh = true;

  double minRadius = 316;//mm

  //VBFH
  const unsigned nEvts = 7;//1000;
  unsigned event[nEvts]={4,5,6,7,9,11,12};//,6,12};//103,659,875};
  //Htautau
  //const unsigned nEvts = 7;//1000;
  //unsigned event[nEvts]={1,2,5,8,10,11,12};//,6,12};//103,659,875};
  //Htautau 125
  //const unsigned nEvts = 2;//1000;
  //unsigned event[nEvts]={10,14};//,6,12};//103,659,875};
  //gluons
  //const unsigned nEvts = 3;//1000;
  //unsigned event[nEvts]={7,16,22};//,6,12};//103,659,875};


  //for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on events
  //  event[ievt] = ievt;
  // }

  const unsigned nLayers = 52;
  const unsigned nEcalLayers = 28;

  double minZ=3170,maxZ=5400;
  unsigned nZ=(maxZ-minZ)/2.;
  //double minX=-150,maxX=150;

  //VBFH
  double minX[nEvts] = {-700,-600,-400,-400,0,0,0};
  double maxX[nEvts] = {-100,0,200,100,700,700,500};
  double minY[nEvts] = {100,-600,-1100,-800,-100,-200,-700};
  double maxY[nEvts] = {480,0,0,0,400,300,0};

  double minXecal[nEvts] = {-600,-490,-400,-300,450,320,50};
  double maxXecal[nEvts] = {-350,-290,200,-100,650,520,350};
  double minYecal[nEvts] = {200,-490,-1100,-650,100,-100,-580};
  double maxYecal[nEvts] = {450,-290,0,-450,300,100,-280};

  //Htautau
  /*double minX[nEvts] = {0,-500,0,-400,0,-500,-700};
  double maxX[nEvts] = {600,0,700,500,700,0,600};
  double minY[nEvts] = {0,-700,-100,-700,-400,0,-400};
  double maxY[nEvts] = {500,0,400,900,200,700,700};

  double minXecal[nEvts] = {350,-400,450,-300,400,-400,-700};
  double maxXecal[nEvts] = {550,-200,650,-100,600,-200,600};
  double minYecal[nEvts] = {180,-550,20,600,-200,400,-400};
  double maxYecal[nEvts] = {380,-300,220,800,0,600,700};*/

  //Htautau 125
  /*double minX[nEvts] = {-300,-850,-200,-450,-550};
  double maxX[nEvts] = {100,-450,300,0,-150};
  double minY[nEvts] = {-900,-280,200,-950,-800};
  double maxY[nEvts] = {-500,150,650,-550,-400};

  double minXecal[nEvts] = {-190,-770,-30,-330,-460};
  double maxXecal[nEvts] = {10,-570,170,-130,-260};
  double minYecal[nEvts] = {-800,-160,320,-840,-720};
  double maxYecal[nEvts] = {-600,40,520,-640,-520};

  double minX[nEvts] = {-450,-550};
  double maxX[nEvts] = {50,-150};
  double minY[nEvts] = {-1100,-950};
  double maxY[nEvts] = {-700,-550};

  double minXecal[nEvts] = {-350,-480};
  double maxXecal[nEvts] = {-150,-280};
  double minYecal[nEvts] = {-900,-760};
  double maxYecal[nEvts] = {-700,-560};*/


  //pp to gg
  /*double minX[nEvts] = {-600,-250};
  double maxX[nEvts] = {0,350};
  double minY[nEvts] = {100,250};
  double maxY[nEvts] = {700,850};

  double minXecal[nEvts] = {-600,-250};
  double maxXecal[nEvts] = {0,350};
  double minYecal[nEvts] = {100,250};
  double maxYecal[nEvts] = {700,850};

  double minX[nEvts] = {100,-100,100};
  double maxX[nEvts] = {700,500,700};
  double minY[nEvts] = {100,200,100};
  double maxY[nEvts] = {700,800,700};
 
  double minXecal[nEvts] = {100,-100,100};
  double maxXecal[nEvts] = {700,500,700};
  double minYecal[nEvts] = {100,200,100};
  double maxYecal[nEvts] = {700,800,700};
  */
  //double minX=-700,maxX=700;
  //double minY=420,maxY=720;
  //double minY=1150,maxY=1450;
  //double minY=-1100,maxY=500;

  //unsigned nX[nEvts],nY[nEvts];

  bool doAll = false;

  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on events
    if (doAll){
      minX[ievt]=-1700;maxX[ievt]=1700;
      minY[ievt]=-1700;maxY[ievt]=1700;
    }
    //nX[ievt]=(maxX[ievt]-minX[ievt])/10;
    //nY[ievt]=(maxY[ievt]-minY[ievt])/10;
  }


  TCanvas *mycECAL = new TCanvas("mycECAL","mycECAL",1200,800);
  TCanvas *mycECALzoom = new TCanvas("mycECALzoom","mycECALzoom",1200,800);
  TCanvas *mycHCAL = new TCanvas("mycHCAL","mycHCAL",1200,800);
  //TCanvas *mycAll = new TCanvas("mycAll","mycAll",1200,800);
  //const unsigned nPads = nEvts>10 ? 10 : nEvts;
  //mycAll->Divide(static_cast<unsigned>(nPads/2.+0.5),nEvts/2<1 ? 1 : 2);

  const unsigned nCanvas = nS;  
  TCanvas *myc[nCanvas];
  for (unsigned iC(0);iC<nCanvas;++iC){
    std::ostringstream lName;
    lName << "myc" << iC;
    myc[iC] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1300,800);
  }

  std::ostringstream saveName;
  
  for (unsigned iV(0); iV<nV;++iV){//loop on versions
    for (unsigned iS(0); iS<nS;++iS){//loop on scenarios
      
      TString inputDir = "/afs/cern.ch/work/a/amagnan/public/HGCalGeant4/git_V00-03-06/version_"+version[iV]+"/model_2/"+scenario[iS]+"/BON/";
      
      unsigned evtcounter = 0;
      std::ostringstream lName;
      lName << inputDir << "DigiPFcal.root";//CalibHistos_E200_evt" << event[ievt] << ".root";
      TFile *inputFile = TFile::Open(lName.str().c_str());
      if (!inputFile) {
	std::cout << " -- Error, input file " << lName.str() << " cannot be opened." << std::endl;
	return 1;
      }
      inputFile->cd();
      TTree *tree = (TTree*)gDirectory->Get("RecoTree");
      if (!tree) {
	std::cout << " -- Error, tree not found." << std::endl;
	return 1;
      }
      TString plotDir = "PLOTS/version_"+version[iV]+"/"+scenario[iS]+"/";
      if (doAll) plotDir += "Overview/";

      for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on events
	std::cout << " -- Processing event " << event[ievt] << std::endl;
	
	std::ostringstream lcutSmall;
	lcutSmall << "(HGCSSRecoHitVec.energy_*absWeight(HGCSSRecoHitVec.layer_))*(event_==" << event[ievt] << " && (TMath::Sqrt(TMath::Power(HGCSSRecoHitVec.xpos_,2)+TMath::Power(HGCSSRecoHitVec.ypos_,2))<minR(HGCSSRecoHitVec.layer_)) && (TMath::Sqrt(TMath::Power(HGCSSRecoHitVec.xpos_,2)+TMath::Power(HGCSSRecoHitVec.ypos_,2))>" << minRadius << ") && (HGCSSRecoHitVec.energy_*absWeight(HGCSSRecoHitVec.layer_))>" << mipThresh << ")";
	std::ostringstream lcut;
	lcut << "(HGCSSRecoHitVec.energy_*absWeight(HGCSSRecoHitVec.layer_))*(event_==" << event[ievt] << " && (TMath::Sqrt(TMath::Power(HGCSSRecoHitVec.xpos_,2)+TMath::Power(HGCSSRecoHitVec.ypos_,2))>=minR(HGCSSRecoHitVec.layer_)) && (HGCSSRecoHitVec.energy_*absWeight(HGCSSRecoHitVec.layer_))>" << mipThresh << ")";

	TH2F *p_xy[nLayers];
	TH2F *p_xySmall[nLayers];
	double EmaxEcal = 0;
	double EmaxHcal = 0;

	unsigned nBins = 339;
	double min=-1695;
	double max=1695;
	std::ostringstream lvar;
	
	for (unsigned iL(0); iL<nLayers; ++iL){
	  lName.str("");
	  lName << "p_xy_" << iL;
	  nBins = 339;
	  p_xy[iL] =  new TH2F(lName.str().c_str(),
			       ";x(mm);y(mm)",
			       nBins,min,max,
			       nBins,min,max);
	  lvar.str("");
	  lvar << "HGCSSRecoHitVec.ypos_:HGCSSRecoHitVec.xpos_>>" << lName.str();
	  std::ostringstream lcutLayer;
	  lcutLayer << "(HGCSSRecoHitVec.energy_*absWeight(HGCSSRecoHitVec.layer_))*(event_==" << event[ievt] << " && HGCSSRecoHitVec.layer_==" << iL << " && (TMath::Sqrt(TMath::Power(HGCSSRecoHitVec.xpos_,2)+TMath::Power(HGCSSRecoHitVec.ypos_,2))>=minR(HGCSSRecoHitVec.layer_)) && (HGCSSRecoHitVec.energy_*absWeight(HGCSSRecoHitVec.layer_))>" << mipThresh << ")";
	  tree->Draw(lvar.str().c_str(),lcutLayer.str().c_str());

	  nBins = 452;
	  lName.str("");
	  lName << "p_xySmall_" << iL;
	  p_xySmall[iL] =  new TH2F(lName.str().c_str(),
				    ";x(mm);y(mm)",
				    nBins,min,max,
				    nBins,min,max);
	  
	  lvar.str("");
	  lvar << "HGCSSRecoHitVec.ypos_:HGCSSRecoHitVec.xpos_>>" << lName.str();
	  std::ostringstream lcutLayerSmall;
	  lcutLayer << "(HGCSSRecoHitVec.energy_*absWeight(HGCSSRecoHitVec.layer_))*(event_==" << event[ievt] << " && HGCSSRecoHitVec.layer_==" << iL << " && (TMath::Sqrt(TMath::Power(HGCSSRecoHitVec.xpos_,2)+TMath::Power(HGCSSRecoHitVec.ypos_,2))<minR(HGCSSRecoHitVec.layer_)) && (TMath::Sqrt(TMath::Power(HGCSSRecoHitVec.xpos_,2)+TMath::Power(HGCSSRecoHitVec.ypos_,2))>" << minRadius << ") && (HGCSSRecoHitVec.energy_*absWeight(HGCSSRecoHitVec.layer_))>" << mipThresh << ")";
	  tree->Draw(lvar.str().c_str(),lcutLayerSmall.str().c_str());

	  double Etot = p_xy[iL]->GetMaximum();
	  if (Etot > EmaxEcal && iL<nEcalLayers) EmaxEcal = Etot;
	  if (Etot > EmaxHcal && iL>=nEcalLayers) EmaxHcal = Etot;
	  
	}//loop on layers

	std::cout << " -- max E found: ECAL: " <<  EmaxEcal << ", HCAL: " << EmaxHcal << std::endl;

	nBins = 339;
	TH2F *p_xytot = new TH2F("p_xytot",";x(mm);y(mm)",
				 nBins,min,max,
				 nBins,min,max);

	lvar.str("");
	lvar << "HGCSSRecoHitVec.ypos_:HGCSSRecoHitVec.xpos_>>p_xytot";
	tree->Draw(lvar.str().c_str(),lcut.str().c_str());

	TH2F *p_xztot = new TH2F("p_xztot",";x(mm);z(mm)",
				 nBins,min,max,
				 nZ,minZ,maxZ);
	lvar.str("");
	lvar << "HGCSSRecoHitVec.zpos_:HGCSSRecoHitVec.xpos_>>p_xztot";
	tree->Draw(lvar.str().c_str(),lcut.str().c_str());

	TH2F *p_zytot = new TH2F("p_zytot",";z(mm);y(mm)",
				 nZ,minZ,maxZ,
				 nBins,min,max);
	lvar.str("");
	lvar << "HGCSSRecoHitVec.ypos_:HGCSSRecoHitVec.zpos_>>p_zytot";
	tree->Draw(lvar.str().c_str(),lcut.str().c_str());

	nBins = 452;
	TH2F *p_xytotSmall = new TH2F("p_xytotSmall",";x(mm);y(mm)",
				      nBins,min,max,
				      nBins,min,max);

	lvar.str("");
	lvar << "HGCSSRecoHitVec.ypos_:HGCSSRecoHitVec.xpos_>>p_xytotSmall";
	tree->Draw(lvar.str().c_str(),lcutSmall.str().c_str());

	TH2F *p_xztotSmall = new TH2F("p_xztotSmall",";x(mm);z(mm)",
				      nBins,min,max,
				      nZ,minZ,maxZ);
	lvar.str("");
	lvar << "HGCSSRecoHitVec.zpos_:HGCSSRecoHitVec.xpos_>>p_xztotSmall";
	tree->Draw(lvar.str().c_str(),lcutSmall.str().c_str());

	TH2F *p_zytotSmall = new TH2F("p_zytotSmall",";z(mm);y(mm)",
				      nZ,minZ,maxZ,
				      nBins,min,max);
	lvar.str("");
	lvar << "HGCSSRecoHitVec.ypos_:HGCSSRecoHitVec.zpos_>>p_zytotSmall";
	tree->Draw(lvar.str().c_str(),lcutSmall.str().c_str());


	set_plot_style();


	gStyle->SetOptStat(0);
	std::ostringstream ltitle;
	ltitle << ltitleBase.str() << mipThresh << " Mips" ;
	
	myc[iS]->Clear();
	myc[iS]->cd();
	gPad->SetLogz(1);
	p_xztot->GetXaxis()->SetLabelSize(0.05);
	p_xztot->GetYaxis()->SetLabelSize(0.05);
	p_xztot->GetXaxis()->SetTitleSize(0.05);
	p_xztot->GetYaxis()->SetTitleSize(0.05);
	p_xztot->GetZaxis()->SetTitleOffset(-0.5);
	p_xztot->GetZaxis()->SetTitle("E(MIP)");
	p_xztot->GetXaxis()->SetRangeUser(minX[ievt],maxX[ievt]);
	p_xztot->GetYaxis()->SetRangeUser(minZ,maxZ);
	p_xztot->Draw("colz");
	p_xztotSmall->Draw("colsame");
	//p_xztot->Draw("lego2z");
	//p_xztot->Draw("");
	myc[iS]->Update();
	saveName.str("");
	saveName << plotDir << "/xzSimHits_evt" << event[ievt] << "_mipThresh" << mipThresh;
	myc[iS]->Print((saveName.str()+".png").c_str());
	myc[iS]->Print((saveName.str()+".pdf").c_str());
	
	myc[iS]->cd();
	gPad->SetLogz(1);
	//p_zytot->RebinX(2);
	//p_zytot->RebinY(4);
	p_zytot->GetXaxis()->SetLabelSize(0.05);
	p_zytot->GetYaxis()->SetLabelSize(0.05);
	p_zytot->GetXaxis()->SetTitleSize(0.05);
	p_zytot->GetYaxis()->SetTitleSize(0.05);
	p_zytot->GetZaxis()->SetTitleOffset(-0.5);
	p_zytot->GetZaxis()->SetTitle("E(MIP)");
	p_zytot->GetXaxis()->SetRangeUser(minZ,maxZ);
	p_zytot->GetYaxis()->SetRangeUser(minY[ievt],maxY[ievt]);
	p_zytot->Draw("colz");
	p_zytotSmall->Draw("colsame");
	//p_zytot->Draw("lego2z");
	//p_zytot->Draw("");
	myc[iS]->Update();
	saveName.str("");
	saveName << plotDir << "/zySimHits_evt" << event[ievt] << "_mipThresh" << mipThresh;
	myc[iS]->Print((saveName.str()+".png").c_str());
	myc[iS]->Print((saveName.str()+".pdf").c_str());
	
	myc[iS]->cd();
	gPad->SetLogz(1);
	p_xytot->GetXaxis()->SetLabelSize(0.05);
	p_xytot->GetYaxis()->SetLabelSize(0.05);
	p_xytot->GetXaxis()->SetTitleSize(0.05);
	p_xytot->GetYaxis()->SetTitleSize(0.05);
	p_xytot->GetZaxis()->SetTitleOffset(-0.5);
	p_xytot->GetZaxis()->SetTitle("E(MIP)");
	p_xytot->GetXaxis()->SetRangeUser(minX[ievt],maxX[ievt]);
	p_xytot->GetYaxis()->SetRangeUser(minY[ievt],maxY[ievt]);
	//p_xytot->Draw("colz");
	p_xytot->Draw("lego2z");
	p_xytotSmall->Draw("lego2same");
	
	myc[iS]->Update();
	saveName.str("");
	saveName << plotDir << "/xySimHits_evt" << event[ievt] << "_mipThresh" << mipThresh;
	myc[iS]->Print((saveName.str()+".png").c_str());
	myc[iS]->Print((saveName.str()+".pdf").c_str());
	
	mycECAL->Clear();
	unsigned counter = 1;
	mycECAL->Divide(3,3);
	mycECALzoom->Clear();
	mycECALzoom->Divide(3,3);
	
	mycHCAL->Clear();
	mycHCAL->Divide(5,3);
	
	unsigned layId = 0;
	for (unsigned iL(1); iL<nLayers; ++iL){//loop on layers
	  //std::cout << " -- Processing layer " << iL << std::endl;
	  if ((iL-1)%3==2 && counter < 10 && iL<nEcalLayers) {
	    std::cout << " -- Ecal layer " << iL << " counter = " << counter << std::endl;
	    mycECAL->cd(counter);
	    counter++;
	    layId = iL;
	  }
	  else if (counter<25 && iL>=nEcalLayers)  {
	    std::cout << " -- Hcal layer " << iL << " counter = " << counter << std::endl;
	    mycHCAL->cd(counter-9);
	    layId = nEcalLayers+counter-10;
	    if (iL>=nEcalLayers+12) layId = nEcalLayers+15+counter-10;
	    counter++;
	  }
	  else if (iL==nEcalLayers+12+3){
	    saveName.str("");
	    saveName << plotDir << "/xySimHits_HCAL_evt" << event[ievt] << "_subset1";
	    mycHCAL->Update();
	    mycHCAL->Print((saveName.str()+".png").c_str());
	    mycHCAL->Print((saveName.str()+".pdf").c_str());
	    counter = 10;
	    std::cout << " -- Hcal layer " << iL << " counter = " << counter << std::endl;
	    mycHCAL->Clear();
	    mycHCAL->Divide(3,3);
	    mycHCAL->cd(counter-9);
	    layId = nEcalLayers+15+counter-10;
	    counter++;
	  }
	  else continue;
	  gPad->SetLogz(1);
	  p_xy[iL]->SetMaximum(iL<nEcalLayers ? EmaxEcal : EmaxHcal);
	  p_xy[iL]->GetXaxis()->SetLabelSize(0.05);
	  p_xy[iL]->GetYaxis()->SetLabelSize(0.05);
	  p_xy[iL]->GetXaxis()->SetTitleSize(0.05);
	  p_xy[iL]->GetYaxis()->SetTitleSize(0.05);
	  p_xy[iL]->GetZaxis()->SetTitleOffset(-0.5);
	  p_xy[iL]->GetZaxis()->SetTitle("E(MIP)");

	  //p_xy[iL]->Draw("colz");
	  char buf[500];
	  sprintf(buf,"Layer %d",layId-1);
	  p_xy[iL]->SetTitle(buf);

	  //p_xy[iL]->Draw("LEGO2z");
	  p_xy[iL]->GetXaxis()->SetRangeUser(minX[ievt],maxX[ievt]);
	  p_xy[iL]->GetYaxis()->SetRangeUser(minY[ievt],maxY[ievt]);

	  gPad->SetTheta(60); // default is 30
	  gPad->SetPhi(40); // default is 30
	  gPad->Update();
	  p_xy[iL]->Draw("lego2z");
	  p_xySmall[iL]->Draw("lego2same");

	}//loop on layers
	saveName.str("");
	saveName << plotDir << "/xySimHits_ECAL_evt" << event[ievt] << "_subset";
	mycECAL->Update();
	mycECAL->Print((saveName.str()+".png").c_str());
	mycECAL->Print((saveName.str()+".pdf").c_str());
	saveName.str("");
	saveName << plotDir << "/xySimHits_HCAL_evt" << event[ievt] << "_subset2";
	mycHCAL->Update();
	mycHCAL->Print((saveName.str()+".png").c_str());
	mycHCAL->Print((saveName.str()+".pdf").c_str());


	//plot ECAL zoom
	counter=1;
	for (unsigned iL(1); iL<nEcalLayers; ++iL){//loop on ecal layers
          //std::cout << " -- Processing layer " << iL << std::endl;
          if ((iL-1)%3==2 && counter < 10) {
	    std::cout << " -- ZOOM Ecal layer " << iL << " counter = " << counter << std::endl;
            mycECALzoom->cd(counter);
            counter++;
            layId = iL;
          }
          else continue;
          gPad->SetLogz(1);
          p_xy[iL]->SetMaximum(iL<nEcalLayers ? EmaxEcal : EmaxHcal);
          p_xy[iL]->GetXaxis()->SetLabelSize(0.05);
          p_xy[iL]->GetYaxis()->SetLabelSize(0.05);
          p_xy[iL]->GetXaxis()->SetTitleSize(0.05);
          p_xy[iL]->GetYaxis()->SetTitleSize(0.05);
          p_xy[iL]->GetZaxis()->SetTitleOffset(-0.5);
          p_xy[iL]->GetZaxis()->SetTitle("E(MIP)");

          //p_xy[iL]->Draw("colz");
          char buf[500];
          sprintf(buf,"Layer %d",layId-1);
          p_xy[iL]->SetTitle(buf);

          //p_xy[iL]->Draw("LEGO2z");
          p_xy[iL]->GetXaxis()->SetRangeUser(minXecal[ievt],maxXecal[ievt]);
          p_xy[iL]->GetYaxis()->SetRangeUser(minYecal[ievt],maxYecal[ievt]);
          p_xy[iL]->Draw("colz");
	  p_xySmall[iL]->Draw("colsame");

        }//loop on layers
        saveName.str("");
        saveName << plotDir << "/xySimHits_ECAL_evt" << event[ievt] << "_subset";
	saveName << "_zoom";
        mycECALzoom->Update();
        mycECALzoom->Print((saveName.str()+".png").c_str());
        mycECALzoom->Print((saveName.str()+".pdf").c_str());


	evtcounter++;

	return 1;

      }//loop on events
 
    }//loop on scenarios

  }//loop on versions
  
  return 0;


}//main
