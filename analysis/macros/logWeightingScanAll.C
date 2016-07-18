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
#include "TProfile.h"

#include "TDRStyle.h"

double calibratedE(const double Etot, const double eta){
  //calibration for signal region 2: 3*3 cm^2
  double pars[3] = {76.9,3.5,-0.53};
  double paro[3] = {-5.3,-12.8,-6.9};
  double offset = paro[0] + paro[1]*eta + paro[2]*eta*eta;
  double slope = pars[0] + pars[1]*eta + pars[2]*eta*eta;
  return (Etot-offset)/slope;
};

double absWeight(const unsigned layer, const double eta){
  double weight = 1.0;
  if (layer == 0) weight = 0.0378011;
  else if (layer == 1) weight = 1;
  else if (layer == 2) weight = 0.646989;
  else if (layer == 3) weight = 0.617619;
  else if (layer == 4) weight = 0.646989;
  else if (layer == 5) weight = 0.617619;
  else if (layer == 6) weight = 0.646989;
  else if (layer == 7) weight = 0.617619;
  else if (layer == 8) weight = 0.646989;
  else if (layer == 9) weight = 0.617619;
  else if (layer == 10) weight = 0.646989;
  else if (layer == 11) weight = 0.942829;
  else if (layer == 12) weight = 0.859702;
  else if (layer == 13) weight = 0.942829;
  else if (layer == 14) weight = 0.859702;
  else if (layer == 15) weight = 0.942829;
  else if (layer == 16) weight = 0.859702;
  else if (layer == 17) weight = 0.942829;
  else if (layer == 18) weight = 0.859702;
  else if (layer == 19) weight = 0.942829;
  else if (layer == 20) weight = 0.859702;
  else if (layer == 21) weight = 1.37644;
  else if (layer == 22) weight = 1.30447;
  else if (layer == 23) weight = 1.37644;
  else if (layer == 24) weight = 1.30447;
  else if (layer == 25) weight = 1.37644;
  else if (layer == 26) weight = 1.30447;
  else if (layer == 27) weight = 1.37644;
  else if (layer == 28) weight = 1.30447;
  else if (layer == 29) weight = 1.79662;
  else weight = 1;
  return weight/tanh(eta);
};


double zpos(const unsigned l){
  if (l==0) return 3173.9;
  if (l==1) return 3183.65;
  if (l==2) return 3190.95;
  if (l==3) return 3200.4;
  if (l==4) return 3207.7;
  if (l==5) return 3217.15;
  if (l==6) return 3224.45;
  if (l==7) return 3233.9;
  if (l==8) return 3241.2;
  if (l==9) return 3250.65;
  if (l==10) return 3257.95;
  if (l==11) return 3268.45;
  if (l==12) return 3276.85;
  if (l==13) return 3287.35;
  if (l==14) return 3295.75;
  if (l==15) return 3306.25;
  if (l==16) return 3314.65;
  if (l==17) return 3325.15;
  if (l==18) return 3333.55;
  if (l==19) return 3344.05;
  if (l==20) return 3352.45;
  if (l==21) return 3364.35;
  if (l==22) return 3375.05;
  if (l==23) return 3386.95;
  if (l==24) return 3397.65;
  if (l==25) return 3409.55;
  if (l==26) return 3420.25;
  if (l==27) return 3432.15;
  if (l==28) return 3442.85;
  if (l==29) return 3454.75;
  return 0;
};

int logWeightingScanAll(){//main
  SetTdrStyle();

  //TString plotDir = "../PLOTS/gitV00-02-12/version12/gamma/200um/";
  TString plotDir = "/afs/cern.ch/work/a/amagnan/PFCalEEAna/PLOTS/gitV00-02-12/version12/gamma/200um/";

  bool useFit = true;

  //const unsigned nPu = 2;
  //unsigned pu[nPu] = {0,140};
  const unsigned nPu = 1;//2;
  unsigned pu[nPu] = {0};//,140};

  const unsigned nScans = 50;
  const double wStart = 1.;
  const double wStep = (6.-wStart)/nScans;

  const unsigned neta = 4;//7;
  const unsigned npt = 3;//13;

  const unsigned nLayers = 30;

  unsigned eta[neta] = {17,21,25,29};
  //unsigned pt[npt] = {20,30,40,50,60,70,80,90,100,125,150,175,200};
  unsigned pt[npt] = {20,50,100};

  //TF1 *cauchy = new TF1("cauchy","1/(TMath::Pi()*[2]*(1+pow((x-[0])/[2],2)))",-15,15);
  //cauchy->SetParameters(0,0,2.5);


  double wxminall[nPu][nLayers][neta];
  double wyminall[nPu][nLayers][neta];
  double etaval[neta];
       
  const unsigned nCanvas = 5;  
  TCanvas *mycx[nCanvas];
  TCanvas *mycy[nCanvas];
  for (unsigned iC(0);iC<nCanvas;++iC){
    std::ostringstream lName;
    lName << "mycx" << iC;
    mycx[iC] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1500,1000);
    mycx[iC]->Divide(3,2);
    lName.str("");
    lName << "mycy" << iC;
    mycy[iC] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1500,1000);
    mycy[iC]->Divide(3,2);
  }
  
  TCanvas *mycW = new TCanvas("mycW","mycW",1500,1000);
  TCanvas *mycFit = new TCanvas("mycFit","mycFit",1);

  for (unsigned ieta(0); ieta<neta;++ieta){
    etaval[ieta] = eta[ieta]/10.;
    bool savePoint = (eta[ieta] == 17);// &&  pt[ipt]==50);
    std::ostringstream pteta;
    pteta << "eta" << eta[ieta];// << "_et" << pt[ipt];
      
    TFile *fin[nPu];
      
    TFile *fout = 0;
    if (savePoint) {
      fout = TFile::Open(("PLOTS/LogWeightingStudy_"+pteta.str()+".root").c_str(),"RECREATE");
      fout->mkdir("scan");
      fout->mkdir("scan/xpos");
      fout->mkdir("scan/ypos");
      for (unsigned iS(0); iS<nScans;++iS){
	std::ostringstream lName;
	lName << "scan/xpos/scan_" << wStart+iS*wStep;
	fout->mkdir(lName.str().c_str());
	lName.str("");
	lName << "scan/ypos/scan_" << wStart+iS*wStep;
	fout->mkdir(lName.str().c_str());
      }
	
      fout->cd();
    }
    TH1F *p_xt;
    TH1F *p_yt; 
    TH1F *p_intercalibSigmaSquare[nPu];
    TH1F *p_chi2ndf; 
    if (savePoint) {
      p_xt = new TH1F("p_xt",";x truth (mm)",100,-5,5); 
      p_yt = new TH1F("p_yt",";y truth (mm)",100,-5,5);//200,1170,1370); 
	
      p_intercalibSigmaSquare[0] = new TH1F("p_intercalibSigmaSquare_0",";#sigma_{E}^{2} (2% intercalib) (GeV^2)",3000,0,30000);
      p_intercalibSigmaSquare[1] = new TH1F("p_intercalibSigmaSquare_140",";#sigma_{E}^{2} (2% intercalib) (GeV^2)",3000,0,30000);
      p_chi2ndf = new TH1F("p_chi2ndf",";#chi^{2}/N",100,0,20);
    }
      
    TH1F *p_posx[nPu][nLayers][nScans];
    TH1F *p_posy[nPu][nLayers][nScans];
      
    TH1F *p_wx[nPu][nLayers][3];
    TH1F *p_wy[nPu][nLayers][3];
      
    TH2F *p_Exy[nPu][nLayers];
    TH2F *p_deltavsreco_x[nPu][nLayers];
    TH2F *p_deltavsreco_y[nPu][nLayers];
    TH2F *p_recovstruth_x[nPu][nLayers];
    TH2F *p_recovstruth_y[nPu][nLayers];
      
    if (savePoint){
      mycFit->Print("PLOTS/fits_x.pdf[");
      mycFit->Print("PLOTS/fits_y.pdf[");
    }
      
    gStyle->SetOptStat("eMRuo");
      
    TGraphErrors *grX[nPu][nLayers];
    TGraphErrors *grY[nPu][nLayers];
    TGraphErrors *grXrms[nPu][nLayers];
    TGraphErrors *grYrms[nPu][nLayers];
      
    double resxmin[nPu][nLayers];
    double resymin[nPu][nLayers];
    double lay[nLayers];
    double wxmin[nPu][nLayers];
    double wymin[nPu][nLayers];

    for (unsigned ipu(0); ipu<nPu; ++ipu){//loop on pu


      std::vector<std::vector<double> > Exy;
      std::vector<double> init;
      init.resize(25,0);
      Exy.resize(nLayers,init);
      std::vector<double> truthPosX;
      truthPosX.resize(nLayers,0);
      std::vector<double> truthPosY;
      truthPosY.resize(nLayers,0);
	  
	  
      std::ostringstream label;
      label << "pu" << pu[ipu];
      if (savePoint){
	fout->mkdir(label.str().c_str());
	fout->cd(label.str().c_str());
      }

      for (unsigned iL(0);iL<nLayers;++iL){
	lay[iL] = iL;
	resxmin[ipu][iL] = 100;
	wxmin[ipu][iL] = 10;
	resymin[ipu][iL] = 100;
	wymin[ipu][iL] = 10;
	label.str("");   
	label << "pu" << pu[ipu];
	if (savePoint) fout->cd(label.str().c_str());
	label.str("");   
	label << "grX_pu" << pu[ipu] << "_" << iL;
	grX[ipu][iL] = new TGraphErrors();
	grX[ipu][iL]->SetName(label.str().c_str());
	label.str("");   
	label << "grY_pu" << pu[ipu] << "_" << iL;
	grY[ipu][iL] = new TGraphErrors();
	grY[ipu][iL]->SetName(label.str().c_str());
	label.str("");   
	label << "grXrms_pu" << pu[ipu] << "_" << iL;
	grXrms[ipu][iL] = new TGraphErrors();
	grXrms[ipu][iL]->SetName(label.str().c_str());
	label.str("");   
	label << "grYrms_pu" << pu[ipu] << "_" << iL;
	grYrms[ipu][iL] = new TGraphErrors();
	grYrms[ipu][iL]->SetName(label.str().c_str());
	if (savePoint) {
	  label.str("");   
	  label << "Exy_pu"<< pu[ipu] << "_" << iL;
	  p_Exy[ipu][iL] = new TH2F(label.str().c_str(),";x idx;y idx; E (mips)",
				    5,0,5,5,0,5);
	}
	label.str("");   
	label << "deltavsreco_x_pu"<< pu[ipu] << "_" << iL;
	p_deltavsreco_x[ipu][iL] = new TH2F(label.str().c_str(),";x reco (mm);x_{reco}-x_{truth} (mm);",
					    30,-15,15,100,-10,10);
	label.str("");
	label << "deltavsreco_y_pu"<< pu[ipu] << "_" << iL;
	p_deltavsreco_y[ipu][iL] = new TH2F(label.str().c_str(),";y reco (mm);y_{reco}-y_{truth} (mm);",
					    30,-15,15,100,-10,10);
	if (savePoint) {
	  label.str("");   
	  label << "recovstruth_x_pu"<< pu[ipu] << "_" << iL;
	  p_recovstruth_x[ipu][iL] = new TH2F(label.str().c_str(),";x_{truth} (mm);x reco (mm)",
					      30,-15,15,30,-15,15);
	  label.str("");
	  label << "recovstruth_y_pu"<< pu[ipu] << "_" << iL;
	  p_recovstruth_y[ipu][iL] = new TH2F(label.str().c_str(),";y_{truth} (mm);y reco (mm)",
					      30,-15,15,//200,1170,1370,
					      30,-15,15
					      ); 
	    
	}


	if (savePoint) {
	  if (savePoint) {
	    label.str("");   
	    label << "pu" << pu[ipu];
	    fout->cd(label.str().c_str());
	  }
	  for (unsigned i(0);i<5;++i){
	    label.str("");     
	    label << "wx_pu" << pu[ipu] << "_" << iL << "_" << i;
	    p_wx[ipu][iL][i] = new TH1F(label.str().c_str(),
					";wx;events",
					100,-10,0);
	    label.str("");     
	    label << "wy_pu" << pu[ipu] << "_" << iL << "_" << i;
	    p_wy[ipu][iL][i] = new TH1F(label.str().c_str(),
					";wy;events",
					100,-10,0);
	  }
	}
	for (unsigned iS(0); iS<nScans;++iS){
	  label.str("");
	  label << "scan/xpos/scan_" << wStart+iS*wStep;
	  if (savePoint) fout->cd(label.str().c_str());
	  label.str("");
	  label << "posx_pu" << pu[ipu] << "_" << iL << "_" << iS;
	  p_posx[ipu][iL][iS] = new TH1F(label.str().c_str(),
					 ";x-x_{truth} (mm);events",
					 100,-10,10);
	  label.str("");
	  label << "scan/ypos/scan_" << wStart+iS*wStep;
	  if (savePoint) fout->cd(label.str().c_str());
	  label.str("");
	  label << "posy_pu" << pu[ipu] << "_" << iL << "_" << iS;
	  p_posy[ipu][iL][iS] = new TH1F(label.str().c_str(),
					 ";y-y_{truth} (mm);events",
					 100,-10,10);
	}
      }//loop on layers
	  
      for (unsigned ipt(0); ipt<npt;++ipt){
	
	std::ostringstream linputStr;
	linputStr << plotDir << "/" << "eta" << eta[ieta] << "_et" << pt[ipt] << "_pu" << pu[ipu] ;
	//linputStr << "_logweight";
	linputStr << ".root";
	fin[ipu] = TFile::Open(linputStr.str().c_str());
	if (!fin[ipu]) {
	  std::cout << " -- Error, input file " << linputStr.str() << " cannot be opened. Skipping..." << std::endl;
	  continue;
	}
	else std::cout << " -- File " << linputStr.str() << " successfully opened." << std::endl;
	    
	TTree *tree = (TTree*)gDirectory->Get("EcellsSR2");
	if (!tree) {
	  std::cout << " Tree not found! " << std::endl;
	  continue;
	  //return 1;
	}

	for (unsigned iL(0);iL<nLayers;++iL){
	  label.str("");     
	  label << "TruthPosX_" << iL;
	  tree->SetBranchAddress(label.str().c_str(),&truthPosX[iL]);
	  label.str("");     
	  label << "TruthPosY_" << iL;
	  tree->SetBranchAddress(label.str().c_str(),&truthPosY[iL]);
	      
	  for (unsigned iy(0);iy<5;++iy){
	    for (unsigned ix(0);ix<5;++ix){
	      unsigned idx = 5*iy+ix;
	      label.str("");     
	      label << "E_" << iL << "_" << idx;
	      tree->SetBranchAddress(label.str().c_str(),&Exy[iL][idx]);
	    }
	  }
	}//loop on layers
	
	unsigned nEvts = tree->GetEntries();
	for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
	      
	  if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;
	      
	  tree->GetEntry(ievt);
	      
	  double Etotsq = 0;
	      
	  for (unsigned iL(0);iL<nLayers;++iL){
	    double Etot = 0;
	    double Ex[5] = {0,0,0,0,0};
	    double Ey[5] = {0,0,0,0,0};
	    for (unsigned idx(0);idx<25;++idx){
	      if (iL>22) Etot += Exy[iL][idx];
	      else if ((idx>5 && idx<9)||
		       (idx>10 && idx<14)||
		       (idx>15 && idx<19)) Etot += Exy[iL][idx];
	    }
	    Etotsq += pow(calibratedE(Etot*absWeight(iL,eta[ieta]/10.),eta[ieta]/10.),2);

	    if (iL>22){
	      Ex[0] = Exy[iL][0]+Exy[iL][5]+Exy[iL][10]+Exy[iL][15]+Exy[iL][20];
	      Ex[1] = Exy[iL][1]+Exy[iL][6]+Exy[iL][11]+Exy[iL][16]+Exy[iL][21];
	      Ex[2] = Exy[iL][2]+Exy[iL][7]+Exy[iL][12]+Exy[iL][17]+Exy[iL][22];
	      Ex[3] = Exy[iL][3]+Exy[iL][8]+Exy[iL][13]+Exy[iL][18]+Exy[iL][23];
	      Ex[4] = Exy[iL][4]+Exy[iL][9]+Exy[iL][14]+Exy[iL][19]+Exy[iL][24];
	      Ey[0] = Exy[iL][0]+Exy[iL][1]+Exy[iL][2]+Exy[iL][3]+Exy[iL][4];
	      Ey[1] = Exy[iL][5]+Exy[iL][6]+Exy[iL][7]+Exy[iL][8]+Exy[iL][9];
	      Ey[2] = Exy[iL][10]+Exy[iL][11]+Exy[iL][12]+Exy[iL][13]+Exy[iL][14];
	      Ey[3] = Exy[iL][15]+Exy[iL][16]+Exy[iL][17]+Exy[iL][18]+Exy[iL][19];
	      Ey[4] = Exy[iL][20]+Exy[iL][21]+Exy[iL][22]+Exy[iL][23]+Exy[iL][24];
	    }
	    else {
	      Ex[0] = Exy[iL][6]+Exy[iL][11]+Exy[iL][16];
	      Ex[1] = Exy[iL][7]+Exy[iL][12]+Exy[iL][17];
	      Ex[2] = Exy[iL][8]+Exy[iL][13]+Exy[iL][18];
	      Ey[0] = Exy[iL][6]+Exy[iL][7]+Exy[iL][8];
	      Ey[1] = Exy[iL][11]+Exy[iL][12]+Exy[iL][13];
	      Ey[2] = Exy[iL][16]+Exy[iL][17]+Exy[iL][18];
	    }
	    
	    double simplex = 0;
	    double simpley = 0;
	    if (Etot!=0) {
	      if (iL>22) {
		simplex = 10*(2*Ex[4]+Ex[3]-Ex[1]-2*Ex[0])/Etot;
		simpley = 10*(2*Ey[4]+Ey[3]-Ey[1]-2*Ey[0])/Etot;
	      }
	      else {
		simplex = 10*(Ex[2]-Ex[0])/Etot;
		simpley = 10*(Ey[2]-Ey[0])/Etot;
	      }
	    }
		
	    double xt = truthPosX[iL];
	    unsigned cellCenter = static_cast<unsigned>((truthPosY[iL]+5)/10.)*10;
	    double yt = 0;
	    //if (cellCenter>truthPosY[iL]) yt = cellCenter-truthPosY[iL];
	    yt=truthPosY[iL]-cellCenter;
	    p_deltavsreco_x[ipu][iL]->Fill(simplex,simplex-xt);
	    p_deltavsreco_y[ipu][iL]->Fill(simpley,simpley-yt);
	    if (savePoint) {
	      p_xt->Fill(xt);
	      p_yt->Fill(yt);
	      p_recovstruth_x[ipu][iL]->Fill(xt,simplex);
	      p_recovstruth_y[ipu][iL]->Fill(yt,simpley);
	
	      for (unsigned idx(0);idx<9;++idx){
		p_Exy[ipu][iL]->Fill(idx%5,idx/5,Exy[iL][idx]/Etot);
	      }
	    }
	    double wx[6][nScans];
	    double wy[6][nScans];
	  
	    for (unsigned i(0);i<6;++i){
	      for (unsigned iS(0); iS<nScans;++iS){
		wx[i][iS] = 0;
		wy[i][iS] = 0;
	      }
	    }
	    for (unsigned i(0);i<5;++i){
	      if (savePoint) p_wx[ipu][iL][i]->Fill(log(Ex[i]/Etot));
	      if (savePoint) p_wy[ipu][iL][i]->Fill(log(Ey[i]/Etot));
	      for (unsigned iS(0); iS<nScans;++iS){
		double w0 = wStart+iS*wStep;
		wx[i][iS] = std::max(0.,log(Ex[i]/Etot)+w0);
		wy[i][iS] = std::max(0.,log(Ey[i]/Etot)+w0);
		// if (log(Ex[i]/Etot)+w0<0)
		//   std::cout << " - iL= " << iL << " i=" << i << " w0=" << w0 
		// 		<< " logEx=" << log(Ex[i]/Etot)
		// 		<< " wx " << wx[i][iS] 
		// 		<< std::endl;
		// if (log(Ey[i]/Etot)+w0<0) 
		//   std::cout << " - iL= " << iL << " i=" << i << " w0=" << w0 
		// 		<< " logEy=" << log(Ey[i]/Etot)
		// 		<< " wy " << wy[i][iS] 
		// 		<< std::endl;
		wx[5][iS] += wx[i][iS];
		wy[5][iS] += wy[i][iS];
	      }
	    }
	    for (unsigned iS(0); iS<nScans;++iS){
	      double x = 0;//10*(wx[2][iS]-wx[0][iS])/wx[3][iS];
	      if (wx[5][iS]!=0) {
		if (iL>22) x = 10*(2*wx[4][iS]+wx[3][iS]-wx[1][iS]-2*wx[0][iS])/wx[5][iS];
		else x = 10*(wx[2][iS]-wx[0][iS])/wx[5][iS];
	      }
	      double y = 0;//10*(wy[2][iS]-wy[0][iS])/wy[3][iS];
	      if (wy[5][iS]!=0) {
		if (iL>22) y = 10*(2*wy[4][iS]+wy[3][iS]-wy[1][iS]-2*wy[0][iS])/wy[5][iS];
		else y = 10*(wy[2][iS]-wy[0][iS])/wy[5][iS];
	      }

	      //if (fabs(y-yt)>5) std::cout << " --- iL=" << iL << " iS=" << iS 
	      //<< " x=" << x << " xt=" << xt 
	      //<< " y=" << y << " yt=" << yt 
	      //<< std::endl;
	      p_posx[ipu][iL][iS]->Fill(x-xt);
	      p_posy[ipu][iL][iS]->Fill(y-yt);
	    }

	  }//loop on layers
	  if (savePoint) p_intercalibSigmaSquare[ipu]->Fill(Etotsq);
	}//loop on entries
    
      }//loop on pt

      //fill first point with linear weighting
      TLatex lat;
      char buf[500];

      for (unsigned iL(0);iL<nLayers;++iL){
	mycFit->cd();
	TH1D *projy = p_deltavsreco_x[ipu][iL]->ProjectionY();
	projy->Draw();
	projy->Fit("gaus","0+Q");
	TF1 *fitx = projy->GetFunction("gaus");
	fitx->SetLineColor(6);
	fitx->Draw("same");
	sprintf(buf,"Layer %d, linear weighting, pu=%d",iL,pu[ipu]);
	lat.DrawLatexNDC(0.1,0.96,buf);
	
	grX[ipu][iL]->SetPoint(0,0.5,fitx->GetParameter(2));
	grX[ipu][iL]->SetPointError(0,0,fitx->GetParError(2));
	grXrms[ipu][iL]->SetPoint(0,0.5,projy->GetRMS());
	grXrms[ipu][iL]->SetPointError(0,0,projy->GetRMSError());
	mycFit->Update();
	if (savePoint) mycFit->Print("PLOTS/fits_x.pdf");
	
	mycFit->cd();
	projy = p_deltavsreco_y[ipu][iL]->ProjectionY();
	projy->Draw();
	projy->Fit("gaus","0+Q");
	TF1 *fity = projy->GetFunction("gaus");
	fitx->SetLineColor(6);
	fitx->Draw("same");
	sprintf(buf,"Layer %d, linear weighting, pu=%d",iL,pu[ipu]);
	lat.DrawLatexNDC(0.1,0.96,buf);
	
	grY[ipu][iL]->SetPoint(0,0.5,fity->GetParameter(2));
	grY[ipu][iL]->SetPointError(0,0,fity->GetParError(2));
	grYrms[ipu][iL]->SetPoint(0,0.5,projy->GetRMS());
	grYrms[ipu][iL]->SetPointError(0,0,projy->GetRMSError());
	mycFit->Update();
	if (savePoint) mycFit->Print("PLOTS/fits_y.pdf");
      }

      //fit vs w0, get w0min
      for (unsigned iL(0);iL<nLayers;++iL){
	for (unsigned iS(0); iS<nScans;++iS){
	  mycFit->cd();
	  p_posx[ipu][iL][iS]->Draw();
	  double w0 = wStart+iS*wStep;
	  //myGaus->SetParameters();
	  p_posx[ipu][iL][iS]->Fit("gaus","0+Q","",-2.,2.);
	  TF1 *fitx = p_posx[ipu][iL][iS]->GetFunction("gaus");
	  fitx->SetLineColor(6);
	  fitx->Draw("same");
	  sprintf(buf,"Layer %d, w0=%3.1f, pu=%d",iL,w0,pu[ipu]);
	  lat.DrawLatexNDC(0.1,0.96,buf);

	  mycFit->Update();
	  if (savePoint) mycFit->Print("PLOTS/fits_x.pdf");

	  mycFit->cd();
	  p_posy[ipu][iL][iS]->Draw();
	  p_posy[ipu][iL][iS]->Fit("gaus","0+Q","",-2.,2.);
	  TF1 *fity = p_posy[ipu][iL][iS]->GetFunction("gaus");
	  fity->SetLineColor(6);
	  fity->Draw("same");
	  sprintf(buf,"Layer %d, w0=%3.1f, pu=%d",iL,w0,pu[ipu]);
	  lat.DrawLatexNDC(0.1,0.96,buf);
	  mycFit->Update();
	  if (savePoint) mycFit->Print("PLOTS/fits_y.pdf");
	  //grX[ipu][iL]->SetPoint(iS,w0,p_posx[ipu][iL][iS]->GetRMS());
	  //grX[ipu][iL]->SetPointError(iS,0,p_posx[ipu][iL][iS]->GetRMSError());
	  //grY[ipu][iL]->SetPoint(iS,w0,p_posy[ipu][iL][iS]->GetRMS());
	  //grY[ipu][iL]->SetPointError(iS,0,p_posy[ipu][iL][iS]->GetRMSError());
	  double xval = useFit? fitx->GetParameter(2) : p_posx[ipu][iL][iS]->GetRMS();
	  if (savePoint) p_chi2ndf->Fill(fitx->GetChisquare()/fitx->GetNDF());
	
	  grX[ipu][iL]->SetPoint(iS+1,w0,fitx->GetParameter(2));
	  grX[ipu][iL]->SetPointError(iS+1,0,fitx->GetParError(2));
	  grXrms[ipu][iL]->SetPoint(iS+1,w0,p_posx[ipu][iL][iS]->GetRMS());
	  grXrms[ipu][iL]->SetPointError(iS+1,0,p_posx[ipu][iL][iS]->GetRMSError());
	  if (xval < resxmin[ipu][iL]){
	    resxmin[ipu][iL] = xval;
	    wxminall[ipu][iL][ieta] = w0;
	    wxmin[ipu][iL] = w0;
	  }
	  double yval = useFit? fity->GetParameter(2) : p_posy[ipu][iL][iS]->GetRMS();
	  if (savePoint) p_chi2ndf->Fill(fity->GetChisquare()/fity->GetNDF());
	  grY[ipu][iL]->SetPoint(iS+1,w0,fity->GetParameter(2));
	  grY[ipu][iL]->SetPointError(iS+1,0,fity->GetParError(2));
	  grYrms[ipu][iL]->SetPoint(iS+1,w0,p_posy[ipu][iL][iS]->GetRMS());
	  grYrms[ipu][iL]->SetPointError(iS+1,0,p_posy[ipu][iL][iS]->GetRMSError());
	
	  if (yval < resymin[ipu][iL]){
	    resymin[ipu][iL] = yval;
	    wyminall[ipu][iL][ieta] = w0;
	    wymin[ipu][iL] = w0;
	  }
	}
	grX[ipu][iL]->SetTitle(";W0; #sigma(x-xt) (mm)");
	grY[ipu][iL]->SetTitle(";W0; #sigma(y-yt) (mm)");
	grX[ipu][iL]->SetLineColor(ipu+1);
	grX[ipu][iL]->SetMarkerColor(ipu+1);
	grX[ipu][iL]->SetMarkerStyle(ipu+21);
	grY[ipu][iL]->SetLineColor(ipu+1);
	grY[ipu][iL]->SetMarkerColor(ipu+1);
	grY[ipu][iL]->SetMarkerStyle(ipu+21);

	grXrms[ipu][iL]->SetLineColor(ipu+3);
	grXrms[ipu][iL]->SetMarkerColor(ipu+3);
	grXrms[ipu][iL]->SetMarkerStyle(ipu+23);
	grYrms[ipu][iL]->SetLineColor(ipu+3);
	grYrms[ipu][iL]->SetMarkerColor(ipu+3);
	grYrms[ipu][iL]->SetMarkerStyle(ipu+23);

	mycx[iL/6]->cd(iL%6+1);
	grX[ipu][iL]->Draw(ipu==0?"APL":"PLsame");
	grXrms[ipu][iL]->Draw("PLsame");
	//gStyle->SetStatX(0.4);
	//gStyle->SetStatY(1.0);
	//p_Exy[ipu][iL]->Draw("colztext");
	sprintf(buf,"Layer %d",iL);
	lat.DrawLatexNDC(0.4,0.85,buf);
	mycy[iL/6]->cd(iL%6+1);
	grY[ipu][iL]->Draw(ipu==0?"APL":"PLsame");
	grYrms[ipu][iL]->Draw("PLsame");
	lat.DrawLatexNDC(0.4,0.85,buf);

      }//loop on layers

      //return 1;
    }//loop on pu

    if (savePoint) {
      mycFit->Print("PLOTS/fits_x.pdf]");
      mycFit->Print("PLOTS/fits_y.pdf]");
    }

    return 1;

    TLegend *leg = new TLegend(0.6,0.6,0.94,0.94);
    leg->SetFillColor(0);
    if (grX[0][10]) leg->AddEntry(grX[0][10],"fit pu=0","P");
    if (grXrms[0][10]) leg->AddEntry(grXrms[0][10],"RMS pu=0","P");
    if (grX[1][10]) leg->AddEntry(grX[1][10],"fit pu=140","P");
    if (grXrms[1][10]) leg->AddEntry(grXrms[1][10],"RMS pu=140","P");
    for (unsigned iC(0);iC<nCanvas;++iC){
      mycx[iC]->cd(1);
      leg->Draw("same");
      mycx[iC]->Update();
      std::ostringstream lsave;
      lsave << "PLOTS/logWeighted_x_" << 6*iC << "_" << 6*iC+5 << "_" << pteta.str() ;
      lsave << ".pdf";
      mycx[iC]->Print(lsave.str().c_str());
      mycx[iC]->Update();
	
      mycy[iC]->cd(1);
      leg->Draw("same");
      lsave.str("");
      lsave << "PLOTS/logWeighted_y_" << 6*iC << "_" << 6*iC+5 << "_" << pteta.str();
      lsave << ".pdf";
      mycy[iC]->Print(lsave.str().c_str());
    }
      
    //plot w0min vs layer
    TGraph *grW[4] = {0,0,0,0};
    for (unsigned ipu(0); ipu<nPu; ++ipu){//loop on pu
      std::cout << " --Processing pu " << pu[ipu] << std::endl;
	
      grW[2*ipu] = new TGraph(nLayers,lay,wxmin[ipu]);
      grW[2*ipu+1] = new TGraph(nLayers,lay,wymin[ipu]);
      for (unsigned iL(0);iL<nLayers;++iL){
	std::cout << " if (layer==" << iL 
		  << ") return " << (wxmin[ipu][iL]+wymin[ipu][iL])/2. <<";"
	  //		<< " minimum x= " << grX[ipu][iL]->GetYaxis()->GetXmin()
	  //		<< " minimum y= " << grY[ipu][iL]->GetYaxis()->GetXmin()
		  << std::endl;
      }
      mycW->cd();
      grW[2*ipu]->SetTitle(";layer;W0");
      grW[2*ipu]->SetMaximum(6);
      grW[2*ipu]->SetMinimum(0);
      grW[2*ipu]->SetLineColor(2*ipu+1);
      grW[2*ipu]->SetMarkerColor(2*ipu+1);
      grW[2*ipu]->SetMarkerStyle(20+2*ipu+1);
      grW[2*ipu+1]->SetLineColor(2*ipu+2);
      grW[2*ipu+1]->SetMarkerColor(2*ipu+2);
      grW[2*ipu+1]->SetMarkerStyle(20+2*ipu+2);
      if (ipu==0) {
	grW[2*ipu]->Draw("APL");
	grW[2*ipu+1]->Draw("PLsame");
      }
      else {
	grW[2*ipu]->Draw("PLsame");
	grW[2*ipu+1]->Draw("PLsame");
      }
    }
    std::ostringstream lsave;
    lsave << "PLOTS/w0minvsLayers_" << pteta.str();
    if (useFit) lsave << "_fit";
    else lsave << "_rms";
    lsave << ".pdf";
    mycW->Update();
    mycW->Print(lsave.str().c_str());

    if (savePoint) fout->Write();
    else {
      for (unsigned i(0);i<4;++i){
	if (grW[i]) grW[i]->Delete();
      }
      for (unsigned ipu(0); ipu<nPu; ++ipu){//loop on pu
	for (unsigned iL(0); iL<nLayers;++iL){//loop on layers
	  grX[ipu][iL]->Delete();
	  grXrms[ipu][iL]->Delete();
	  grY[ipu][iL]->Delete();
	  grYrms[ipu][iL]->Delete();
	  p_deltavsreco_x[ipu][iL]->Delete();
	  p_deltavsreco_y[ipu][iL]->Delete();
	  for (unsigned iS(0); iS<nScans;++iS){
	    p_posx[ipu][iL][iS]->Delete();
	    p_posy[ipu][iL][iS]->Delete();
	  }
	}
      }
    }

    std::cout << " -- eta point finished successfully" << std::endl;

  }//loop on eta

  //plot w0min vs pt for all eta
  TGraph *grW[4] = {0,0,0,0};
  mycW->Clear();
  for (unsigned iL(0); iL<nLayers;++iL){//loop on layers
    for (unsigned ipu(0); ipu<nPu; ++ipu){//loop on pu
      std::cout << " --Processing pu " << pu[ipu] << std::endl;
      grW[2*ipu] = new TGraph(neta,etaval,wxminall[ipu][iL]);
      grW[2*ipu+1] = new TGraph(neta,etaval,wyminall[ipu][iL]);
      mycW->cd();//ieta+1);
      grW[2*ipu]->SetTitle(";#eta;W0");
      grW[2*ipu]->SetMaximum(6);
      grW[2*ipu]->SetMinimum(0);
      grW[2*ipu]->SetLineColor(2*ipu+1);
      grW[2*ipu]->SetMarkerColor(2*ipu+1);
      grW[2*ipu]->SetMarkerStyle(20+2*ipu+1);
      grW[2*ipu+1]->SetLineColor(2*ipu+2);
      grW[2*ipu+1]->SetMarkerColor(2*ipu+2);
      grW[2*ipu+1]->SetMarkerStyle(20+2*ipu+2);
      if (ipu==0) {
	grW[2*ipu]->Draw("APL");
	grW[2*ipu+1]->Draw("PLsame");
      }
      else {
	grW[2*ipu]->Draw("PLsame");
	grW[2*ipu+1]->Draw("PLsame");
      }
    }//loop on pu
    std::ostringstream lsave;
    lsave << "PLOTS/w0minvseta_layer" << iL;
    if (useFit) lsave << "_fit";
    else lsave << "_rms";
    lsave << ".pdf";
    
    char buf[500];
    TLatex lat;
    sprintf(buf,"Layer %d",iL);
    lat.DrawLatexNDC(0.1,0.96,buf);
    
    mycW->Update();
    mycW->Print(lsave.str().c_str());
  }//loop on layers
  
  
  return 0;
}//main
