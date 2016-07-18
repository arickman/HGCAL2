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
#include "TRandom3.h"

#include "TDRStyle.h"

struct FitResult{
  double chi2;
  unsigned ndf;
  double mean;
  double sigma;
  double meanerr;
  double sigmaerr;
};

double E(const unsigned pT, const unsigned eta){
  return pT*cosh(eta/10.);
};

double pT(const unsigned E, const unsigned eta){
  return E/cosh(eta/10.);
};

double absWeight(const unsigned layer, const double eta){
  //already added to tree variables
  return 1.;
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


TPad* plot_ratio(TPad *canv, bool up){
  canv->SetFillColor      (0);
  canv->SetBorderMode     (0);
  canv->SetBorderSize     (10);
  // Set margins to reasonable defaults
  canv->SetLeftMargin     (0.18);
  canv->SetLeftMargin     (0.17);
  canv->SetRightMargin    (0.05);
  canv->SetTopMargin      (0.12);
  canv->SetBottomMargin   (0.18);
  // Setup a frame which makes sense
  canv->SetFrameFillStyle (0);
  canv->SetFrameLineStyle (0);
  canv->SetFrameBorderMode(0);
  canv->SetFrameBorderSize(10);
  canv->SetFrameFillStyle (0);
  canv->SetFrameLineStyle (0);
  canv->SetFrameBorderMode(0);
  canv->SetFrameBorderSize(10);      

  canv->cd();
  TPad *pad = 0;
  if (up){
    pad = new TPad("upper","pad",0, 0.26 ,1 ,1);
    pad->SetBottomMargin(0.05);
    pad->SetTopMargin(0.09);
    pad->Draw();
    pad->cd();
    return pad;
  }
  else {
    pad = new TPad("lower","pad",0, 0   ,1 ,0.26);  
    pad->SetTopMargin(0.05);
    pad->SetBottomMargin(0.24);
    pad->Draw();
    return pad;
  }

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

  char buf[500];
  TLatex lat;
  double latx = hist->GetXaxis()->GetXmin()+(hist->GetXaxis()->GetXmax()-hist->GetXaxis()->GetXmin())/20.;
  double laty = hist->GetMaximum();
  sprintf(buf,"<Efit> = %3.3f +/- %3.3f %s",fitResult->GetParameter(1),fitResult->GetParError(1),unitStr.c_str());
  lat.DrawLatex(latx,laty*0.9,buf);
  sprintf(buf,"RMSfit = %3.3f +/- %3.3f %s",fitResult->GetParameter(2),fitResult->GetParError(2),unitStr.c_str());
  lat.DrawLatex(latx,laty*0.8,buf);
  sprintf(buf,"RMS/meanfit = %3.3f",fitResult->GetParameter(2)/fitResult->GetParameter(1));
  lat.DrawLatex(latx,laty*0.7,buf);
  
  sprintf(buf,"#chi^{2}/N = %3.3f/%d = %3.3f",fitResult->GetChisquare(),fitResult->GetNDF(),fitResult->GetChisquare()/fitResult->GetNDF());
  lat.DrawLatex(latx,laty*0.6,buf);
  
  lres.chi2 = fitResult->GetChisquare();
  lres.ndf = fitResult->GetNDF();
  lres.mean = fitResult->GetParameter(1);
  lres.meanerr = fitResult->GetParError(1);
  lres.sigma = fitResult->GetParameter(2);
  lres.sigmaerr = fitResult->GetParError(2);

  return 0;
};


TPad* plotCalibration(TGraphErrors *gr,TPad *pad,bool doRatio, TGraphErrors *grDelta,std::string unit, double & calib,double & calibErr, double & offset, double & offsetErr,const unsigned eta, const bool dovsE){

  TPad *upper = 0;
  TPad *lower = 0;

  if (!doRatio) pad->cd();
  else {
    pad->Clear();
    upper = plot_ratio(pad, true);
    lower = plot_ratio(pad, false);
    upper->cd();
  }
  gr->GetXaxis()->SetLabelSize(0.0);
  gr->GetXaxis()->SetTitleSize(0.0);
  gr->GetYaxis()->SetLabelSize(0.06);
  gr->GetYaxis()->SetTitleSize(0.06);
  gr->GetXaxis()->SetTitleOffset(7);
  gr->GetYaxis()->SetTitleOffset(0.9);
  
  gr->Draw("ap");

  if (!doRatio){
    if (!dovsE) gr->GetXaxis()->SetTitle("E_{T} (GeV)");
    else gr->GetXaxis()->SetTitle("E (GeV)");
  }
  else gr->GetXaxis()->SetTitle("");
  
  gr->GetYaxis()->SetTitle(("Average energy deposited ("+unit+")").c_str()); 
  char buf[500];
  TF1 *fitFunc=new TF1("calib","[0]+[1]*x",gr->GetXaxis()->GetXmin(),gr->GetXaxis()->GetXmax());
  fitFunc->SetLineColor(6);

  //if (dovsE) gr->Fit(fitFunc,"RIME","same",0,200);
  //else gr->Fit(fitFunc,"RIME","same",0,pT(200,eta));
  if (dovsE) gr->Fit(fitFunc,"IME","same");
  else gr->Fit(fitFunc,"IME","same");
  TLatex lat;
  lat.SetTextColor(6);
  lat.SetTextSize(0.1);
  if (!dovsE) sprintf(buf,"<E> #propto a + b #times E_{T} ");
  else sprintf(buf,"<E> #propto a + b #times E ");
  lat.DrawLatexNDC(0.2,0.85,buf);
  sprintf(buf,"a = %3.3f #pm %3.3f %s",fitFunc->GetParameter(0),fitFunc->GetParError(0),unit.c_str());
  lat.DrawLatexNDC(0.2,0.7,buf);
  sprintf(buf,"b = %3.3f #pm %3.3f %s/GeV",fitFunc->GetParameter(1),fitFunc->GetParError(1),unit.c_str());
  lat.DrawLatexNDC(0.2,0.55,buf);
  sprintf(buf,"#chi^{2}/N = %3.3f/%d = %3.3f",fitFunc->GetChisquare(),fitFunc->GetNDF(),fitFunc->GetChisquare()/fitFunc->GetNDF());
  lat.DrawLatexNDC(0.2,0.4,buf);

  calib = fitFunc->GetParameter(1);
  offset = fitFunc->GetParameter(0);
  calibErr = fitFunc->GetParError(1);
  offsetErr = fitFunc->GetParError(0);

  if (doRatio){
    //draw deltaE/E vs E
    lower->cd();
    gPad->SetLogx(0);
    gPad->SetGridx(1);
    gPad->SetGridy(1);

    double loffset = fitFunc->GetParameter(0);
    double lslope = fitFunc->GetParameter(1);
    double range = 0.05;
    if (unit=="GeV") {
      loffset=0;
      lslope=1;
      range = 0.1;
    }

    //fill delta
    for (int ip(0);ip<gr->GetN();++ip){
      double x=0;
      double y=0;
      gr->GetPoint(ip,x,y);
      grDelta->SetPoint(ip,x,((y-loffset)/lslope-x)/x);
      double err = gr->GetErrorY(ip)/lslope*1./x;
      grDelta->SetPointError(ip,0,err);
      std::cout << "Calib " << ip << " Egen=" << x << " Erec=" << y << " delta=" << ((y-loffset)/lslope-x)/x << std::endl;
    }
    grDelta->SetTitle("");
    grDelta->SetMinimum(-1.*range);
    grDelta->SetMaximum(range);
    grDelta->GetXaxis()->SetLabelSize(0.15);
    grDelta->GetXaxis()->SetTitleSize(0.15);
    grDelta->GetYaxis()->SetLabelSize(0.12);
    grDelta->GetYaxis()->SetTitleSize(0.15);
    grDelta->GetXaxis()->SetTitleOffset(0.5);
    grDelta->GetYaxis()->SetTitleOffset(0.3);
    
    grDelta->Draw("ap");
    //grDelta->GetYaxis()->SetRangeUser(0,grDelta->GetYaxis()->GetXmax());
    if (!dovsE) grDelta->GetXaxis()->SetTitle("E_{T} (GeV)");
    else grDelta->GetXaxis()->SetTitle("E (GeV)");
    grDelta->GetYaxis()->SetTitle("(#Delta E)/E");

    TLine *line = new TLine(grDelta->GetXaxis()->GetXmin(),0,grDelta->GetXaxis()->GetXmax(),0);
    line->SetLineColor(2);//kYellow+4);
    line->Draw();
    

    lower->Update();

  }

  return upper;
};


int plotSmearingEffect(){//main

  SetTdrStyle();

  bool dovsE = true;

  const unsigned nIC = 1;
  const unsigned ICval[nIC] = {10};//,1,2,3,4,5,10,15,20,50};

  const unsigned nPts = 900;
  const unsigned startIdx = 100;

  const unsigned nPu = 1;//4;
  unsigned pu[nPu] = {0};//,140,200};

  const unsigned nS = 1;
  std::string scenario[nS] = {
    "gamma/200um/"
  };

  const unsigned neta = 1;
  unsigned eta[neta]={21};
  //const unsigned neta = 7;
  //unsigned eta[neta]={17,19,21,23,25,27,29};

  double etaval[neta];
  double etaerr[neta];

  const unsigned nEvtMin = 150;

  TString pSuffix = "";

  
  const unsigned nV = 1;
  TString version[nV] = {"12"};//,"0"};
  
  const unsigned nLayers = 30;
  const unsigned nSR = 1;

  double calib[nPu][neta][nSR];
  double calibErr[nPu][neta][nSR];
  double offset[nPu][neta][nSR];
  double offsetErr[nPu][neta][nSR];

  std::ostringstream saveName;

  //unsigned genEnAll[]={3,5,7,10,20,30,40,50,60,70,80,90,100,125,150,175,200};
  unsigned genEnAll[]={5,10,20,50,100};
  const unsigned nGenEnAll=sizeof(genEnAll)/sizeof(unsigned);

  double fitQual = 30;

  //canvas so they are created only once
  TCanvas *mycE[nGenEnAll];
  for (unsigned iE(0); iE<nGenEnAll;++iE){
    std::ostringstream lName;
    lName << "mycE" << genEnAll[iE];
    mycE[iE] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1500,1000);
    //mycE[iE]->Divide(4,2);
  }

  TCanvas *mycEtot = new TCanvas("mycEtot","mycEtot",1500,1000);
  mycEtot->Divide(2,3);
  const unsigned nCanvas = 1;  
  TCanvas *myc[nCanvas];
  for (unsigned iC(0);iC<nCanvas;++iC){
    std::ostringstream lName;
    lName << "myc" << iC;
    myc[iC] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1500,1000);
    //myc[iC]->Divide(2,4);
  }
  for (unsigned ip(startIdx);ip<(startIdx+nPts);++ip){//loop on tries
 
    for (unsigned ic(0);ic<nIC;++ic){//loop on intercalib
      
      TFile *fcalib;
      std::ostringstream label;
      label << "PLOTS/CalibTree";
      label << "_SR7_IC" << ICval[ic];
      label << "_try" << ip;
      label << ".root";
      fcalib = TFile::Open(label.str().c_str(),"RECREATE");
      fcalib->cd();
      TTree *outtree[nGenEnAll];
      std::vector<double> Etot;
      Etot.resize(nGenEnAll,0);
      std::vector<double> Esmear;
      Esmear.resize(nGenEnAll,0);
      double calibSlope = 0;
      double calibOffset = 0;
      for (unsigned iE(0); iE<nGenEnAll;++iE){
	label.str("");     
	label << "TreeEtot_" << genEnAll[iE];
	outtree[iE] = new TTree(label.str().c_str(),"Etot calibrated and smeared");
	label.str("");     
	label << "Etot";
	outtree[iE]->Branch(label.str().c_str(),&Etot[iE]);
	label.str("");     
	label << "Esmear";
	outtree[iE]->Branch(label.str().c_str(),&Esmear[iE]);
	outtree[iE]->Branch("calibSlope",&calibSlope);
	outtree[iE]->Branch("calibOffset",&calibOffset);
      }

      TRandom3 lrndm;
      lrndm.SetSeed(0);
      double smearFact[nLayers];
      for (unsigned iL(0); iL<nLayers;++iL){
	smearFact[iL]=0;
	while (1){
	  smearFact[iL] = lrndm.Gaus(1.,ICval[ic]/100.);
	  if (smearFact[iL]>0) break;
	}
      }

      for (unsigned iV(0); iV<nV;++iV){//loop on versions
	for (unsigned iS(0); iS<nS;++iS){//loop on scenarios
	  
	  TString plotDir = "/afs/cern.ch/work/a/amagnan/PFCalEEAna/PLOTS/gitV00-02-12/version"+version[iV]+"/"+scenario[iS]+"/";
	  TTree *ltree[neta][nPu][nGenEnAll];
	  TFile *inputFile[neta][nPu][nGenEnAll];
	  
	  for (unsigned ieta(0);ieta<neta;++ieta){//loop on eta
	    
	    etaval[ieta] = eta[ieta]/10.;
	    etaerr[ieta] = 0;
	    
	    for (unsigned ipu(0); ipu<nPu; ++ipu){//loop on pu
	      unsigned puOption = pu[ipu];
	      
	      std::string unit = "MIPS";
	      if (ipu>0) unit = "GeV";
	      
	      //identify valid energy values
	      bool skip[nGenEnAll];
	      unsigned nValid = 0;
	      for (unsigned iE(0); iE<nGenEnAll; ++iE){
		ltree[ieta][ipu][iE] = 0;
		skip[iE] = false;
		inputFile[ieta][ipu][iE] = 0;
		std::ostringstream linputStr;
		if (pu[ipu]!=200) linputStr << plotDir ;
		else linputStr << "/afs/cern.ch/work/a/amagnan/PFCalEEAna/PLOTS/gitV00-02-12/version"+version[iV]+"/"+scenario[iS]+"/";
		linputStr << "eta" << eta[ieta] << "_et" << genEnAll[iE] << "_pu" << pu[ipu] << "_IC0";// << ICval[ic];
		linputStr << ".root";
		inputFile[ieta][ipu][iE] = TFile::Open(linputStr.str().c_str());
		if (!inputFile[ieta][ipu][iE]) {
		  std::cout << " -- Error, input file " << linputStr.str() << " cannot be opened. Skipping..." << std::endl;
		  //	    return 1;
		  skip[iE] = true;
		}
		else {
		  inputFile[ieta][ipu][iE]->cd("Energies");
		  ltree[ieta][ipu][iE] = (TTree*)gDirectory->Get("Ereso");
		  
		  if (!ltree[ieta][ipu][iE]){
		    std::cout << " -- File " << inputFile[ieta][ipu][iE]->GetName() << " sucessfully opened but tree Ereso not found! Skipping." << std::endl;
		    skip[iE] = true;
		  } else { 
		    std::cout << " -- File " << inputFile[ieta][ipu][iE]->GetName() << " sucessfully opened and tree found." << std::endl;
		    nValid++;
		  }
		}
	      }
	      
	      unsigned newidx = 0;
	      const unsigned nGenEn = nValid;
	      unsigned genEn[nGenEn];
	      unsigned oldIdx[nGenEn];
	      for (unsigned iE(0); iE<nGenEnAll; ++iE){	  
		if (!skip[iE]) {
		  genEn[newidx]=genEnAll[iE];
		  oldIdx[newidx]=iE;
		  newidx++;
		}
	      }
	      

	      TH1F *p_Ereco[nGenEn][nSR];
	      TGraphErrors *calibRecoFit[nSR];
	      TGraphErrors *calibRecoDelta[nSR];
	      unsigned iSR = 7;
	      TString srStr = "";
	      srStr += iSR;
	      srStr += "eta";
	      srStr += eta[ieta];
	      srStr += "pu";
	      srStr += ipu;
	      calibRecoFit[0] = new TGraphErrors();
	      calibRecoFit[0]->SetName("calibRecoFit"+srStr);
	      calibRecoFit[0]->SetTitle("");
	      calibRecoFit[0]->SetMarkerStyle(20);
	      calibRecoFit[0]->SetMarkerColor(1);
	      calibRecoFit[0]->SetLineColor(1);
	      calibRecoDelta[0] = (TGraphErrors *) calibRecoFit[0]->Clone("calibRecoDelta"+srStr);

	      for (unsigned iE(0); iE<nGenEn; ++iE){
		
		std::cout << "- Processing energy : " << genEn[iE] 
			  << std::endl;
		
		inputFile[ieta][ipu][oldIdx[iE]]->cd("Energies");
		
		std::cout << " -- Tree entries for eta=" << eta[ieta] << " pu=" << pu[ipu] << " : " << ltree[ieta][ipu][oldIdx[iE]]->GetEntries() << std::endl;
		
		
		gStyle->SetOptStat(0);
		gStyle->SetOptFit(0);

		mycE[iE]->cd();

		std::ostringstream lName;
		lName.str("");
		lName << std::setprecision(6);
		//lName << "wgtEtotal";///" << tanh(etaval[ieta]);
		for (unsigned iL(0);iL<nLayers;++iL){
		  double fact = absWeight(iL,etaval[ieta]);
		  if (iL==0) lName << fact << "*subtractedenergy_" << iL << "_SR4";
		  else lName << "+" << fact << "*subtractedenergy_" << iL << "_SR4";
		}

		ltree[ieta][ipu][oldIdx[iE]]->Draw(lName.str().c_str(),"","");
		lName.str("");
		lName << "energy" << genEn[iE] << "_SR" << iSR ;
		p_Ereco[iE][0] = (TH1F*)(gPad->GetPrimitive("htemp"))->Clone(lName.str().c_str()); // 1D
		//if (iSR==7) p_Ereco[iE][0] = (TH1F*)gDirectory->Get("p_wgtEtotal");
		if (!p_Ereco[iE][0]){
		  std::cout << " -- ERROR, pointer for histogram " << lName.str() << " is null." << std::endl;
		  return 1;
		}
		std::cout << " --- Reco E = entries " << p_Ereco[iE][0]->GetEntries() 
			  << " mean " << p_Ereco[iE][0]->GetMean() 
			  << " rms " << p_Ereco[iE][0]->GetRMS() 
			  << " overflows " << p_Ereco[iE][0]->GetBinContent(p_Ereco[iE][0]->GetNbinsX()+1)
			  << std::endl;
		
		
		p_Ereco[iE][0]->SetTitle((";E ("+unit+");events").c_str());
		
		//take min 20 bins
		//if(p_Ereco[iE][0]->GetNbinsX()>40) p_Ereco[iE][0]->Rebin(2);
		
		//skip data with too little stat
		if (p_Ereco[iE][0]->GetEntries()<nEvtMin) {
		  gPad->Clear();
		  continue;
		}
		
		TPad *lpad = (TPad*)(mycE[iE]->cd());
		FitResult lres;
		if (fitEnergy(p_Ereco[iE][0],lpad,unit,lres,iSR)!=0) return 1;
		lpad->cd();
		char buf[500];
		sprintf(buf,"#gamma E_{T}=%d GeV + PU %d",genEn[iE],pu[ipu]);
		TLatex lat;
		lat.SetTextSize(0.05);
		lat.DrawLatexNDC(0.25,0.965,buf);
		sprintf(buf,"#eta=%3.1f, SR %d",etaval[ieta],iSR);
		lat.SetTextSize(0.06);
		lat.DrawLatexNDC(0.15,0.87,buf);
		lat.DrawLatexNDC(0.01,0.01,"HGCAL G4 standalone");
		
		if (ipu>0){
		  mycEtot->cd(iE+1);
		  p_Ereco[iE][0]->SetLineColor(ipu+1);
		  p_Ereco[iE][0]->Draw(ipu==1?"":"same");
		  TF1 *fite = p_Ereco[iE][0]->GetFunction("fitResult");
		  fite->SetLineColor(ipu+1);
		  fite->Draw("same");
		  sprintf(buf,"#gamma E_{T}=%d GeV + PU %d",genEn[iE],pu[ipu]);
		  lat.SetTextSize(0.05);
		  lat.DrawLatexNDC(0.25,0.965,buf);
		  sprintf(buf,"#eta=%3.1f, IC %d %%",etaval[ieta],ICval[ic]);
		  lat.SetTextSize(0.06);
		  lat.DrawLatexNDC(0.15,0.87,buf);
		  lat.SetTextColor(ipu+1);
		  sprintf(buf,"#sigma/E = %3.1f/%3.1f=%3.3f",fite->GetParameter(2),fite->GetParameter(1),fite->GetParameter(2)/fite->GetParameter(1));
		  lat.DrawLatexNDC(0.15,0.85-0.1*ipu,buf);
		  lat.SetTextColor(1);
		  
		}
		
		
		//filter out bad points
		if (lres.chi2/lres.ndf > fitQual) {
		  std::cout << " --- INFO! Point Egen=" 
			    << genEn[iE] 
			    << " eta=" << etaval[ieta]
			    << " pu=" << pu[ipu]
			    << " skipped, chi2/ndf = "
			    << lres.chi2/lres.ndf
			    << std::endl;
		  continue;
		}

		Int_t np=calibRecoFit[0]->GetN();
		if (!dovsE) calibRecoFit[0]->SetPoint(np,genEn[iE],lres.mean);
		else calibRecoFit[0]->SetPoint(np,E(genEn[iE],eta[ieta]),lres.mean);
		calibRecoFit[0]->SetPointError(np,0.0,lres.meanerr);
		
		
		saveName.str("");
		saveName << plotDir << "/Ereco_eta" << eta[ieta] << "_pu" << puOption;
		if (ipu==0) saveName << "raw";
		saveName << "_E" << genEn[iE] << "_IC" << ICval[ic];
		mycE[iE]->Update();
		mycE[iE]->Print((saveName.str().c_str()+pSuffix)+".pdf");
		
		
	      }//loop on energies


	      //plot and fit calib
	      TPad *lpad = (TPad*)(myc[0]->cd());
	      TPad *upper = plotCalibration(calibRecoFit[0],lpad,
					    true,calibRecoDelta[0],
					    unit,
					    calib[ipu][ieta][0],
					    calibErr[ipu][ieta][0],
					    offset[ipu][ieta][0],
					    offsetErr[ipu][ieta][0],
					    eta[ieta],dovsE);
	      upper->cd();
	      char buf[500];
	      sprintf(buf,"#gamma #eta=%3.1f + PU %d",etaval[ieta],pu[ipu]);
	      TLatex lat;
	      lat.SetTextSize(0.1);
	      lat.DrawLatexNDC(0.7,0.15,buf);
	      sprintf(buf,"SR %d",iSR);
	      lat.DrawLatexNDC(0.7,0.3,buf);
	      lat.SetTextSize(0.06);
	      lat.DrawLatexNDC(0.01,0.94,"HGCAL G4 standalone");
	      
	      
	      myc[0]->Update();
	      std::ostringstream lsave;
	      lsave.str("");
	      lsave << plotDir << "/";
	      if (ipu==0) lsave << "CalibMipToGeV";
	      else lsave << "Calib";
	      lsave << "_eta" << eta[ieta] << "_pu" << puOption << pSuffix << "_IC" << ICval[ic];
	      if (dovsE) lsave << "_vsE";
	      myc[0]->Print((lsave.str()+".pdf").c_str());
	      
	      calibOffset = offset[0][ieta][0];
	      calibSlope = calib[0][ieta][0];


	      for (unsigned iE(0); iE<nGenEn; ++iE){
		
		std::cout << "- Processing energy : " << genEn[iE] 
			  << std::endl;
		
		inputFile[ieta][ipu][oldIdx[iE]]->cd("Energies");
		const unsigned nEvts = ltree[ieta][ipu][oldIdx[iE]]->GetEntries();
		//loop over tree entries
		std::vector<std::vector<double> > energySR;
		std::vector<std::vector<double> > subtractedenergySR;
		std::vector<double> emptyvec;
		emptyvec.resize(5,0);
		energySR.resize(nLayers,emptyvec);
		subtractedenergySR.resize(nLayers,emptyvec);
		for (unsigned iL(0); iL<nLayers;++iL){
		  for (unsigned i(0);i<5;++i){
		    label.str("");
		    label << "energy_" << iL << "_SR" << i;
		    ltree[ieta][ipu][oldIdx[iE]]->SetBranchAddress(label.str().c_str(),&energySR[iL][i]);
		    label.str("");
		    label << "subtractedenergy_" << iL << "_SR" << i;
		    ltree[ieta][ipu][oldIdx[iE]]->SetBranchAddress(label.str().c_str(),&subtractedenergySR[iL][i]);
		  }
		}

		for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
		  
		  ltree[ieta][ipu][oldIdx[iE]]->GetEntry(ievt);
		  
		  Etot[oldIdx[iE]] = 0;
		  Esmear[oldIdx[iE]]=0;
		  for (unsigned iL(0);iL<nLayers;++iL){
		    double fact = absWeight(iL,etaval[ieta]);
		    Etot[oldIdx[iE]] += fact*subtractedenergySR[iL][4];
		    Esmear[oldIdx[iE]] += fact*subtractedenergySR[iL][4]*smearFact[iL];
		  }
		  outtree[oldIdx[iE]]->Fill();
		}
		inputFile[ieta][ipu][oldIdx[iE]]->Close();
	      }//loop on energies

		
	    }//loop on pu


	    saveName.str("");
	    saveName << plotDir << "/Ereco_eta" << eta[ieta];
	    saveName << "_SR7" << "_IC" << ICval[ic] << "_try" << ip;
	    mycEtot->Update();
	    mycEtot->Print((saveName.str().c_str()+pSuffix)+".pdf");
	    
	    
	  }//loop on eta

	  
	}//loop on scenarios
	
      }//loop on versions
      fcalib->cd();

      for (unsigned iE(0); iE<nGenEnAll; ++iE){
	outtree[iE]->Write();
      }
      //fcalib->Write();
      fcalib->Close();
      
    }//loop on IC vals
    
  }//loop on points
  
  return 0;
  
}//main
