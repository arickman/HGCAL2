#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>

#include "TFile.h"
#include "TVirtualPad.h"
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



void extractMeanEnergy(int *color, int *marker,
		       TTree *tree,
		       unsigned nPts,
		       TString *Enames,unsigned genEn,
		       float *meanE,float *rmsE,
		       float *meanEerr,float *rmsEerr,
		       float *meanFitE,float *rmsFitE,
		       float *meanFitEerr,float *rmsFitEerr,
		       std::string plotDir,
		       double *offset,
		       double *slope,
		       double * a, double * b, double * c){
  TCanvas *mycE = new TCanvas("mycE","mycE",1500,1000);

  TH1F *hrefit[nPts];
  double lmax = 0;
  for (unsigned iP(0); iP<nPts; ++iP){
    //if (iP==0) 
    mycE->cd();
    TString varname;
    if (iP<2){
      varname = "("+Enames[iP]+"-";
      varname += offset[iP];
      varname += ")/";
      varname += slope[iP];
    }
    else {
      varname = Enames[iP];
      varname += "*(G4XT2d5Rand1156N3Noise12-";
      varname += offset[iP];
      varname += ")/";
      varname += slope[iP];
    }
    TString full = varname;
    full += "*(";
    full += a[iP];
    full += "+";
    full += b[iP];
    full += "*";
    full += varname;
    full += "+";
    full += c[iP];
    full += "*";
    full += varname;
    full += "*";
    full += varname;
    full += ")";
    tree->Draw(full,"","");

    std::ostringstream varstr;
    varstr << Enames[iP] << "_e" << genEn;
    hrefit[iP] = (TH1F*)(gPad->GetPrimitive("htemp"))->Clone(varstr.str().c_str()); // 1D
    //hrefit[iP]->Rebin(2);
    hrefit[iP]->SetLineColor(color[iP]);
    hrefit[iP]->SetMarkerColor(color[iP]);
    hrefit[iP]->SetMarkerStyle(marker[iP]);
    hrefit[iP]->Draw("same");
    hrefit[iP]->Sumw2();

    meanE[iP] = hrefit[iP]->GetMean();
    rmsE[iP] = hrefit[iP]->GetRMS();
    meanEerr[iP] = hrefit[iP]->GetMeanError();
    rmsEerr[iP] = hrefit[iP]->GetRMSError();
    
    std::cout << " --- " << varstr.str() << " = entries " << hrefit[iP]->GetEntries()
	      << " nbins " << hrefit[iP]->GetNbinsX()
	      << " " << hrefit[iP]->GetBinLowEdge(1)
	      << " " << hrefit[iP]->GetBinLowEdge(hrefit[iP]->GetNbinsX())
	      << " mean " << meanE[iP]
	      << " rms " << rmsE[iP]
	      << " underflows " << hrefit[iP]->GetBinContent(0)
	      << " overflows " << hrefit[iP]->GetBinContent(hrefit[iP]->GetNbinsX()+1) 
	      << std::endl;
    varstr.str("");
    varstr << "Egen = " << genEn << " GeV";
    hrefit[iP]->SetTitle(varstr.str().c_str());
    hrefit[iP]->GetXaxis()->SetTitle("Energy");
    if (hrefit[iP]->GetMaximum()>lmax) lmax = hrefit[iP]->GetMaximum();

  }
  //fit
  mycE->cd();
  for (unsigned iP(0); iP<nPts; ++iP){
    hrefit[iP]->SetMaximum(lmax);
    if (iP==0) hrefit[iP]->Draw("PE");
    else hrefit[iP]->Draw("PEsame");
    //fit
    hrefit[iP]->Fit("gaus","LR0","",
		    meanE[iP]-2*rmsE[iP],
		    meanE[iP]+2*rmsE[iP]);
    TF1 *fitResult = hrefit[iP]->GetFunction("gaus");
    hrefit[iP]->Fit("gaus","LR+","same",
		    fitResult->GetParameter(1)-2*fitResult->GetParameter(2),
		    fitResult->GetParameter(1)+2*fitResult->GetParameter(2));
    fitResult = hrefit[iP]->GetFunction("gaus");
    fitResult->SetLineColor(color[iP]);
    
    meanFitE[iP] = fitResult->GetParameter(1);
    rmsFitE[iP] = fitResult->GetParameter(2);                                                                                  
    meanFitEerr[iP] = fitResult->GetParError(1);
    rmsFitEerr[iP] = fitResult->GetParError(2);                                                                                  
    
    fitResult->Draw("same");
    
  }//loop on pts
  std::ostringstream saveName;
  saveName.str("");
  saveName << plotDir << "/FitEtotal_e" << genEn;
  mycE->Update();
  mycE->Print((saveName.str()+".png").c_str());
  mycE->Print((saveName.str()+".pdf").c_str());

};//method


TPad* plot_ratio(TCanvas *canv, bool up){
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
    pad = new TPad("lower","pad",0, 0,1 ,0.26);  
    pad->SetTopMargin(0.05);
    pad->SetBottomMargin(0.24);
    pad->Draw();
    pad->cd();
    return pad;
  }

};

int plotEfromTree(){//main

  std::string baseDir = "../PLOTS/gitV00-02-04/version23/e-/";
  std::string baseDirPi = "../PLOTS/gitV00-02-04/version23/pi-/";

  bool isEM = baseDir.find("e-")!=baseDir.npos;

  TString pSuffix = "";

  const bool doVsE = true;

  unsigned genEn[]={10,15,18,20,25,30,35,40,45,50,60,80};
  const unsigned nGenEn=sizeof(genEn)/sizeof(unsigned);

  // unsigned nx=0,ny=0;

  // if (nGenEn>12) {nx=5;ny=3;}
  // else if (nGenEn > 10)
  //   {nx=4;ny=3;}
  // else if (nGenEn > 6)
  //   {nx=5;ny=2;}
  // else if (nGenEn > 4)
  //   {nx=3;ny=2;}
  // else if (nGenEn > 2)
  //   {nx=2;ny=2;}
  // else 
  //   {nx=nGenEn;ny=1;}
  
  // mycE->Divide(nx,ny);
  // std::cout << " Divide: " << nx <<  " " << ny << std::endl;
  const unsigned nLimits = 10;//5;
  const double pElim[nLimits] = {3,4,5,6,7,8,9,10,11,12};
  const unsigned idxRef = 2;

  const unsigned nEPts = 2;//14;
  TString EnamesE[nEPts] = {
  "G4",
  //"G4mipcut",
  //"G4Noise12",
  //"G4Noise15",
  //"G4Noise20",
  //"G4XT2d5Noise",
  //"G4XT3d5Noise",
  //"G4XT5Noise",
  //"G4Rand1156N3Noise",
  //"G4Rand1156N6Noise",
  //"G4Rand925N3Noise",
  //"G4Rand925N6Noise",
  "G4XT2d5Rand1156N3Noise12",
  //"G4XT3d5Rand925N6Noise15"
  };
  TString EnamesC[nLimits];

  for (unsigned ilim(0);ilim<nLimits;++ilim){
    std::ostringstream bname;
    bname << "Cglobal_" << static_cast<unsigned>(pElim[ilim]) << "mip";
    EnamesC[ilim] = bname.str();
  }
  unsigned nPts = nEPts+nLimits;

  TString Enames[nPts];
  for (unsigned i(0);i<nPts;++i){
    if (i<nEPts) Enames[i] = EnamesE[i];
    else Enames[i] = EnamesC[i-nEPts];
  }


  int marker[14] = {20,28,21,22,23,21,22,23,24,24,25,25,26,27};
  int color[14] = {1,1,kRed-2,kRed,kRed+2,kGreen-2,kGreen,kGreen+2,kCyan,kBlue,kYellow,kYellow+2,kViolet-1,kViolet+1};

  const unsigned nCanvas = nPts+6;
  TCanvas *myc[nCanvas];
  TCanvas *mycinv[nPts];
  for (unsigned iC(0);iC<nCanvas;++iC){
    std::ostringstream lName;
    lName << "myc" << iC;
    myc[iC] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1);
    lName.str("");
    if (iC<nPts){
      lName << "mycinv" << iC;
      mycinv[iC] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1);
    }
  }

  double offsetVal[nPts];
  double slopeVal[nPts];
  double aVal[nPts];
  double bVal[nPts];
  double cVal[nPts];
  for (unsigned iP(0); iP<nPts; ++iP){
    offsetVal[iP] = 0;
    slopeVal[iP] = 1;
    aVal[iP] = 1;
    bVal[iP] = 0;
    cVal[iP] = 0;
  }

  for (unsigned iC(0); iC<4;++iC){//calib

    TLegend *leg = new TLegend(0,0,1,1);
    leg->SetFillColor(10);

    float meanE[nGenEn][nPts];
    float rmsE[nGenEn][nPts];
    float meanEerr[nGenEn][nPts];
    float rmsEerr[nGenEn][nPts];
    float meanFitE[nGenEn][nPts];
    float rmsFitE[nGenEn][nPts];
    float meanFitEerr[nGenEn][nPts];
    float rmsFitEerr[nGenEn][nPts];
    
    gStyle->SetOptStat(0);
    std::string plotDir;
    if (iC==1 || iC==3) plotDir = baseDir+"calib/";
    else if (iC==2) {
      baseDir = baseDirPi;
      isEM = false;
      plotDir = baseDir;
      for (unsigned iP(0); iP<nEPts; ++iP){
	slopeVal[iP] = slopeVal[iP]/1.19;
      }
    }
    else plotDir = baseDir;
    
    for (unsigned iE(0); iE<nGenEn; ++iE){
      std::cout << "- Processing energy : " << genEn[iE] 
		<< std::endl;
      
      TFile *inputFile = 0;
      std::ostringstream linputStr;
      linputStr << baseDir << "validateCalice_e" << genEn[iE] << pSuffix << ".root";
      inputFile = TFile::Open(linputStr.str().c_str());
      if (!inputFile) {
	std::cout << " -- Error, input file " << linputStr.str() << " cannot be opened. Exiting..." << std::endl;
	return 1;
      }
      else std::cout << " -- File " << inputFile->GetName() << " successfully opened." << std::endl;
      
      TTree *tree = (TTree*)gDirectory->Get("Estudy");
      if (!tree) {
	std::cout << " Tree not found." << std::endl;
	return 1;
      }
      else std::cout << " tree found." << std::endl;
      
      //extract histos
      extractMeanEnergy(color,marker,tree,nPts,
			Enames,genEn[iE],
			meanE[iE],rmsE[iE],
			meanEerr[iE],rmsEerr[iE],
			meanFitE[iE],rmsFitE[iE],
			meanFitEerr[iE],rmsFitEerr[iE],
			plotDir,
			offsetVal,slopeVal,
			aVal,bVal,cVal);
      
      inputFile->Close();
      
    }//loop on energies
    
    std::ostringstream saveName;
    saveName.str("");
    
    TGraphErrors *calib[nPts];
    TGraphErrors *sigma[nPts];
    TGraphErrors *reso[nPts];
    TGraphErrors *calibinv[nPts];
    TGraphErrors *calibFit[nPts];
    TGraphErrors *deltaFit[nPts];
    TGraphErrors *sigmaFit[nPts];
    TGraphErrors *resoFit[nPts];
    
    TGraphErrors * gr[nPts];
    TGraphErrors * grDelta[nPts];
    
    TGraphErrors *groffset = new TGraphErrors();
    TGraphErrors *grslope = new TGraphErrors();

    //draw all in one plot
    TPad *upperTot = plot_ratio(myc[nPts], true);
    TPad *lowerTot = plot_ratio(myc[nPts], false);
    if (!upperTot || !lowerTot){
      std::cout << " Pb..." << upperTot << " " << lowerTot << std::endl;
      return 1;
    }

    //extract linearity
    for (unsigned iP(0); iP<nPts; ++iP){
      TString lSuf = "_";
      lSuf += iP;
      //draw calibration curves
      calib[iP] = new TGraphErrors();
      calib[iP]->SetName("calib"+lSuf);
      calib[iP]->SetMarkerStyle(marker[iP]);
      calib[iP]->SetLineColor(color[iP]);
      calib[iP]->SetLineStyle(3);
      calib[iP]->SetMarkerColor(color[iP]);
      calib[iP]->SetTitle("");
      sigma[iP] = (TGraphErrors *) calib[iP]->Clone("sigma"+lSuf);
      reso[iP] = (TGraphErrors *) calib[iP]->Clone("reso"+lSuf);
      calibFit[iP] = (TGraphErrors *) calib[iP]->Clone("calibFit"+lSuf);
      calibFit[iP]->SetMarkerStyle(marker[iP]);
      calibFit[iP]->SetLineStyle(1);
      calibFit[iP]->SetLineColor(color[iP]);
      calibFit[iP]->SetMarkerColor(color[iP]);
      calibinv[iP] = (TGraphErrors *) calibFit[iP]->Clone("calibinv"+lSuf);
      deltaFit[iP] = (TGraphErrors *) calibFit[iP]->Clone("deltaFit"+lSuf);
      sigmaFit[iP] = (TGraphErrors *) calibFit[iP]->Clone("sigmaFit"+lSuf);
      
      resoFit[iP] = (TGraphErrors *) calibFit[iP]->Clone("resoFit"+lSuf);
      
      for (unsigned iE(0); iE<nGenEn; ++iE){
	Int_t np=calib[iP]->GetN();
	calib[iP]->SetPoint(np,genEn[iE],meanE[iE][iP]);
	calib[iP]->SetPointError(np,0.0,meanEerr[iE][iP]);
	sigma[iP]->SetPoint(np,genEn[iE],rmsE[iE][iP]);
	sigma[iP]->SetPointError(np,0.0,rmsEerr[iE][iP]);
	reso[iP]->SetPoint(np,doVsE?genEn[iE] : 1/sqrt(genEn[iE]),rmsE[iE][iP]/meanE[iE][iP]);
	double err = rmsE[iE][iP]/meanE[iE][iP]*sqrt(pow(rmsEerr[iE][iP]/rmsE[iE][iP],2)+pow(meanEerr[iE][iP]/meanE[iE][iP],2));
	reso[iP]->SetPointError(np,0,err);
	
	calibinv[iP]->SetPoint(np,meanFitE[iE][iP],genEn[iE]);
	calibinv[iP]->SetPointError(np,meanFitEerr[iE][iP],0.01*genEn[iE]);
	calibFit[iP]->SetPoint(np,genEn[iE],meanFitE[iE][iP]);
	calibFit[iP]->SetPointError(np,0.0,meanFitEerr[iE][iP]);
	sigmaFit[iP]->SetPoint(np,genEn[iE],rmsFitE[iE][iP]);
	sigmaFit[iP]->SetPointError(np,0.0,rmsFitEerr[iE][iP]);
	resoFit[iP]->SetPoint(np,doVsE?genEn[iE] : 1/sqrt(genEn[iE]),rmsFitE[iE][iP]/meanFitE[iE][iP]);
	double errFit = rmsFitE[iE][iP]/meanFitE[iE][iP]*sqrt(pow(rmsFitEerr[iE][iP]/rmsFitE[iE][iP],2)+pow(meanFitEerr[iE][iP]/meanFitE[iE][iP],2));
	resoFit[iP]->SetPointError(np,0,errFit);
	
      }//loop on energies
      
      
      TPad *upper = plot_ratio(myc[iP], true);
      TPad *lower = plot_ratio(myc[iP], false);
      if (!upper || !lower){
	std::cout << " Pb..." << upper << " " << lower << std::endl;
	return 1;
      }

      //draw inverted calib
      mycinv[iP]->cd();
      gPad->SetGridx(1);
      gPad->SetGridy(1);
      gr[iP] = calibinv[iP];
      
      //gr[iP]->GetXaxis()->SetLabelSize(0.06);
      //gr[iP]->GetXaxis()->SetTitleSize(0.06);
      //gr[iP]->GetYaxis()->SetLabelSize(0.06);
      //gr[iP]->GetYaxis()->SetTitleSize(0.06);
      //gr[iP]->GetXaxis()->SetTitleOffset(0.7);
      //gr[iP]->GetYaxis()->SetTitleOffset(0.8);
      
      gr[iP]->SetTitle(Enames[iP]);
      gr[iP]->Draw("ap");
      gr[iP]->GetYaxis()->SetTitle("Beam energy [GeV]");
      if (iC==0) gr[iP]->GetXaxis()->SetTitle("Corrected energy deposited [MIPs]");
      else gr[iP]->GetXaxis()->SetTitle("Corrected energy deposited [GeV]");
      //gr[iP]->GetXaxis()->SetTitle(doVsE?"Beam energy [GeV]" :"1/#sqrt{Beam energy} [1/#sqrt{GeV}]"); 
      
      char buf[500];
      //TF1 *fitFunc=new TF1("calib","[0]+[1]*x",gr->GetXaxis()->GetXmin(),gr->GetXaxis()->GetXmax());
      TF1 *fitFuncinv=new TF1("calibinv","[0]*x+[1]*x*x+[2]*x*x*x",gr[iP]->GetXaxis()->GetXmin(),gr[iP]->GetXaxis()->GetXmax());
      fitFuncinv->SetLineColor(1);
      gr[iP]->Fit(fitFuncinv,"RME");
      gr[iP]->GetFunction("calibinv")->SetLineColor(color[iP]);
      TLatex lat;
      lat.SetTextColor(1);
      sprintf(buf,"<E> #propto E #times (a + b #times E +c #times E^{2})");
      lat.DrawLatex(genEn[0],gr[iP]->GetYaxis()->GetXmax()*0.9,buf);
      sprintf(buf,"a = %3.3f #pm %3.3f %s",fitFuncinv->GetParameter(0),fitFuncinv->GetParError(0),iC==0?"MIPs":"GeV");
      lat.DrawLatex(genEn[0],gr[iP]->GetYaxis()->GetXmax()*(0.8),buf);
      sprintf(buf,"b = %3.3f #pm %3.3f %s/GeV",fitFuncinv->GetParameter(1),fitFuncinv->GetParError(1),iC==0?"MIPs":"GeV");
      lat.DrawLatex(genEn[0],gr[iP]->GetYaxis()->GetXmax()*(0.7),buf);
      sprintf(buf,"c = %3.6f #pm %3.6f %s/GeV",fitFuncinv->GetParameter(2),fitFuncinv->GetParError(2),iC==0?"MIPs":"GeV");
      lat.DrawLatex(genEn[0],gr[iP]->GetYaxis()->GetXmax()*(0.6),buf);
      sprintf(buf,"chi2/NDF = %3.3f/%d = %3.3f",fitFuncinv->GetChisquare(),fitFuncinv->GetNDF(),fitFuncinv->GetChisquare()/fitFuncinv->GetNDF());
      lat.DrawLatex(genEn[0],gr[iP]->GetYaxis()->GetXmax()*(0.5),buf);
      
      if (iP>=nEPts && iC==2) {
	aVal[iP] = fitFuncinv->GetParameter(0);
	bVal[iP] = fitFuncinv->GetParameter(1);
	cVal[iP] = fitFuncinv->GetParameter(2);
      }

      mycinv[iP]->Update();
      saveName.str("");
      saveName << plotDir << "/CalibEshower_" << Enames[iP];
      mycinv[iP]->Update();
      mycinv[iP]->Print((saveName.str()+".png").c_str());
      mycinv[iP]->Print((saveName.str()+".pdf").c_str());
     
      //draw calib
      upper->cd();
      gPad->SetGridx(1);
      gPad->SetGridy(1);
      gr[iP] = calibFit[iP];
      
      gr[iP]->GetXaxis()->SetLabelSize(0.06);
      gr[iP]->GetXaxis()->SetTitleSize(0.06);
      gr[iP]->GetYaxis()->SetLabelSize(0.06);
      gr[iP]->GetYaxis()->SetTitleSize(0.06);
      gr[iP]->GetXaxis()->SetTitleOffset(0.7);
      gr[iP]->GetYaxis()->SetTitleOffset(0.8);
      
      gr[iP]->SetTitle(Enames[iP]);
      gr[iP]->SetMaximum(std::max(calib[iP]->GetMaximum(),calibFit[iP]->GetMaximum()));
      gr[iP]->Draw("ap");
      calib[iP]->Draw("pl");
      gr[iP]->GetXaxis()->SetTitle("");
      if (iC==0) gr[iP]->GetYaxis()->SetTitle("Average energy deposited [MIPs]");
      else gr[iP]->GetYaxis()->SetTitle("Average energy deposited [GeV]");
      //gr[iP]->GetXaxis()->SetTitle(doVsE?"Beam energy [GeV]" :"1/#sqrt{Beam energy} [1/#sqrt{GeV}]"); 
      
      //TF1 *fitFunc=new TF1("calib","[0]+[1]*x",gr->GetXaxis()->GetXmin(),gr->GetXaxis()->GetXmax());
      TF1 *fitFunc=new TF1("calib","[0]+[1]*x",9,51);
      fitFunc->SetLineColor(1);
      gr[iP]->Fit(fitFunc,"RME");
      gr[iP]->GetFunction("calib")->SetLineColor(color[iP]);
      lat.SetTextColor(1);
      sprintf(buf,"<E> #propto a + b #times E ");
      lat.DrawLatex(genEn[0],gr[iP]->GetYaxis()->GetXmax()*0.9,buf);
      sprintf(buf,"a = %3.3f #pm %3.3f %s",fitFunc->GetParameter(0),fitFunc->GetParError(0),iC==0?"MIPs":"GeV");
      lat.DrawLatex(genEn[0],gr[iP]->GetYaxis()->GetXmax()*(0.8),buf);
      sprintf(buf,"b = %3.3f #pm %3.3f %s/GeV",fitFunc->GetParameter(1),fitFunc->GetParError(1),iC==0?"MIPs":"GeV");
      lat.DrawLatex(genEn[0],gr[iP]->GetYaxis()->GetXmax()*(0.7),buf);
      sprintf(buf,"chi2/NDF = %3.3f/%d = %3.3f",fitFunc->GetChisquare(),fitFunc->GetNDF(),fitFunc->GetChisquare()/fitFunc->GetNDF());
      lat.DrawLatex(genEn[0],gr[iP]->GetYaxis()->GetXmax()*(0.6),buf);
      
      groffset->SetPoint(iP,iP,fitFunc->GetParameter(0));
      groffset->SetPointError(iP,0.0,fitFunc->GetParError(0));
      if (iC==0) offsetVal[iP] = fitFunc->GetParameter(0);
      if (iP>=nEPts && iC==0) offsetVal[iP] = offsetVal[1];
      grslope->SetPoint(iP,iP,fitFunc->GetParameter(1));
      grslope->SetPointError(iP,0.0,fitFunc->GetParError(1));
      if (iC==0) slopeVal[iP] = fitFunc->GetParameter(1);
      if (iP>=nEPts && iC==0) slopeVal[iP] = slopeVal[1];
      //groffset->GetXaxis()->SetBinLabel(iP+1,Enames[iP].Data());
      //grslope->GetHistogram()->GetXaxis()->SetBinLabel(iP+1,Enames[iP].Data());
      //draw deltaE/E vs E
      lower->cd();
      gPad->SetLogx(0);
      gPad->SetGridx(1);
      gPad->SetGridy(1);
      for (unsigned iE(0);iE<nGenEn;++iE){
	if (iC<1){
	  deltaFit[iP]->SetPoint(iE,genEn[iE],( ((meanFitE[iE][iP]-fitFunc->GetParameter(0))/fitFunc->GetParameter(1))-genEn[iE])/genEn[iE]);
	  deltaFit[iP]->SetPointError(iE,0.0,meanFitEerr[iE][iP]/fitFunc->GetParameter(1)*1./genEn[iE]);
	}
	else {
	  deltaFit[iP]->SetPoint(iE,genEn[iE],(meanFitE[iE][iP]-genEn[iE])/genEn[iE]);
	  deltaFit[iP]->SetPointError(iE,0.0,meanFitEerr[iE][iP]/genEn[iE]);
	}
      }
      
      grDelta[iP] = deltaFit[iP];
      grDelta[iP]->SetTitle("");
      if (iC<2){
	grDelta[iP]->SetMinimum(-0.03);
	grDelta[iP]->SetMaximum(0.03);
      }
      else {
	grDelta[iP]->SetMinimum(-0.1);
	grDelta[iP]->SetMaximum(0.1);
      }
      grDelta[iP]->GetXaxis()->SetLabelSize(0.15);
      grDelta[iP]->GetXaxis()->SetTitleSize(0.15);
      grDelta[iP]->GetYaxis()->SetLabelSize(0.12);
      grDelta[iP]->GetYaxis()->SetTitleSize(0.15);
      grDelta[iP]->GetXaxis()->SetTitleOffset(0.5);
      grDelta[iP]->GetYaxis()->SetTitleOffset(0.3);
      
      grDelta[iP]->Draw("ap");
      //grDelta[iP]->GetYaxis()->SetRangeUser(0,grDelta->GetYaxis()->GetXmax());
      grDelta[iP]->GetXaxis()->SetTitle("Beam energy [GeV]");
      grDelta[iP]->GetYaxis()->SetTitle("(#Delta E)/E");
      
      myc[iP]->Update();
      saveName.str("");
      saveName << plotDir << "/CalibE_" << Enames[iP];
      myc[iP]->Update();
      myc[iP]->Print((saveName.str()+".png").c_str());
      myc[iP]->Print((saveName.str()+".pdf").c_str());
     
      upperTot->cd();
      gPad->SetGridx(1);
      gPad->SetGridy(1);
      if (iP==0) gr[iP]->Draw("ap");
      else gr[iP]->Draw("psame");
      lowerTot->cd();
      gPad->SetGridx(1);
      gPad->SetGridy(1);
      if (iP==0) grDelta[iP]->Draw("ap");
      else grDelta[iP]->Draw("psame");
    
      myc[nPts+1]->cd();
      gPad->SetGridx(1);
      gPad->SetGridy(1);
      sigmaFit[iP]->SetTitle("");
      sigmaFit[iP]->SetMaximum(std::max(sigma[iP]->GetMaximum(),sigmaFit[iP]->GetMaximum()));
      if (iP==0) sigmaFit[iP]->Draw("ap");
      else sigmaFit[iP]->Draw("p");
      //sigma[iP]->Draw("pl");
      sigmaFit[iP]->GetXaxis()->SetTitle("Beam energy [GeV]");
      if (iC==0) sigmaFit[iP]->GetYaxis()->SetTitle("RMS energy deposited [MIPs]");
      else sigmaFit[iP]->GetYaxis()->SetTitle("RMS energy deposited [GeV]");
      leg->AddEntry(sigmaFit[iP],Enames[iP],"P");
      
      myc[nPts+2]->cd();
      gPad->SetGridx(1);
      gPad->SetGridy(1);
      resoFit[iP]->SetTitle("");
      //resoFit[iP]->SetMaximum(std::max(reso[iP]->GetMaximum(),resoFit[iP]->GetMaximum()));
      if (iP==0) resoFit[iP]->Draw("ap");
      else resoFit[iP]->Draw("p");
      //reso[iP]->Draw("pl");
      resoFit[iP]->GetXaxis()->SetTitle("Beam energy (GeV)");
      resoFit[iP]->GetYaxis()->SetTitle("Relative energy resolution");
      
      TF1* fitref =new TF1("reso","sqrt([0]/x+[1]+[2]/(x*x))",resoFit[iP]->GetXaxis()->GetXmin(),resoFit[iP]->GetXaxis()->GetXmax());
      if (isEM) fitref->SetParameters(0.215*0.215,0.007*0.007,0.06*0.06);
      else {
	if (iP<nEPts) fitref->SetParameters(0.518*0.518,0.04*0.04,0.18*0.18);
	else fitref->SetParameters(0.436*0.436,0.0*0.0,0.18*0.18);
      }
      fitref->SetLineColor(4);
      fitref->SetLineStyle(2);
      fitref->Draw("same");
      lat.SetTextColor(4);
      //sprintf(buf,"CALICE");
      sprintf(buf,"CALICE s=%3.3f, c=%3.3f, n=%3.3f",sqrt(fitref->GetParameter(0)),sqrt(fitref->GetParameter(1)),sqrt(fitref->GetParameter(2)));
      if (iP==0) lat.DrawLatex(30,resoFit[iP]->GetYaxis()->GetXmax()*0.9,buf);


      myc[iP]->Clear();
      myc[iP]->cd();
      gPad->SetGridx(1);
      gPad->SetGridy(1);
      resoFit[iP]->SetMaximum(std::max(reso[iP]->GetMaximum(),resoFit[iP]->GetMaximum()));
      resoFit[iP]->GetXaxis()->SetLabelSize(0.05);
      resoFit[iP]->GetXaxis()->SetTitleSize(0.05);
      resoFit[iP]->GetYaxis()->SetLabelSize(0.05);
      resoFit[iP]->GetYaxis()->SetTitleSize(0.05);
      //resoFit[iP]->GetXaxis()->SetTitleOffset(0.7);
      resoFit[iP]->GetYaxis()->SetTitleOffset(1.2);
      

      resoFit[iP]->Draw("ap");
      //reso[iP]->Draw("l");
      bool addNoiseTerm = false;
      if (iP>0) addNoiseTerm = true;
      TF1 *fitFunc2;
      if (doVsE){
      	fitFunc2 =new TF1("reso","sqrt([0]/x+[1])",resoFit[iP]->GetXaxis()->GetXmin(),resoFit[iP]->GetXaxis()->GetXmax());
      	if (addNoiseTerm) fitFunc2 =new TF1("reso","sqrt([0]/x+[1]+[2]/(x*x))",resoFit[iP]->GetXaxis()->GetXmin(),resoFit[iP]->GetXaxis()->GetXmax());
      }
      else {
      	fitFunc2 =new TF1("reso","sqrt([0]*x*x+[1])",resoFit[iP]->GetXaxis()->GetXmin(),resoFit[iP]->GetXaxis()->GetXmax());
      	if (addNoiseTerm) fitFunc2 =new TF1("reso","sqrt([0]*x*x+[1]+[2]*x*x*x*x)",resoFit[iP]->GetXaxis()->GetXmin(),resoFit[iP]->GetXaxis()->GetXmax());
      }
      fitFunc2->SetParameter(0,0.2);
      fitFunc2->SetParLimits(0,0,1);
      fitFunc2->SetParameter(1,0.01);
      fitFunc2->SetParLimits(1,0,1);
      if (addNoiseTerm) {
	if (!isEM) fitFunc2->SetParameter(2,0.18*0.18);
	else fitFunc2->SetParameter(2,0.06*0.06);
	fitFunc2->SetParLimits(2,0,2);
	if (!isEM) fitFunc2->FixParameter(2,0.18*0.18);
	else fitFunc2->FixParameter(2,0.06*0.06);
      }
      fitref->Draw("same");
      resoFit[iP]->Fit(fitFunc2,"RME+");
      resoFit[iP]->GetFunction("reso")->SetLineColor(color[iP]);
      fitFunc2->SetLineColor(color[iP]);
      fitFunc2->Draw("same");

      double sigmaStoch = sqrt(fitFunc2->GetParameter(0));
      double sigmaStochErr = fitFunc2->GetParError(0)/(2*sigmaStoch);
      double sigmaConst = sqrt(fitFunc2->GetParameter(1));
      double sigmaConstErr = fitFunc2->GetParError(1)/(2*sigmaConst);
      double sigmaNoise = 0;
      double sigmaNoiseErr = 0;
      if (addNoiseTerm) {
	sigmaNoise = sqrt(fitFunc2->GetParameter(2));
	sigmaNoiseErr = fitFunc2->GetParError(2)/(2*sigmaNoise);
      }
      lat.SetTextColor(1);
      sprintf(buf,"#frac{#sigma}{E} #propto #frac{s}{#sqrt{E}} #oplus c #oplus #frac{n}{E}");
      double Emin = doVsE?50 : 1/sqrt(genEn[nGenEn-1]);
      lat.DrawLatex(Emin,resoFit[iP]->GetYaxis()->GetXmax()*0.94,buf);
      sprintf(buf,"s=%3.3f #pm %3.3f",sigmaStoch,sigmaStochErr);
      lat.DrawLatex(Emin,resoFit[iP]->GetYaxis()->GetXmax()*0.84,buf);
      sprintf(buf,"c=%3.3f #pm %3.3f",sigmaConst,sigmaConstErr);
      lat.DrawLatex(Emin,resoFit[iP]->GetYaxis()->GetXmax()*0.74,buf);
      sprintf(buf,"chi2/NDF = %3.3f/%d = %3.3f",fitFunc2->GetChisquare(),fitFunc2->GetNDF(),fitFunc2->GetChisquare()/fitFunc2->GetNDF());
      //lat.DrawLatex(Emin,resoFit[iP]->GetYaxis()->GetXmax()*0.54,buf);
      if (addNoiseTerm) {
	sprintf(buf,"n=%3.2f",sigmaNoise);
	lat.DrawLatex(Emin,resoFit[iP]->GetYaxis()->GetXmax()*0.64,buf);
      }
      lat.SetTextColor(4);
      //sprintf(buf,"CALICE s=%3.3f, c=%3.3f, n=%3.3f",sqrt(fitref->GetParameter(0)),sqrt(fitref->GetParameter(1)),sqrt(fitref->GetParameter(2)));
      sprintf(buf,"CALICE");
      lat.DrawLatex(70,resoFit[iP]->GetYaxis()->GetXmin()*1.35,buf);
      lat.SetTextColor(1);
      lat.DrawLatex(20,resoFit[iP]->GetYaxis()->GetXmax()*1.03,"Geant4 QGSP_BERT, #pi^{-} gun");


      myc[iP]->Update();
      saveName.str("");
      saveName << plotDir << "/ResoE_" << Enames[iP];
      myc[iP]->Update();
      myc[iP]->Print((saveName.str()+".png").c_str());
      myc[iP]->Print((saveName.str()+".pdf").c_str());
      


    }//loop on points
    
    saveName.str("");
    saveName << plotDir << "/Linearity_all";// << Enames[iP];
    myc[nPts]->Update();
    myc[nPts]->Print((saveName.str()+".png").c_str());
    myc[nPts]->Print((saveName.str()+".pdf").c_str());
    
    myc[nPts+1]->Update();
    saveName.str("");
    saveName << plotDir << "/RMSE_all";// << Enames[iP];
    myc[nPts+1]->Update();
    myc[nPts+1]->Print((saveName.str()+".png").c_str());
    myc[nPts+1]->Print((saveName.str()+".pdf").c_str());
    
    myc[nPts+2]->Update();
    saveName.str("");
    saveName << plotDir << "/Ereso_all";// << Enames[iP];
    myc[nPts+2]->Update();
    myc[nPts+2]->Print((saveName.str()+".png").c_str());
    myc[nPts+2]->Print((saveName.str()+".pdf").c_str());
    
    //myc[nPts+1]->Divide(1,2);
    myc[nPts+3]->cd();
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    groffset->SetTitle("");
    groffset->GetXaxis()->SetTitle("Point");
    groffset->GetYaxis()->SetTitle("Offset");
    groffset->SetMarkerStyle(21);
    groffset->Draw("ap");
    saveName.str("");
    saveName << plotDir << "/Offset";// << Enames[iP];
    myc[nPts+3]->Update();
    myc[nPts+3]->Print((saveName.str()+".png").c_str());
    myc[nPts+3]->Print((saveName.str()+".pdf").c_str());
    
    
    myc[nPts+4]->cd();
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    grslope->SetTitle("");
    grslope->GetXaxis()->SetTitle("Point");
    grslope->GetYaxis()->SetTitle("Slope");
    grslope->SetMarkerStyle(21);
    grslope->Draw("ap");
    saveName.str("");
    saveName << plotDir << "/Slope";// << Enames[iP];
    myc[nPts+4]->Update();
    myc[nPts+4]->Print((saveName.str()+".png").c_str());
    myc[nPts+4]->Print((saveName.str()+".pdf").c_str());
    
    myc[nPts+5]->cd();
    leg->Draw();
    saveName.str("");
    saveName << plotDir << "/Legend";// << Enames[iP];
    myc[nPts+5]->Update();
    myc[nPts+5]->Print((saveName.str()+".png").c_str());
    myc[nPts+5]->Print((saveName.str()+".pdf").c_str());
    
  }//calib
  
  return 0;
}
