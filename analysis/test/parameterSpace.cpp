#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <iomanip>
#include <boost/algorithm/string.hpp>

#include "TH2F.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TLatex.h"
#include "parameterScan.hh"

struct SamplingLayer {
  double thick;
  double X0;
  double L0;
};

double getX0(std::string element){
  //unit = mm
  if (element == "Cu") return 14.3558;
  else if (element == "Pb") return 5.61253;
  else if (element == "W") return 3.50418;
  else if (element == "Air") return 303921;
  else if (element == "PCB") return 187.31;
  else if (element == "Si") return 93.6607;
  //  else if (element == "") return ;
  else {
    std::cerr << "Element not known, please add !" << std::endl;
    exit(1);
  }
};

double getL0(std::string element){
  //unit = mm
  if (element == "Cu") return 153.2;
  else if (element == "Pb") return 175.9;
  else if (element == "W") return 99.46;
  else if (element == "Air") return 747749;
  else if (element == "PCB") return 700;
  else if (element == "Si") return 465.2;
  //  else if (element == "") return ;
  else {
    std::cerr << "Element not known, please add !" << std::endl;
    exit(1);
  }
};

void fillStructure(const std::vector<double>& thick,
		   const std::vector<std::string>& elements,
		   std::vector<SamplingLayer> & fullstruct){
 
  for (unsigned i(0); i<thick.size();++i){
    SamplingLayer ltmp;
    ltmp.thick=thick[i];
    ltmp.X0=getX0(elements[i]);
    ltmp.L0=getL0(elements[i]);
    fullstruct.push_back(ltmp);
  }
 
}

void calculateLength(const std::vector<SamplingLayer> & fullstruct,
		     double & length,
		     double & X0tot,
		     double & L0tot){
  length = 0;
  X0tot = 0;
  L0tot = 0;
  for (unsigned i(0); i<fullstruct.size();++i){
    length += fullstruct[i].thick;
    X0tot += fullstruct[i].thick/fullstruct[i].X0;
    L0tot += fullstruct[i].thick/fullstruct[i].L0;
  }
}

void fillFullStruct(std::vector<SamplingLayer> & fullstruct,
		    const std::vector<double> & wThick,
		    const std::vector<double> & pbThick,
		    const unsigned & nLayers){
  
  //all in mm
  double airThick = 2;
  double pcbThick = 1.2;
  
  //first and last layers    
  std::vector<double> lThick1;
  std::vector<std::string> lEle1;
  lThick1.push_back(0);lEle1.push_back("Cu"); 
  lThick1.push_back(0);lEle1.push_back("W"); 
  lThick1.push_back(0.5);lEle1.push_back("Cu"); 
  lThick1.push_back(airThick);lEle1.push_back("Air"); 
  lThick1.push_back(pcbThick);lEle1.push_back("PCB"); 
  lThick1.push_back(0.3);lEle1.push_back("Si"); 
  
  fillStructure(lThick1,lEle1,fullstruct);
  
  std::vector<double> lThickL;
  std::vector<std::string> lEleL;
  lThickL.push_back(3);lEleL.push_back("Cu"); 
  lThickL.push_back(1);lEleL.push_back("Pb"); 
  lThickL.push_back(wThick[0]);lEleL.push_back("W"); 
  lThickL.push_back(0.5);lEleL.push_back("Cu"); 
  lThickL.push_back(airThick);lEleL.push_back("Air"); 
  lThickL.push_back(pcbThick);lEleL.push_back("PCB"); 
  lThickL.push_back(0.3);lEleL.push_back("Si"); 
  
  std::vector<double> lThickR;
  std::vector<std::string> lEleR;
  lThickR.push_back(3);lEleR.push_back("Cu"); 
  lThickR.push_back(pbThick[0]);lEleR.push_back("Pb"); 
  lThickR.push_back(3);lEleR.push_back("Cu"); 
  lThickR.push_back(0.3);lEleR.push_back("Si"); 
  lThickR.push_back(pcbThick);lEleR.push_back("PCB"); 
  lThickR.push_back(airThick);lEleR.push_back("Air");                     
  
  //second layer with Cu/Pb in front
  fillStructure(lThickL,lEleL,fullstruct);
  fillStructure(lThickR,lEleR,fullstruct);
  
  //reset to 0.5 Cu and no lead for following layers
  lThickL[0] = 0.5;
  lThickL[1] = 0;

  for(unsigned i=0; i<(nLayers-4)/2; i++) {
    //if (i>3 && i<9) { wThick[1+i]=2.8;pbThick[1+i]=2.1;}
    //else if (i>3) { wThick[1+i]=4.2;pbThick[1+i]=4.4;}
    lThickL[2] = wThick[1+i];//1.75
    lThickR[1] = pbThick[1+i];//1;    
    fillStructure(lThickL,lEleL,fullstruct);
    fillStructure(lThickR,lEleR,fullstruct);
  }
  
  //last layer: add Cu+W+Cu...

  //wThick[1+(nLayers-4)/2]=4.2;

  lThick1[0] = 0.5;
  lThick1[1] = wThick[1+(nLayers-4)/2];//4.2;
  //add last structure layers
  
  lThick1.push_back(3);lEle1.push_back("Cu"); 
  lThick1.push_back(1);lEle1.push_back("Pb");
  fillStructure( lThick1,lEle1,fullstruct);


}

//int parameterSpace(){
int main(int argc, char** argv){//main  

  if (argc < 2) {
    std::cout << " Usage: " 
	      << argv[0] << " <nLayers>"
	      << std::endl;
    return 1;
  }
  unsigned nLayers = atoi(argv[1]);

  if (nLayers%2 != 0) nLayers+=1;
  if (nLayers<20) nLayers==20;
  if (nLayers>30) nLayers==30;

  parameterScan scan;
  scan.nLayers = nLayers;

  std::vector<double> wThick;
  std::vector<double> pbThick;

  TString plotDir = "PLOTS/ParameterScan/";
  plotDir += nLayers;
  plotDir += "layers/";
  
  scan.maxThick = 300;
  scan.minX0 = 26;
  scan.maxLambda = 1.4;

  scan.stepSizeW = 1.*0.35;//multiples of 0.1X0
  scan.stepSizePb = 1.*0.56;//multiples of 0.1X0

  scan.x0w = getX0("W");
  scan.x0pb = getX0("Pb");
  scan.l0w = getL0("W");
  scan.l0pb = getL0("Pb");


  //TO do: add quick first check based on smaller and larger values :/
  unsigned nSteps = 15;
  //if (nLayers == 30) nSteps = 14;
  scan.nSteps = nSteps;//8;
  
  //  wThick.resize(nLayers/2,1.75);
  double minWthick = 1.5;
  double minPbthick = 1;
  wThick.resize(nLayers/2,minWthick);
  pbThick.resize(nLayers/2,minPbthick);

  std::vector<SamplingLayer>fullstruct;
  fullstruct.reserve(nLayers);
 
  fillFullStruct(fullstruct,
		 wThick,pbThick,
		 nLayers);
  
  double length=0;
  double X0tot=0;
  double L0tot=0;
  calculateLength(fullstruct,length,X0tot,L0tot);

  scan.length = length;
  scan.X0tot = X0tot;
  scan.L0tot = L0tot;

  (scan.iS).resize(nLayers/2,0);
  scan.validModels = 0;
  scan.totModels = 0;


  scan.print();


  std::cout << " -- Starting point:" << std::endl
	    << " ------  Total length = " << length
	    << "mm , nX0 = " << X0tot
	    << ", nL0 = " << L0tot << std::endl;
  
  if (length>scan.maxThick || L0tot>scan.maxLambda){
    std::cout << " ERROR! Minimum model is already out from max length..." << std::endl;
    return 1;
  }
  
  std::cout << " -- Max scanned: " << std::endl;
  double tmptot = length;
  double tmpxtot=X0tot;
  double tmpltot=L0tot;
  for (unsigned iL(0); iL<nLayers/2;++iL){
    double tmpw = (nSteps-1)*scan.stepSizeW;
    double tmppb = (nSteps-1)*scan.stepSizePb;
    tmptot += tmpw+tmppb;
    tmpxtot += tmpw/scan.x0w+tmppb/scan.x0pb;
    tmpltot += tmpw/scan.l0w+tmppb/scan.l0pb;
  }
  std::cout << " ------  Total length = " << tmptot
	    << "mm , nX0 = " << tmpxtot
	    << ", nL0 = " << tmpltot << std::endl;
  
  if (tmpxtot < scan.minX0) {
    std::cout << " ERROR! Maximum model is already out from min X0..." << std::endl;
    return 1;
  }

  //check that middle between both is acceptable
  double midlength = (tmptot+length)/2;
  double midx0 = (tmpxtot+X0tot)/2.;
  double midl0 = (tmpltot+L0tot)/2.;

  std::cout  << " -- Middle scanned: " << std::endl
	     << " ------  Total length = " << midlength
	     << "mm , nX0 = " << midx0
	     << ", nL0 = " << midl0 << std::endl;
  
  if (midlength>scan.maxThick || midl0>scan.maxLambda || midx0<scan.minX0) {
    std::cout << " WARNING! Middle model is already out... Please check it is OK..." << std::endl;
    //return 1;
  }

  //return 1;

  const unsigned nC = 4;
  scan.nC=nC;
  scan.myc.resize(4,0);
  scan.first.resize(4,0);
  scan.counter.resize(4,0);
  for (unsigned iC(0);iC<nC;++iC){
    std::ostringstream label;
    label << "myc" << iC;
    scan.myc[iC] = new TCanvas(label.str().c_str(),
			       label.str().c_str(),
			       1);
    scan.first[iC] = true;
    scan.counter[iC] = 0;
  }

  if (nLayers==30)
    scan.process30layers();
  else if (nLayers==28)
    scan.process28layers();
  else if (nLayers==26)
    scan.process26layers();
  else if (nLayers==24)
    scan.process24layers();
  else if (nLayers==22)
    scan.process22layers();
  else
    scan.process20layers();


  std::cout << " Found " << scan.validModels << " valid models out of " << scan.totModels << std::endl;
  
  std::cout << " -- number of uniq X0tot: " << scan.modelMap.size() << std::endl;
  std::map<double,unsigned>::iterator iter;
  unsigned count = 0;
  for (iter=scan.modelMap.begin(); iter!=scan.modelMap.end();++iter,++count){
    std::cout << count << " " << std::setprecision(10) << iter->first << " " << iter->second << std::endl;
  }

  //  summary->Draw("colz");

  std::cout << " -- Found " << scan.counter[0] << " models for ordered" << std::endl;
  std::cout << " -- Found " << scan.counter[1] << " models for down-up" << std::endl;
  std::cout << " -- Found " << scan.counter[2] << " models for up-down" << std::endl;
  std::cout << " -- Found " << scan.counter[3] << " models for others" << std::endl;

  gStyle->SetOptStat(0);

  TLatex lat;
  char buf[500];
  for (unsigned iC(0);iC<nC;++iC){
    sprintf(buf,"Number of models: %d",scan.counter[iC]);
    scan.myc[iC]->cd();
    lat.SetTextSize(0.05);
    lat.DrawLatexNDC(0.15,0.85,buf);
    scan.myc[iC]->Update();
  }

  scan.myc[0]->Print(plotDir+"summary_ordered.pdf");
  scan.myc[1]->Print(plotDir+"summary_downup.pdf");
  scan.myc[2]->Print(plotDir+"summary_updown.pdf");
  scan.myc[3]->Print(plotDir+"summary_others.pdf");

  return 0;
  
}//main
