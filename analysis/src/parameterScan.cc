#include "parameterScan.hh"
#include "TStyle.h"
#include "TLatex.h"
#include "TGraph.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>


void parameterScan::process30layers(){
  
  const unsigned nBlocks = 5;
  std::vector<unsigned>idx;
  idx.resize(nBlocks,0);
  for (unsigned iB(0); iB<nBlocks;++iB){
    idx[iB] = 3*iB;
  }
  
  for (iS[0]=0; iS[0]<nSteps;++iS[0]){
    std::cout << " Processing iS0 = " << iS[0] << std::endl;
    iS[1] = iS[0];
    iS[2] = iS[0];
    for (iS[3]=0; iS[3]<nSteps;++iS[3]){
      iS[4] = iS[3];
      iS[5] = iS[3];
      for (iS[6]=0; iS[6]<nSteps;++iS[6]){
	iS[7] = iS[6];
	iS[8] = iS[6];
	for (iS[9]=0; iS[9]<nSteps;++iS[9]){
	  iS[10] = iS[9];
	  iS[11] = iS[9];
	  for (iS[12]=0; iS[12]<nSteps;++iS[12]){
	    iS[13] = iS[12];
	    iS[14] = iS[12];
	    
	    processModel(nBlocks,idx);
	    
	  }//loop on S12
	}//loop on S9
      }//loop on S6
    }//loop on S3
  }//loop on S0
  
}//process30lay

void parameterScan::process28layers(){
  
  const unsigned nBlocks = 5;
  std::vector<unsigned>idx;
  idx.resize(nBlocks,0);
  for (unsigned iB(0); iB<nBlocks;++iB){
    idx[iB] = 3*iB;
  }
  
  for (iS[0]=0; iS[0]<nSteps;++iS[0]){
    std::cout << " Processing iS0 = " << iS[0] << std::endl;
    iS[1] = iS[0];
    iS[2] = iS[0];
    for (iS[3]=0; iS[3]<nSteps;++iS[3]){
      iS[4] = iS[3];
      iS[5] = iS[3];
      for (iS[6]=0; iS[6]<nSteps;++iS[6]){
	iS[7] = iS[6];
	iS[8] = iS[6];
	for (iS[9]=0; iS[9]<nSteps;++iS[9]){
	  iS[10] = iS[9];
	  iS[11] = iS[9];
	  for (iS[12]=0; iS[12]<nSteps;++iS[12]){
	    iS[13] = iS[12];
	    
	    processModel(nBlocks,idx);
	    
	  }//loop on S12
	}//loop on S9
      }//loop on S6
    }//loop on S3
  }//loop on S0
  
}//process28lay

void parameterScan::process26layers(){
  
  const unsigned nBlocks = 5;
  std::vector<unsigned>idx;
  idx.resize(nBlocks,0);
  for (unsigned iB(0); iB<nBlocks;++iB){
    idx[iB] = 3*iB;
  }
  
  for (iS[0]=0; iS[0]<nSteps;++iS[0]){
    std::cout << " Processing iS0 = " << iS[0] << std::endl;
    iS[1] = iS[0];
    iS[2] = iS[0];
    for (iS[3]=0; iS[3]<nSteps;++iS[3]){
      iS[4] = iS[3];
      iS[5] = iS[3];
      for (iS[6]=0; iS[6]<nSteps;++iS[6]){
	iS[7] = iS[6];
	iS[8] = iS[6];
	for (iS[9]=0; iS[9]<nSteps;++iS[9]){
	  iS[10] = iS[9];
	  iS[11] = iS[9];
	  for (iS[12]=0; iS[12]<nSteps;++iS[12]){
	    
	    processModel(nBlocks,idx);
	    
	  }//loop on S12
	}//loop on S9
      }//loop on S6
    }//loop on S3
  }//loop on S0
  
}//process26lay

void parameterScan::process24layers(){
  
  const unsigned nBlocks = 4;
  std::vector<unsigned>idx;
  idx.resize(nBlocks,0);
  for (unsigned iB(0); iB<nBlocks;++iB){
    idx[iB] = 3*iB;
  }
  //to not break conditions in processModel method which expects 5 blocks...
  idx.push_back(9);
  
  for (iS[0]=0; iS[0]<nSteps;++iS[0]){
    std::cout << " Processing iS0 = " << iS[0] << std::endl;
    iS[1] = iS[0];
    iS[2] = iS[0];
    for (iS[3]=0; iS[3]<nSteps;++iS[3]){
      iS[4] = iS[3];
      iS[5] = iS[3];
      for (iS[6]=0; iS[6]<nSteps;++iS[6]){
	iS[7] = iS[6];
	iS[8] = iS[6];
	for (iS[9]=0; iS[9]<nSteps;++iS[9]){
	  iS[10] = iS[9];
	  iS[11] = iS[9];

	    processModel(nBlocks,idx);
	    
	}//loop on S9
      }//loop on S6
    }//loop on S3
  }//loop on S0
  
  
}//process24lay

void parameterScan::process22layers(){
  
  const unsigned nBlocks = 5;
  std::vector<unsigned>idx;
  idx.resize(nBlocks,0);
  for (unsigned iB(0); iB<nBlocks;++iB){
    idx[iB] = 2*iB;
  }
  
  for (iS[0]=0; iS[0]<nSteps;++iS[0]){
    std::cout << " Processing iS0 = " << iS[0] << std::endl;
    iS[1] = iS[0];
    for (iS[2]=0; iS[2]<nSteps;++iS[2]){
      iS[3] = iS[2];
      for (iS[4]=0; iS[4]<nSteps;++iS[4]){
	iS[5] = iS[4];
	for (iS[6]=0; iS[6]<nSteps;++iS[6]){
	  iS[7] = iS[6];
	  for (iS[8]=0; iS[8]<nSteps;++iS[8]){
	    iS[9] = iS[8];
	    iS[10] = iS[8];
	    processModel(nBlocks,idx);
	  }//loop on S8
	}//loop on S6
      }//loop on S4
    }//loop on S2
  }//loop on S0
  
}//process22lay

void parameterScan::process20layers(){
  
  const unsigned nBlocks = 5;
  std::vector<unsigned>idx;
  idx.resize(nBlocks,0);
  for (unsigned iB(0); iB<nBlocks;++iB){
    idx[iB] = 2*iB;
  }
  
  for (iS[0]=0; iS[0]<nSteps;++iS[0]){
    std::cout << " Processing iS0 = " << iS[0] << std::endl;
    iS[1] = iS[0];
    for (iS[2]=0; iS[2]<nSteps;++iS[2]){
      iS[3] = iS[2];
      for (iS[4]=0; iS[4]<nSteps;++iS[4]){
	iS[5] = iS[4];
	for (iS[6]=0; iS[6]<nSteps;++iS[6]){
	  iS[7] = iS[6];
	  for (iS[8]=0; iS[8]<nSteps;++iS[8]){
	    iS[9] = iS[8];
	    processModel(nBlocks,idx);
	  }//loop on S8
	}//loop on S6
      }//loop on S4
    }//loop on S2
  }//loop on S0
  
}//process22lay

bool parameterScan::processModel(const unsigned nBlocks,
				 const std::vector<unsigned> & idx){
  totModels++;
  
  //set thick array
  double totthick = length;
  double xtot=X0tot;
  double ltot=L0tot;
  std::vector<double> wThick;
  std::vector<double> pbThick;
  wThick.resize(nLayers/2,0);
  pbThick.resize(nLayers/2,0);
  
  for (unsigned iL(0); iL<nLayers/2;++iL){
    wThick[iL] = iS[iL]*stepSizeW;
    pbThick[iL] = iS[iL]*stepSizePb;
    totthick += wThick[iL]+pbThick[iL];
    xtot += wThick[iL]/x0w+pbThick[iL]/x0pb;
    ltot += wThick[iL]/l0w+pbThick[iL]/l0pb;
  }
  
  //filter
  if (totthick > maxThick || xtot<minX0 || ltot>maxLambda) return false;
  //reject obvious wrong models...
  if (iS[idx[0]]>iS[idx[1]] && iS[idx[1]]>iS[idx[2]] && iS[idx[2]]>iS[idx[3]] && iS[idx[3]]>iS[idx[4]]) return false;
  if (iS[idx[0]]>iS[idx[4]] || iS[idx[1]]>iS[idx[4]]) return false;
  //if (iS[idx[0]]>iS[idx[3]] || iS[idx[1]]>iS[idx[3]]) continue;
  
  //std::cout << " - Model kept: Total length = " << totthick << "mm, nX0 = " << xtot << ", nL0 = " << ltot << std::endl;
  
  //std::cout << " ---- Model Kept ! ---- " << std::endl
  //	      << "Total length = " << totthick
  //	      << "mm, nX0 = " << xtot
  //	      << ", nL0 = " << ltot << std::endl;
  //for (unsigned iL(0); iL<nLayers/2;++iL){
  //std::cout << "Layer pair " << iL << ": " << iS[iL] << " " << wThick[iL] << " " << pbThick[iL] << std::endl;
  //}
  validModels++;
  
  float block[nBlocks];
  float extraW[nBlocks];
  for (unsigned iB(0); iB<nBlocks;++iB){
    block[iB] = iB;
    extraW[iB] = iS[idx[iB]]*stepSizeW/x0w;
  }
  
  TGraph *tmp = new TGraph(nBlocks,block,extraW);
  tmp->SetMaximum(nSteps*stepSizeW/x0w);
  tmp->SetMinimum(0);
  if (iS[idx[0]]<iS[idx[1]] && iS[idx[1]]<iS[idx[2]] && iS[idx[2]]<iS[idx[3]] && iS[idx[3]]<=iS[idx[4]]) {
    tmp->SetTitle("1<2<3<4<5;block;extra thickness (X0)");
    myc[0]->cd();
    tmp->SetLineColor(counter[0]%9+1);
    tmp->Draw(first[0]?"AL":"L");
    first[0] = false;
    counter[0]++;
    /*std::cout << " --- " << xtot << " " << totthick << " " << ltot 
      << " " << iS[0] << " " << iS[3] << " " << iS[6]
      << " " << iS[9] << " " << iS[12]
      << std::endl;*/
  }
  else if (iS[idx[0]]>iS[idx[1]] && iS[idx[1]]>iS[idx[2]] && iS[idx[2]]<iS[idx[3]] && iS[idx[3]]<=iS[idx[4]]) {
    tmp->SetTitle("1>2>3<4<5;block;extra thickness (X0)");
    myc[1]->cd();
    tmp->SetLineColor(counter[1]%9+1);
    tmp->Draw(first[1]?"AL":"L");
    first[1] = false;
    counter[1]++;
    /*std::cout << " --- " << xtot << " " << totthick << " " << ltot 
      << " " << iS[0] << " " << iS[3] << " " << iS[6]
      << " " << iS[9] << " " << iS[12]
      << std::endl;*/
  }
  else if (iS[idx[0]]<iS[idx[1]] && iS[idx[1]]<iS[idx[2]] && iS[idx[2]]>iS[idx[3]] && iS[idx[3]]>=iS[idx[4]]) {
    tmp->SetTitle("1<2<3>4>5;block;extra thickness (X0)");
    myc[2]->cd();
    tmp->SetLineColor(counter[2]%9+1);
    tmp->Draw(first[2]?"AL":"L");
    first[2] = false;
    counter[2]++;
    /*std::cout << " --- " << xtot << " " << totthick << " " << ltot 
      << " " << iS[0] << " " << iS[3] << " " << iS[6]
      << " " << iS[9] << " " << iS[12]
      << std::endl;*/
  }
  else {
    tmp->SetTitle("others;block;extra thickness (X0)");
    myc[nC-1]->cd();
    tmp->SetLineColor(counter[nC-1]%9+1);
    tmp->Draw(first[nC-1]?"AL":"L");
    first[nC-1] = false;
    counter[nC-1]++;
  }
  //tmp->GetXaxis()->SetTitle();
  
  isInserted = modelMap.insert(std::pair<double,unsigned>(static_cast<unsigned>(xtot*100000)/100000.,1));
  if (!isInserted.second) isInserted.first->second += 1;
  //else std::cout << " ---- element inserted: " << xtot << " " << totthick << " " << ltot 
  //	       << " " << iS[0] << " " << iS[3] << " " << iS[6]
  //	       << " " << iS[9] << " " << iS[12]
  //	       << std::endl;
  
  /*if (validModels>500) {
    std::cout << " Warning!!! Too many valid models ! Stopping the scan..." << std::endl;
    std::cout << " nValid = " << validModels << "/" << totModels << std::endl;
    std::cout << " ---- Last Model ---- " << std::endl
    << "Total length = " << totthick
    << "mm, nX0 = " << xtot
    << ", nL0 = " << ltot << std::endl;
    for (unsigned iL(0); iL<nLayers/2;++iL){
    std::cout << "Layer pair " << iL << ": " << iS[iL] << " " << minWthick+wThick[iL] << "mm " << minPbthick+pbThick[iL] << "mm" << std::endl;
    }
    return 1;
    }*/
  return true;
}

void parameterScan::print(){

  std::cout << " ------------------------------------- " << std::endl
	    << " ---- Printout of scan parameters ---- " << std::endl
	    << " ------------------------------------- " << std::endl;
  std::cout <<  "nLayers             = " << nLayers << std::endl;
  std::cout <<  "max thickness (mm)  = " << maxThick << std::endl;
  std::cout <<  "minX0               = " << minX0 << std::endl;
  std::cout <<  "max lambda_nucl_int = " << maxLambda << std::endl;
  
  std::cout <<  "Step size for W (mm / X0)  = " << stepSizeW << " / " << stepSizeW/x0w << std::endl;
  std::cout <<  "Step size for Pb (mm / X0) = " << stepSizePb  << " / " << stepSizePb/x0pb << std::endl;
  
  std::cout <<  "Number of steps = " << nSteps << std::endl;
  
  std::cout << "Number of layer pair: " << iS.size()  << std::endl;

  std::cout << " ------------------------------------- " << std::endl;
  std::cout << " ------------------------------------- " << std::endl;

}
