#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TF1.h"
#include "TString.h"
#include "TProfile.h"

#include "../../userlib/include/HGCSSRecoHit.hh"

void sumPlane(const char *glob)
{
  TFile *fout = TFile::Open("energysums.root","RECREATE");

  TH1D *hsumplane     = new TH1D("hsumplane","#Sigma_{all cells} E(cell);MIPs",55,-50,500); //80,0,800); 
  TH1D *hsum7cellE    = new TH1D("hsum7cellE","#Sigma_{7 central cells} E(cell);MIPs",55,-50,500); //80,0,800);
  TH1D *hcentralcellE = new TH1D("hcentralcellE","E(central cell); MIPs",55,-50,500); //80,0,800);

  TH1I *hnum1cell = new TH1I("hnum1cell","hnum1cell",10,0,10);
  TH1I *hnum7cell = new TH1I("hnum7cell","hnum7cell",10,0,10);
  TH1I *hnumplanecells = new TH1I("hnumplanecells","hnumplanecells",100,0,100);

  TChain *chain = new TChain("RecoTree");

  if (!chain->Add(glob)) {
    std::cerr << "couldn't find RecoTree in" << TString(glob) << std::endl;
    delete chain;
    return;
  }

  HGCSSRecoHitVec * rechits;
  chain->SetBranchAddress("HGCSSRecoHitVec",&rechits);

  unsigned nEvts = chain->GetEntries();
  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries

    if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;
    chain->GetEntry(ievt);

    double planeE=0;
    double centralcellE = 0;
    double sum7cellE = 0;
    int numplanecells=0;
    int num7cells=0;
    int num1cell=0;

    for (int i=0; i<rechits->size(); i++) {
      HGCSSRecoHit& hit=(*rechits)[i];

      //std::cout << "hit " << i << " layer " << rechits[i].layer() << std::endl;

      // only layer 0;
      if (hit.layer()) continue;

      planeE += hit.energy();
      numplanecells++;
      double celldistfromO = hit.get_x()*hit.get_x() + hit.get_y()*hit.get_y();
      std::cout << celldistfromO << std::endl;

      if (celldistfromO < 1) {
	num1cell++;
	centralcellE += hit.energy();
      }

      if (celldistfromO < 3*6.5*6.5) {
	num7cells++;
	sum7cellE += hit.energy();
      }
    }
    double MIPtoKeVsim = 54.8;
    double MIPadjust = MIPtoKeVsim/50.7;
    hsumplane->Fill(MIPadjust*planeE);
    hsum7cellE->Fill(MIPadjust*sum7cellE);
    hcentralcellE->Fill(MIPadjust*centralcellE);

    hnum1cell->Fill(num1cell);
    hnum7cell->Fill(num7cells);
    hnumplanecells->Fill(numplanecells);
      
  }//loop on entries

  fout->Write();

}//main
