#include "TTree.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TMath.h"
#include <string>
#include <cstring>
#include <sstream>
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TLine.h"
#include "TBox.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "THStack.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include "TFitResultPtr.h"
const char* effSigma_cstr;
double effSigma_val;
Double_t effSigmaMacro(TH1 * hist)
{
TAxis *xaxis = hist->GetXaxis();
Int_t nb = xaxis->GetNbins();
// std::cout << nb << std::endl;
if(nb < 10) {
std::cout << "effsigma: Not a valid histo. nbins = " << nb << std::endl;
return 0.;
}
Double_t bwid = xaxis->GetBinWidth(1);
if(bwid == 0) {
std::cout << "effsigma: Not a valid histo. bwid = " << bwid << std::endl;
return 0.;
}
Double_t xmin = xaxis->GetXmin();
Double_t ave = hist->GetMean();
Double_t rms = hist->GetRMS();
Double_t total=0.;
for(Int_t i=0; i<nb+2; i++) {
total+=hist->GetBinContent(i);
}
// if(total < 100.) {
// std::cout << "effsigma: Too few entries " << total << std::endl;
// return 0.;
// }
Int_t ierr=0;
Int_t ismin=999;
Double_t rlim=0.683*total;
Int_t nrms=rms/(bwid); // Set scan size to +/- rms
if(nrms > nb/10) nrms=nb/10; // Could be tuned...
Double_t widmin=9999999.;
for(Int_t iscan=-nrms;iscan<nrms+1;iscan++) { // Scan window centre
Int_t ibm=(ave-xmin)/bwid+1+iscan;
Double_t x=(ibm-0.5)*bwid+xmin;
Double_t xj=x;
Double_t xk=x;
Int_t jbm=ibm;
Int_t kbm=ibm;
Double_t bin=hist->GetBinContent(ibm);
total=bin;
for(Int_t j=1;j<nb;j++){
if(jbm < nb) {
jbm++;
xj+=bwid;
bin=hist->GetBinContent(jbm);
total+=bin;
if(total > rlim) break;
}
else ierr=1;
if(kbm > 0) {
kbm--;
xk-=bwid;
bin=hist->GetBinContent(kbm);
total+=bin;
if(total > rlim) break;
}
else ierr=1;
}
Double_t dxf=(total-rlim)*bwid/bin;
Double_t wid=(xj-xk+bwid-dxf)*0.5;
if(wid < widmin) {
widmin=wid;
ismin=iscan;
}
}
if(ismin == nrms || ismin == -nrms) ierr=3;
if(ierr != 0) std::cout << "effsigma: Error of type " << ierr << std::endl;
//std::cout <<widmin << std::endl;
std::ostringstream ostr1;
ostr1 << "#sigma_{eff} = " << std::fixed << widmin ;
// std::cout << "#sigma_{eff} = " << std::fixed << std::setprecision(3) << widmin << std::endl;
std::string effSigma_str = ostr1.str();
effSigma_cstr = effSigma_str.c_str();
// std::cout << effSigma_str << std::endl <<effSigma_cstr << std::endl;
effSigma_val=widmin;
return widmin;
}
