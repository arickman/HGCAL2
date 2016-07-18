import sys
import ROOT
import tdrstyle
tdrstyle.setTDRStyle()
ROOT.gStyle.SetPadRightMargin(0.10);
ROOT.gStyle.SetOptFit(0 );
ROOT.gStyle.SetPalette(1);
import math
import array


##---------------------------------------
def main():

	fIn = ROOT.TFile("data/analyzed_tuple.root");
	t   = fIn.Get("hadrons");

	t.Print();

	h_layerHFlux  = ROOT.TH1F("h_layerHFlux","; hadron flux; Nevents*Nlayers", 100, 0, 4000);
	h_layerNFlux  = ROOT.TH1F("h_layerNFlux","; neutron flux; Nevents*Nlayers", 100, 0, 4000);
	h2_layerHFlux = ROOT.TH2F("h2_layerHFlux","; Nlayers; hadron flux", 26, 0.5, 26.5, 99, 4, 4000);
	h2_layerNFlux = ROOT.TH2F("h2_layerNFlux","; Nlayers; neutron flux", 26, 0.5, 26.5, 99, 4, 4000);
	h_summedHFlux = ROOT.TH1F("h_summedHFlux","; summed hadron flux; Nevents",100,0,4000);
	h_summedNFlux = ROOT.TH1F("h_summedNFlux","; summed neutron flux; Nevents",100,0,4000);
	h_layerSen    = ROOT.TH1F("h_layerSen","; layer Si Deposition; Nevents*Nlayers", 100,0,25);
	h_summedSen   = ROOT.TH1F("h_summedSen","; total Si Deposition; Nevents", 100,0,100);

	nent = t.GetEntries();
	for i in range(nent):
		
		if i > 10000000: break;
		if(i % (1 * nent/100) == 0):
			sys.stdout.write("\r[" + "="*int(20*i/nent) + " " + str(round(100.*i/nent,0)) + "% done");
			sys.stdout.flush();

		t.GetEntry(i);
	
		nlayers = t.caloLen;
		# print nlayers
		for i in range(nlayers):
			h_layerHFlux.Fill(t.layerHFlux[i]);
			h_layerNFlux.Fill(t.layerNFlux[i]);
			h2_layerHFlux.Fill(i+1, t.layerHFlux[i]);
			h2_layerNFlux.Fill(i+1, t.layerNFlux[i]);
			h_layerSen.Fill(t.layerSen[i]);

		h_summedHFlux.Fill( t.summedHFlux );
		h_summedNFlux.Fill( t.summedNFlux );
		h_summedSen.Fill(t.summedSen);

		# 	print t.layerHFlux[i]

	print "\n";

	makeCanvas1(h_layerHFlux,"plots/h_layerHFlux");
	makeCanvas1(h_layerNFlux,"plots/h_layerNFlux");
	makeCanvas1_2D(h2_layerHFlux,"plots/h2_layerHFlux");
	makeCanvas1_2D(h2_layerNFlux,"plots/h2_layerNFlux");
	makeCanvas1(h_summedHFlux,"plots/h_summedHFlux");
	makeCanvas1(h_summedNFlux,"plots/h_summedNFlux");
	makeCanvas1(h_layerSen,"plots/h_layerSen",False);
	makeCanvas1(h_summedSen,"plots/h_summedSen",False);

##---------------------------------------
def makeCanvas1(h,name,setlog=True):

	c = ROOT.TCanvas("c_"+h.GetName(),"c_"+h.GetName(),1000,800);
	h.SetLineWidth(2);
	h.Draw("hist");
	if setlog: ROOT.gPad.SetLogy();
	c.SaveAs(name+".pdf");

##---------------------------------------
def makeCanvas1_2D(h,name):

	c = ROOT.TCanvas("c_"+h.GetName(),"c_"+h.GetName(),1000,800);
	h.SetLineWidth(2);
	h.Draw("colz");
	ROOT.gPad.SetLogz();
	c.SaveAs(name+".pdf");

##---------------------------------------
if __name__ == '__main__':
        main();
##---------------------------------------        