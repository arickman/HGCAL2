#include "HGCSSDigitisation.hh"
#include <cmath>
#include <sstream>
#include <iostream>

unsigned HGCSSDigitisation::nRandomPhotoElec(const double & aMipE) {
	double mean = aMipE * npe_;
	int result = rndm_.Poisson(mean);
	if (result < 0) {
		std::cout
				<< "WARNING!! HGCSSDigitisation::nRandomPhotoElec Poisson return negative number!! "
				<< aMipE << " " << mean << " " << result << std::endl;

	}
	return static_cast<unsigned>(result);
}

unsigned HGCSSDigitisation::nPixels(const double & aMipE) {
	unsigned npe = nRandomPhotoElec(aMipE);
	double x = exp(-1. * npe / nTotal_);
	double res = nTotal_ * 1.0 * (1 - x) / (1 - crossTalk_ * x);
	unsigned npix = static_cast<unsigned>(res);
	return npix;
}

unsigned HGCSSDigitisation::positiveRandomGaus(const unsigned & mean) {
	double result = rndm_.Gaus(mean, sigmaPix_);
	if (result < 0)
		result = 0; //positiveRandomGaus(mean);
	if (result >= nTotal_)
		result = nTotal_ - 1;
	return static_cast<unsigned>(result);
}

double HGCSSDigitisation::mipCor(const double & aMipE, const double & posx,
		const double & posy, const double & posz) {
	double costheta = fabs(posz)
			/ sqrt(posz * posz + posx * posx + posy * posy);
	if (costheta > 0)
		return aMipE * costheta;
	return aMipE;
}

double HGCSSDigitisation::digiE(const double & aMipE, TH2F * & p_pixvspe,
		TH1F * & p_npixels, TH1F * & p_npixelssmeared, TH2F * & p_outvsnpix) {
	if (aMipE == 0)
		return 0;
	unsigned npix = nPixels(aMipE);
	if (p_pixvspe)
		p_pixvspe->Fill(nRandomPhotoElec(aMipE), npix);
	if (p_npixels)
		p_npixels->Fill(npix);
	unsigned npixsmear = positiveRandomGaus(npix);
	if (p_npixelssmeared)
		p_npixelssmeared->Fill(npixsmear);
	double result = nTotal_ * 1.0 / npe_
			* log((nTotal_ - crossTalk_ * npixsmear) / (nTotal_ - npixsmear));
	if (result < 0) {
		std::cout << "WARNING!! HGCSSDigitisation::digiE negative result!! "
				<< npix << " " << npixsmear << " " << nTotal_ << " " << result
				<< std::endl;
		result = 0;
	}
	p_outvsnpix->Fill(npixsmear, result);
	return result;
}

double HGCSSDigitisation::digiE(const double & aMipE) {
	if (aMipE == 0)
		return 0;
	unsigned npix = nPixels(aMipE);
	unsigned npixsmear = positiveRandomGaus(npix);
	double result = nTotal_ * 1.0 / npe_
			* log((nTotal_ - crossTalk_ * npixsmear) / (nTotal_ - npixsmear));
	if (result < 0) {
		std::cout << "WARNING!! HGCSSDigitisation::digiE negative result!! "
				<< npix << " " << npixsmear << " " << nTotal_ << " " << result
				<< std::endl;
		result = 0;
	}
	return result;
}

double HGCSSDigitisation::ipXtalk(const std::vector<double> & aSimEvec) {
	double result = 0;
	const unsigned nEdges = aSimEvec.size() - 1;
	//give away X% per edge
	result = aSimEvec[0] * (1 - ipXtalk_ * 4); //nEdges);
	//get X% back from neighbours
	for (unsigned i(1); i < nEdges + 1; ++i) {
		result += ipXtalk_ * aSimEvec[i];
	}

	return result;
}

void HGCSSDigitisation::addNoise(double & aDigiE, const unsigned & alay,
		TH1F * & hist) {
	bool print = false;
	//if (aDigiE>0) print = true;
	if (print)
		std::cout << "HGCSSDigitisation::addNoise " << aDigiE << " ";
	double lNoise = rndm_.Gaus(0, noise_[alay]);
	if (hist)
		hist->Fill(lNoise);
	aDigiE += lNoise;
	if (aDigiE < 0)
		aDigiE = 0;
	if (print)
		std::cout << lNoise << " " << aDigiE << std::endl;
}

unsigned HGCSSDigitisation::adcConverter(double eMIP, DetectorEnum adet) {
	if (eMIP < 0)
		eMIP = 0;
	double eADC = static_cast<unsigned>(eMIP * mipToADC_[adet]);
	if (eADC > maxADC_[adet])
		eADC = maxADC_[adet];
	return eADC;
}

double HGCSSDigitisation::adcToMIP(const unsigned adcCounts, DetectorEnum adet,
		const bool smear) {
	double lE = adcCounts * 1.0 / mipToADC_[adet];
	if (!smear)
		return lE;
	return rndm_.Gaus(lE, gainSmearing_[adet] * lE);
}

double HGCSSDigitisation::MIPtoGeV(const HGCSSSubDetector & adet,
		const double & aMipE) {
	double lE = aMipE * adet.absWeight * adet.gevWeight - adet.gevOffset;
	return lE;
}

double HGCSSDigitisation::sumBins(const std::vector<TH2D *> & aHistVec,
		const double & aMipThresh) {
	double energy = 0;
	for (unsigned iL(0); iL < aHistVec.size(); ++iL) {
		for (int ix(1); ix < aHistVec[iL]->GetNbinsX() + 1; ++ix) {
			for (int iy(1); iy < aHistVec[iL]->GetNbinsY() + 1; ++iy) {
				double eTmp = aHistVec[iL]->GetBinContent(ix, iy);
				if (eTmp > aMipThresh)
					energy += eTmp;
			}
		}
	}
	return energy;
}
void HGCSSDigitisation::Print(std::ostream & aOs) const {
	aOs << "====================================" << std::endl
			<< "=== INIT DIGITISATION PARAMETERS ===" << std::endl
			<< "====================================" << std::endl
			<< " = Random seed: " << rndm_.GetSeed() << std::endl
			<< " = Nphoto-electrons: " << npe_ << std::endl << " = cross-talk: "
			<< crossTalk_ << std::endl << " = Npixels total: " << nTotal_
			<< std::endl << " = sigmaPixel: " << sigmaPix_ << std::endl
			//<< " = sigmaNoise: ECAL " << noise_.find(DetectorEnum::ECAL)->second << ", FHCAL " << noise_.find(DetectorEnum::FHCAL)->second << ", BHCAL " << noise_.find(DetectorEnum::BHCAL)->second << std::endl
			<< " = MIPtoADC conversions: ECAL "
			<< mipToADC_.find(DetectorEnum::FECAL)->second << ", FHCAL "
			<< mipToADC_.find(DetectorEnum::FHCAL)->second << std::endl
			<< " = Time cut: ECAL "
			<< timeCut_.find(DetectorEnum::FECAL)->second << ", FHCAL "
			<< timeCut_.find(DetectorEnum::FHCAL)->second << ", BHCAL "
			<< timeCut_.find(DetectorEnum::BHCAL1)->second << std::endl
			<< " = Intercalibration: ECAL "
			<< gainSmearing_.find(DetectorEnum::FECAL)->second << ", FHCAL "
			<< gainSmearing_.find(DetectorEnum::FHCAL)->second << ", BHCAL "
			<< gainSmearing_.find(DetectorEnum::BHCAL1)->second << std::endl
			<< "====================================" << std::endl;
}
;
