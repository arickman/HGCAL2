#ifndef _hgcsssamplingsection_hh_
#define _hgcsssamplingsection_hh_
#include <iomanip>
#include <vector>
#include "Rtypes.h"
#include <sstream>
#include <map>

class HGCSSSamplingSection {

public:
	HGCSSSamplingSection() :
			volNb_(0), volX0trans_(0), voldEdx_(0), volLambdatrans_(0), sensDep_(
					0), absorberDep_(0), totalNonIonE_(0), totalRawE_(0),gamFrac_(0), eleFrac_(0), muFrac_(
					0), neutronFrac_(0), hadFrac_(0), muKin_(0), neutronKin_(0), hadKin_(
					0), avgTime_(0), nSiHits_(0), eleShowerSize_(0), hadShowerSize_(0), neutShowerSize_(
					0), muShowerSize_(0),neutronCount_(0),hadCount_(0),muCount_(0),eleCount_(0),gamCount_(0) {

	}
	;

	~HGCSSSamplingSection() {
	}
	;
	//getters
	inline unsigned volNb() const {
		return volNb_;
	}
	;

	inline double volX0trans() const {
		return volX0trans_;
	}
	;

	inline double voldEdx() const {
		return voldEdx_;
	}
	;

	inline double volLambdatrans() const {
		return volLambdatrans_;
	}
	;

	inline double sensDep() const {
		return sensDep_;
	}
	;
	inline double absorberDep() const {
		return absorberDep_;
	}
	;
	inline double totalDep() const {
		return totalRawE_;
	}
	;
	inline double totalNonIonDep() {
		return totalNonIonE_;
	}
	;

	inline double gamDepFrac() const {
		return gamFrac_;
	}
	;
	inline double eleDepFrac() const {
		return eleFrac_;
	}
	;
	inline double muDepFrac() const {
		return muFrac_;
	}
	;

	inline double neutronDepFrac() const {
		return neutronFrac_;
	}
	;

	inline double hadDepFrac() const {
		return hadFrac_;
	}
	;

	inline double neutronKinFlux() const {
		return neutronKin_;
	}
	;
	inline double muKinFlux() const {
		return muKin_;
	}
	;
	inline double hadKinFlux() const {
		return hadKin_;
	}
	;
	inline double eleKinFlux() const {
		return eleKin_;
	}
	;
	inline double gamKinFlux() const {
		return gamKin_;
	}
	;

	inline unsigned muCount() {
		return muCount_;
	}
	;
	inline unsigned hadCount() {
		return hadCount_;
	}
	;
	inline unsigned neutronCount() const {
		return neutronCount_;
	}
	;
	inline unsigned gamCount() {
		return gamCount_;
	}
	;
	inline unsigned eleCount() const {
		return eleCount_;
	}
	;

	inline double avgTime() const {
		return avgTime_;
	}
	;
	inline double eleShowerSize() {
		return eleShowerSize_;
	}
	;
	inline double hadronShowerSize() {
		return hadShowerSize_;
	}
	;
	inline double neutronShowerSize() {
		return neutShowerSize_;
	}
	;
	inline double muShowerSize() {
		return muShowerSize_;
	}
	;
	inline unsigned nSiHits() const {
		return nSiHits_;
	}
	;
	//setters
	inline void volNb(const unsigned & aVal) {
		volNb_ = aVal;
	}
	;
	inline void volX0trans(const double & aVal) {
		volX0trans_ = aVal;
	}
	;
	inline void voldEdx(const double & aVal) {
		voldEdx_ = aVal;
	}
	;
	inline void volLambdatrans(const double & aVal) {
		volLambdatrans_ = aVal;
	}
	;
	inline void sensDep(const double & aVal) {
		sensDep_ = aVal;
	}
	;
	inline void absorberDep(const double & aVal) {
		absorberDep_ = aVal;
	}
	;
	inline void totalDep(const double & aVal) {
		totalRawE_ = aVal;
	}
	;
	inline void totalNonIonDep(const double & aVal) {
		totalNonIonE_ = aVal;
	}
	;
	inline void gamDepFrac(const double & aVal) {
		gamFrac_ = aVal;
	}
	;
	inline void eleDepFrac(const double & aVal) {
		eleFrac_ = aVal;
	}
	;
	inline void muDepFrac(const double & aVal) {
		muFrac_ = aVal;
	}
	;
	inline void neutronDepFrac(const double & aVal) {
		neutronFrac_ = aVal;
	}
	;
	inline void hadDepFrac(const double & aVal) {
		hadFrac_ = aVal;
	}
	;
	inline void avgTime(const double & aVal) {
		avgTime_ = aVal;
	}
	;
	inline void nSiHits(const unsigned & aVal) {
		nSiHits_ = aVal;
	}
	;
	inline void muKinFlux(const double & aVal) {
		muKin_ = aVal;
	}
	;
	inline void neutronKinFlux(const double & aVal) {
		neutronKin_ = aVal;
	}
	;
	inline void hadKinFlux(const double & aVal) {
		hadKin_ = aVal;
	}
	;
	inline void gamKinFlux(const double & aVal) {
		gamKin_ = aVal;
	}
	;
	inline void eleKinFlux(const double & aVal) {
		eleKin_ = aVal;
	}
	;
	inline void eleShowerSize(const double & aVal) {
		eleShowerSize_ = aVal;
	}
	;
	inline void hadronShowerSize(const double & aVal) {
		neutShowerSize_ = aVal;
	}
	;
	inline void neutronShowerSize(const double & aVal) {
		hadShowerSize_ = aVal;
	}
	;
	inline void muShowerSize(const double & aVal) {
		muShowerSize_ = aVal;
	}
	;

	inline void muCount(const double & aVal) {
		muCount_ = aVal;
	}
	;
	inline void hadCount(const double & aVal) {
		hadCount_ = aVal;
	}
	;
	inline void neutronCount(const double & aVal) {
		neutronCount_ = aVal;
	}
	;
	inline void gamCount(const double & aVal) {
		gamCount_ = aVal;
	}
	;
	inline void eleCount(const double & aVal) {
		eleCount_ = aVal;
	}
	;
private:
	unsigned volNb_;
	double volX0trans_;
	double voldEdx_;
	double volLambdatrans_;
	double sensDep_;
	double absorberDep_;
	double totalNonIonE_;
	double totalRawE_;


	double gamFrac_;
	double eleFrac_;
	double muFrac_;
	double neutronFrac_;
	double hadFrac_;


	double neutronKin_;
	double muKin_;
	double hadKin_;
	double eleKin_;
	double gamKin_;


	unsigned muCount_;
	unsigned hadCount_;
	unsigned neutronCount_;
	unsigned gamCount_;
	unsigned eleCount_;


	double eleShowerSize_;
	double hadShowerSize_;
	double neutShowerSize_;
	double muShowerSize_;


	double avgTime_;
	unsigned nSiHits_;

ClassDef(HGCSSSamplingSection,1)
	;

};

typedef std::vector<HGCSSSamplingSection> HGCSSSamplingSectionVec;

#endif

