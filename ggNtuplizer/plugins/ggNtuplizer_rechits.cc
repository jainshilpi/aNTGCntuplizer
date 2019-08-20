#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"

#include "limits"




/////////////////////////////////REMOVE////////////////
/*UShort_t necalSC_;
std::vector<UShort_t> ecalSCindex_;
std::vector<Float_t> ecalSCeta_;
std::vector<Float_t> ecalSCphi_;
std::vector<Float_t> ecalSCEn_;
std::vector<Float_t> ecalSCRawEn_;
std::vector<Float_t> ecalSCetaWidth_;
std::vector<Float_t> ecalSCphiWidth_;
std::vector<Float_t> ecalSC_LICTD_;
std::vector<UChar_t> ecalSC_nL1Spike_;
std::vector<UChar_t> ecalSC_nDiweird_;
std::vector<UChar_t> ecalSC_nWeird_;
std::vector<UChar_t> ecalSC_nSaturated_;
std::vector<UChar_t> ecalSC_nOutOfTime_;
std::vector<UChar_t> ecalSC_nXtals_;
std::vector<Float_t> ecalSC_maxEnXtalTime_;
std::vector<Float_t> ecalSC_maxEnXtalSwissCross_;
std::vector<UChar_t> ecalSC_maxEnXtalBits_;
*/
////////////////////////////////////////////

Int_t nhbheRH_;
Int_t nHBRH_;
Int_t nHERH_;

vector<Float_t> hbherhE_;
vector<Int_t> hbherhDepth_;
vector<Int_t> hbherhiEta_;
vector<Int_t> hbherhiPhi_;
vector<UChar_t> hbherhSubDetID_;
vector<Float_t> hbherhEta_;
vector<Float_t> hbherhPhi_;
vector<Float_t> hbherhX_;
vector<Float_t> hbherhY_;
vector<Float_t> hbherhZ_;

Int_t nhoRH_;
vector<Float_t> horhE_;
vector<Int_t> horhDepth_;
vector<Int_t> horhiEta_;
vector<Int_t> horhiPhi_;
vector<Float_t> horhEta_;
vector<Float_t> horhPhi_;
vector<Float_t> horhX_;
vector<Float_t> horhY_;
vector<Float_t> horhZ_;


Int_t nhfRH_;
vector<Float_t> hfrhE_;
vector<Int_t> hfrhDepth_;
vector<Int_t> hfrhiEta_;
vector<Int_t> hfrhiPhi_;
vector<Float_t> hfrhEta_;
vector<Float_t> hfrhPhi_;
vector<Float_t> hfrhX_;
vector<Float_t> hfrhY_;
vector<Float_t> hfrhZ_;



/////EB
Int_t nebRH_;
vector<Float_t> ebrhE_;
vector<Int_t> ebrhiEta_;
vector<Int_t> ebrhiPhi_;
vector<Int_t> ebrhZside_;
vector<Float_t> ebrhEta_;
vector<Float_t> ebrhPhi_;
vector<Float_t> ebrhX_;
vector<Float_t> ebrhY_;
vector<Float_t> ebrhZ_;

/////EE
Int_t neeRH_;
vector<Float_t> eerhE_;
vector<Int_t> eerhiEta_;
vector<Int_t> eerhiPhi_;
vector<Int_t> eerhZside_;
vector<Float_t> eerhEta_;
vector<Float_t> eerhPhi_;
vector<Float_t> eerhX_;
vector<Float_t> eerhY_;
vector<Float_t> eerhZ_;


/////ES
Int_t nesRH_;
vector<Float_t> esrhE_;
vector<Int_t> esrhiEta_;
vector<Int_t> esrhiPhi_;
vector<Int_t> esrhZside_;
vector<Int_t> esrhPlane_;
vector<Int_t> esrhStrip_;
vector<Float_t> esrhEta_;
vector<Float_t> esrhPhi_;
vector<Float_t> esrhX_;
vector<Float_t> esrhY_;
vector<Float_t> esrhZ_;

///CSC
Int_t nCSCSeg_;
vector<int> cscSegEE_;
vector<int> cscSegStation_;
vector<int> cscSegRing_;
vector<Float_t> cscSegDx_;
vector<Float_t> cscSegDy_;
vector<Float_t> cscSegDz_;
vector<Float_t> cscSegX_;
vector<Float_t> cscSegY_;
vector<Float_t> cscSegZ_;
vector<Float_t> cscSegEta_;
vector<Float_t> cscSegPhi_;
vector<int> cscSegDim_;
vector<int> cscSegnRH_;

/*
Float_t ggNtuplizer::ECALrecHitE(const DetId & id, const EcalRecHitCollection *recHits, int di, int dj){
	DetId nid;
	if(di == 0 && dj == 0) nid = id;
	else if(id.subdetId() == EcalBarrel) nid = EBDetId::offsetBy(id, di, dj);
	else if(id.subdetId() == EcalEndcap) nid = EEDetId::offsetBy(id, di, dj);
	return ((recHits->find(nid) != recHits->end()) ? recHits->find(nid)->energy() : 0.);
};



Float_t ggNtuplizer::swissCross(const DetId & id, noZS::EcalClusterLazyTools & ltNoZS){
	const EcalRecHitCollection * recHits = nullptr;
	if (id.subdetId() == EcalBarrel) {
		recHits = ltNoZS.getEcalEBRecHitCollection();
		EBDetId ebId(id);
	    // avoid recHits at |iEta|=85 where one side of the neighbours is missing
		if (std::abs(ebId.ieta())==85) return -100.;
		Float_t e1 = ECALrecHitE(id, recHits);
		Float_t s4 = 0.;
		if (!(e1 > 0.)) return -900.;
		s4 += ECALrecHitE(id, recHits,  1,  0);
		s4 += ECALrecHitE(id, recHits, -1,  0);
		s4 += ECALrecHitE(id, recHits,  0,  1);
		s4 += ECALrecHitE(id, recHits,  0, -1);
		return (1 - s4 / e1);
	} else if (id.subdetId() == EcalEndcap) {
		recHits = ltNoZS.getEcalEERecHitCollection();
		EEDetId eeId(id);
		Float_t e1 = ECALrecHitE(id, recHits);
		Float_t s4 = 0.;
		if (!(e1 > 0.)) return -920.;
		s4 += ECALrecHitE(id, recHits,  1,  0);
		s4 += ECALrecHitE(id, recHits, -1,  0);
		s4 += ECALrecHitE(id, recHits,  0,  1);
		s4 += ECALrecHitE(id, recHits,  0, -1);
		return (1 - s4 / e1);
	}
	return -9999.;
};
*/


void ggNtuplizer::branchesRecHit(TTree* tree) {
  /*tree->Branch("necalSC",                    		&necalSC_);
	tree->Branch("ecalSCindex",                    	&ecalSCindex_);
	tree->Branch("ecalSCeta",                    	&ecalSCeta_);
	tree->Branch("ecalSCphi",                    	&ecalSCphi_);
	tree->Branch("ecalSCEn",                    	&ecalSCEn_);
	tree->Branch("ecalSCRawEn",                    	&ecalSCRawEn_);
	tree->Branch("ecalSCetaWidth",                  &ecalSCetaWidth_);
	tree->Branch("ecalSCphiWidth",                  &ecalSCphiWidth_);
	tree->Branch("ecalSC_LICTD",                    &ecalSC_LICTD_);
	tree->Branch("ecalSC_nL1Spike",                 &ecalSC_nL1Spike_);
	tree->Branch("ecalSC_nDiweird",                 &ecalSC_nDiweird_);
	tree->Branch("ecalSC_nWeird",					&ecalSC_nWeird_);
	tree->Branch("ecalSC_nSaturated",				&ecalSC_nSaturated_);
	tree->Branch("ecalSC_nOutOfTime",				&ecalSC_nOutOfTime_);
	tree->Branch("ecalSC_nXtals",					&ecalSC_nXtals_);
	tree->Branch("ecalSC_maxEnXtalTime",			&ecalSC_maxEnXtalTime_);
	tree->Branch("ecalSC_maxEnXtalSwissCross",		&ecalSC_maxEnXtalSwissCross_);
	tree->Branch("ecalSC_maxEnXtalBits",			&ecalSC_maxEnXtalBits_);
  */

  ///HBHE
  tree->Branch("nhbheRH",         &nhbheRH_);
  tree->Branch("nHBRH",         &nHBRH_);
  tree->Branch("nHERH",         &nHERH_);
  tree->Branch("hbherhE",         &hbherhE_);
  tree->Branch("hbherhDepth",         &hbherhDepth_);
  tree->Branch("hbherhiEta",         &hbherhiEta_);
  tree->Branch("hbherhiPhi",         &hbherhiPhi_);
  tree->Branch("hbherhSubDetID",         &hbherhSubDetID_);
  
  tree->Branch("hbherhEta",         &hbherhEta_);
  tree->Branch("hbherhPhi",         &hbherhPhi_);
  tree->Branch("hbherhX",         &hbherhX_);
  tree->Branch("hbherhY",         &hbherhY_);
  tree->Branch("hbherhZ",         &hbherhZ_);

  ///HO
  tree->Branch("nhoRH",         &nhoRH_);
  tree->Branch("horhE",         &horhE_);
  tree->Branch("horhDepth",         &horhDepth_);
  tree->Branch("horhiEta",         &horhiEta_);
  tree->Branch("horhiPhi",         &horhiPhi_);

  tree->Branch("horhEta",         &horhEta_);
  tree->Branch("horhPhi",         &horhPhi_);
  tree->Branch("horhX",         &horhX_);
  tree->Branch("horhY",         &horhY_);
  tree->Branch("horhZ",         &horhZ_);


  ///HF
  tree->Branch("nhfRH",         &nhfRH_);
  tree->Branch("hfrhE",         &hfrhE_);
  tree->Branch("hfrhDepth",         &hfrhDepth_);
  tree->Branch("hfrhiEta",         &hfrhiEta_);
  tree->Branch("hfrhiPhi",         &hfrhiPhi_);

  tree->Branch("hfrhEta",         &hfrhEta_);
  tree->Branch("hfrhPhi",         &hfrhPhi_);
  tree->Branch("hfrhX",         &hfrhX_);
  tree->Branch("hfrhY",         &hfrhY_);
  tree->Branch("hfrhZ",         &hfrhZ_);

  ///EB
  tree->Branch("nebRH",         &nebRH_);
  tree->Branch("ebrhE",         &ebrhE_);
  tree->Branch("ebrhiEta",         &ebrhiEta_);
  tree->Branch("ebrhiPhi",         &ebrhiPhi_);
  tree->Branch("ebrhZside",         &ebrhZside_);
  tree->Branch("ebrhEta",         &ebrhEta_);
  tree->Branch("ebrhPhi",         &ebrhPhi_);
  tree->Branch("ebrhX",         &ebrhX_);
  tree->Branch("ebrhY",         &ebrhY_);
  tree->Branch("ebrhZ",         &ebrhZ_);

  ///EE
  tree->Branch("neeRH",         &neeRH_);
  tree->Branch("eerhE",         &eerhE_);
  tree->Branch("eerhiEta",         &eerhiEta_);
  tree->Branch("eerhiPhi",         &eerhiPhi_);
  tree->Branch("eerhZside",         &eerhZside_);
  tree->Branch("eerhEta",         &eerhEta_);
  tree->Branch("eerhPhi",         &eerhPhi_);
  tree->Branch("eerhX",         &eerhX_);
  tree->Branch("eerhY",         &eerhY_);
  tree->Branch("eerhZ",         &eerhZ_);

  ///ES
  tree->Branch("nesRH",         &nesRH_);
  tree->Branch("esrhE",         &esrhE_);
  tree->Branch("esrhiEta",         &esrhiEta_);
  tree->Branch("esrhiPhi",         &esrhiPhi_);
  tree->Branch("esrhZside",         &esrhZside_);
  tree->Branch("esrhPlane",         &esrhPlane_);
  tree->Branch("esrhStrip",         &esrhStrip_);
  tree->Branch("esrhEta",         &esrhEta_);
  tree->Branch("esrhPhi",         &esrhPhi_);
  tree->Branch("esrhX",         &esrhX_);
  tree->Branch("esrhY",         &esrhY_);
  tree->Branch("esrhZ",         &esrhZ_);

  ///CSC
  tree->Branch("nCSCSeg",         &nCSCSeg_);
  tree->Branch("cscSegEE",         &cscSegEE_);
  tree->Branch("cscSegStation",         &cscSegStation_);
  tree->Branch("cscSegRing",         &cscSegRing_);
  tree->Branch("cscSegDx",         &cscSegDx_);
  tree->Branch("cscSegDy",         &cscSegDy_);
  tree->Branch("cscSegDz",         &cscSegDz_);
  tree->Branch("cscSegX",         &cscSegX_);
  tree->Branch("cscSegY",         &cscSegY_);
  tree->Branch("cscSegZ",         &cscSegZ_);
  tree->Branch("cscSegDim",         &cscSegDim_);
  tree->Branch("cscSegnRH",         &cscSegnRH_);
  tree->Branch("cscSegEta",         &cscSegEta_);
  tree->Branch("cscSegPhi",         &cscSegPhi_);



};


/*
Float_t ggNtuplizer::getLICTD(const reco::SuperCluster *sc, noZS::EcalClusterLazyTools & ltNoZS, Float_t &_maxEnXtalTime, UChar_t & _nL1Spike, UChar_t & _nDiweird, UChar_t & _nWeird, UChar_t & _nSaturated, UChar_t & _nOutOfTime, UChar_t & _nXtals, UChar_t & _maxEnXtalBits, Float_t & _maxEnXtalSwissCross){

	_maxEnXtalTime = -999.;
	_nL1Spike = 0;
	_nDiweird = 0;
	_nWeird = 0;
	_nSaturated = 0;
	_nOutOfTime = 0;
	_nXtals = 0;
	_maxEnXtalBits= 0;
	_maxEnXtalSwissCross = -999.;

	if(sc->clustersSize()<1) return -999999.;
	if(!sc->clusters().isAvailable()) return -99999.;

	DetId seedDetID = sc->seed()->seed();
	Bool_t isEB = (seedDetID.subdetId() == EcalBarrel);
	if(!isEB && !(seedDetID.subdetId() == EcalEndcap)) return -999.;

	const EcalRecHitCollection *recHitsEB = ltNoZS.getEcalEBRecHitCollection();
	const EcalRecHitCollection *recHitsEE = ltNoZS.getEcalEERecHitCollection();

	std::vector<uint32_t> _detIDs_nL1Spike;
	std::vector<uint32_t> _detIDs_nDiweird;
	std::vector<uint32_t> _detIDs_nWeird;
	std::vector<uint32_t> _detIDs_nSaturated;
	std::vector<uint32_t> _detIDs_nOutOfTime;
	std::vector<uint32_t> _detIDs_nXtals;

	const EcalRecHit *maxEnergyXtalRecHit_ 	= nullptr;
	Float_t maxEnergyXtal_ 					= std::numeric_limits<Float_t>::min();
	Float_t maxTime_						= std::numeric_limits<Float_t>::min();
	Float_t minTime_						= std::numeric_limits<Float_t>::max();
	Bool_t  _timeValid						= 0;

	for(reco::CaloCluster_iterator cluster = sc->clustersBegin(); cluster != sc->clustersEnd(); cluster++){
		if((*cluster)->size()<1) continue;

		for(const std::pair<DetId, Float_t> & _xtal : (*cluster)->hitsAndFractions()){

			// Find rechit
			const EcalRecHit * _xtalHit = nullptr;
			if(isEB){
				if(recHitsEB->find(_xtal.first) != recHitsEB->end()) _xtalHit = &(*(recHitsEB->find (_xtal.first)));
				else if(recHitsEE->find(_xtal.first) != recHitsEE->end()) _xtalHit = &(*(recHitsEE->find (_xtal.first)));
			} else{
				if(recHitsEE->find(_xtal.first) != recHitsEE->end()) _xtalHit = &(*(recHitsEE->find (_xtal.first)));
				else if(recHitsEB->find(_xtal.first) != recHitsEB->end()) _xtalHit = &(*(recHitsEB->find (_xtal.first)));
			}

			if(!_xtalHit) continue;

			// skip xtal if energy deposit is < 1 GeV
			if(_xtalHit->energy() > 1.) {
				// Get time range
				if(_xtalHit->time() > maxTime_) maxTime_ = _xtalHit->time();
				if(_xtalHit->time() < minTime_) minTime_ = _xtalHit->time();
				_timeValid = 1.;
			}

			// Get max energy xtal
			if(_xtalHit->energy() > maxEnergyXtal_){
				maxEnergyXtal_ = _xtalHit->energy();
				maxEnergyXtalRecHit_ = _xtalHit;
			}

			uint32_t _xtalDetID = _xtal.first.rawId();
			if(_xtalHit->checkFlag(EcalRecHit::kL1SpikeFlag) && (std::find(_detIDs_nL1Spike.begin(), _detIDs_nL1Spike.end(), _xtalDetID) == _detIDs_nL1Spike.end())){
				_detIDs_nL1Spike.push_back(_xtalDetID);
				_nL1Spike++;
			}
			if(_xtalHit->checkFlag(EcalRecHit::kDiWeird) && (std::find(_detIDs_nDiweird.begin(), _detIDs_nDiweird.end(), _xtalDetID) == _detIDs_nDiweird.end())){
				_detIDs_nDiweird.push_back(_xtalDetID);
				_nDiweird++;
			}
			if(_xtalHit->checkFlag(EcalRecHit::kWeird) && (std::find(_detIDs_nWeird.begin(), _detIDs_nWeird.end(), _xtalDetID) == _detIDs_nWeird.end())){
				_detIDs_nWeird.push_back(_xtalDetID);
				_nWeird++;
			}
			if(_xtalHit->checkFlag(EcalRecHit::kSaturated) && (std::find(_detIDs_nSaturated.begin(), _detIDs_nSaturated.end(), _xtalDetID) == _detIDs_nSaturated.end())){
				_detIDs_nSaturated.push_back(_xtalDetID);
				_nSaturated++;
			}
			if(_xtalHit->checkFlag(EcalRecHit::kOutOfTime) && (std::find(_detIDs_nOutOfTime.begin(), _detIDs_nOutOfTime.end(), _xtalDetID) == _detIDs_nOutOfTime.end())){
				_detIDs_nOutOfTime.push_back(_xtalDetID);
				_nOutOfTime++;
			}
			if(std::find(_detIDs_nXtals.begin(), _detIDs_nXtals.end(), _xtalDetID) == _detIDs_nXtals.end()){
				_detIDs_nXtals.push_back(_xtalDetID);
				_nXtals++;
			}
		}
	};

	if(!maxEnergyXtalRecHit_) return -9999.;
	_maxEnXtalTime 		= maxEnergyXtalRecHit_->time();
	_maxEnXtalBits 		= 0;
	if(maxEnergyXtalRecHit_->checkFlag(EcalRecHit::kDiWeird)) setbit(_maxEnXtalBits, 0);
	if(maxEnergyXtalRecHit_->checkFlag(EcalRecHit::kWeird)) setbit(_maxEnXtalBits, 1);
	if(maxEnergyXtalRecHit_->checkFlag(EcalRecHit::kSaturated)) setbit(_maxEnXtalBits, 2);
	if(maxEnergyXtalRecHit_->checkFlag(EcalRecHit::kOutOfTime)) setbit(_maxEnXtalBits, 3);

	_maxEnXtalSwissCross = swissCross(maxEnergyXtalRecHit_->detid(), ltNoZS);
	if(maxEnergyXtalRecHit_->detid().rawId() != seedDetID.rawId()) setbit(_maxEnXtalBits, 4);

	Float_t LICTD_ = (_timeValid) ? std::max(std::abs(maxTime_ - _maxEnXtalTime), std::abs(_maxEnXtalTime - minTime_)) : -9999999.;

	return LICTD_;
};
*/

///SJ
void ggNtuplizer::fillRecHits(const edm::Event& e, const edm::EventSetup& es){

  //bool debugRH = true;
  bool debugRH = false;
  ///HCAL rechits: https://cmssdt.cern.ch/lxr/source/CalibCalorimetry/CaloMiscalibTools/plugins/HcalRecHitRecalib.cc#0063
  //http://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_9_4_14/doc/html/d9/d45/HcalHBHEMuonAnalyzer_8cc_source.html#l00237
  //http://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_9_4_14/doc/html/d2/d7a/HCALRecHitAnalyzer_8cc_source.html#l00452 ---> mainly this one
  //https://cmssdt.cern.ch/lxr/source/Calibration/IsolatedParticles/src/FindDistCone.cc?v=CMSSW_10_2_0#0193

  if(debugRH){

    std::cout<<"Getting Calogeom"<<std::endl;
  }


  es.get<CaloGeometryRecord>().get(pG_);
  const CaloGeometry* geo = pG_.product();

  if(debugRH){

    std::cout<<"Got Calogeom"<<std::endl;
  }

  //HBHE
  nhbheRH_ = 0;  
  nHBRH_ = 0;
  nHERH_ = 0;

  hbherhE_.clear();                 
  hbherhDepth_.clear();
  hbherhiEta_.clear();
  hbherhiPhi_.clear();
  hbherhSubDetID_.clear();

  hbherhEta_.clear();
  hbherhPhi_.clear();
  hbherhX_.clear();
  hbherhY_.clear();
  hbherhZ_.clear();

  //HO
  nhoRH_ = 0;
  horhE_.clear();                 
  horhDepth_.clear();
  horhiEta_.clear();
  horhiPhi_.clear();


  horhEta_.clear();
  horhPhi_.clear();
  horhX_.clear();
  horhY_.clear();
  horhZ_.clear();

  ///HF
  nhfRH_ = 0;
  hfrhE_.clear();                 
  hfrhDepth_.clear();
  hfrhiEta_.clear();
  hfrhiPhi_.clear();


  hfrhEta_.clear();
  hfrhPhi_.clear();
  hfrhX_.clear();
  hfrhY_.clear();
  hfrhZ_.clear();


  //EB
  nebRH_ = 0;  

  ebrhE_.clear();                 
  ebrhiEta_.clear();
  ebrhiPhi_.clear();
  ebrhZside_.clear();

  ebrhEta_.clear();
  ebrhPhi_.clear();
  ebrhX_.clear();
  ebrhY_.clear();
  ebrhZ_.clear();

  //EE
  neeRH_ = 0;  

  eerhE_.clear();                 
  eerhiEta_.clear();
  eerhiPhi_.clear();
  esrhZside_.clear();

  eerhEta_.clear();
  eerhPhi_.clear();
  eerhX_.clear();
  eerhY_.clear();
  eerhZ_.clear();

  //ES
  nesRH_ = 0;  

  esrhE_.clear();                 
  esrhiEta_.clear();
  esrhiPhi_.clear();
  esrhZside_.clear();
  esrhPlane_.clear();
  esrhStrip_.clear();

  esrhEta_.clear();
  esrhPhi_.clear();
  esrhX_.clear();
  esrhY_.clear();
  esrhZ_.clear();
  

  ///CSC
  nCSCSeg_ = 0;
  cscSegEE_.clear();
  cscSegStation_.clear();
  cscSegRing_.clear();
  cscSegDx_.clear();
  cscSegDy_.clear();
  cscSegDz_.clear();
  cscSegX_.clear();
  cscSegY_.clear();
  cscSegZ_.clear();
  cscSegDim_.clear();
  cscSegnRH_.clear();
  cscSegEta_.clear();
  cscSegPhi_.clear();

  if(debugRH){

    std::cout<<"Getting HB hangles"<<std::endl;
  }


  edm::Handle<HBHERecHitCollection> HBHERecHitsHandle;
  edm::Handle<HFRecHitCollection> HFRecHitsHandle;
  edm::Handle<HORecHitCollection> HORecHitsHandle;

  const HBHERecHitCollection* HBHERecHits = nullptr;
  const HFRecHitCollection* HFRecHits = nullptr;
  const HORecHitCollection* HORecHits = nullptr;

  e.getByToken(hbheRecHitCollection_, HBHERecHitsHandle);
  if (!HBHERecHitsHandle.isValid()) {
    LogDebug("") << "ggNtuplizer_rechits: Error! HBHE can't get product!" << std::endl;
  } else {
    HBHERecHits = HBHERecHitsHandle.product();  // get a ptr to the product
  }

  e.getByToken(hoRecHitCollection_, HORecHitsHandle);
  if (!HORecHitsHandle.isValid()) {
    LogDebug("") << "ggNtuplizer_rechits: Error! HO can't get product!" << std::endl;
  } else {
    HORecHits = HORecHitsHandle.product();  // get a ptr to the product
  }

  e.getByToken(hfRecHitCollection_, HFRecHitsHandle);
  if (!HFRecHitsHandle.isValid()) {
    LogDebug("") << "ggNtuplizer_rechits: Error! HF can't get product!" << std::endl;
  } else {
    HFRecHits = HFRecHitsHandle.product();  // get a ptr to the product
  }



  if(debugRH){

    std::cout<<"Got HBHandles"<<std::endl;
  }


  // Loop over HBHERecHit's
  HBHERecHitCollection::const_iterator hbherechit;
    
  for (hbherechit = HBHERecHits->begin(); hbherechit != HBHERecHits->end(); hbherechit++) {

    HcalDetId det = hbherechit->id();
    double Energy = hbherechit->energy();
    Int_t depth = det.depth();
    Int_t ieta = det.ieta();
    Int_t iphi = det.iphi();

    hbherhE_.push_back(Energy);
    hbherhDepth_.push_back(depth);
    hbherhiEta_.push_back(ieta);
    hbherhiPhi_.push_back(iphi);

    const GlobalPoint & rechitPoint = geo->getPosition(det);
    hbherhEta_.push_back(rechitPoint.eta());
    hbherhPhi_.push_back(rechitPoint.phi());
    hbherhX_.push_back(rechitPoint.x());
    hbherhY_.push_back(rechitPoint.y());
    hbherhZ_.push_back(rechitPoint.z());

    /*double eta = hHCAL_ieta_iphi_etaMap->getBinContent(EtaRing+1,iphi);
    double phi = hHCAL_ieta_iphi_phiMap->getBinContent(EtaRing+1,iphi);
    double theta = 2*TMath::ATan(exp(-1*eta));
    double ET = Energy*TMath::Sin(theta);
    */
    HcalSubdetector HcalNum = det.subdet();

    UShort_t tmpdetIDbit = 0;
    if(HcalNum==HcalEmpty) setbit(tmpdetIDbit, 0);
    if(HcalNum==HcalBarrel) setbit(tmpdetIDbit, 1);
    if(HcalNum==HcalEndcap) setbit(tmpdetIDbit, 2);
    if(HcalNum==HcalOuter) setbit(tmpdetIDbit, 3);
    if(HcalNum==HcalForward) setbit(tmpdetIDbit, 4);
    if(HcalNum==HcalTriggerTower) setbit(tmpdetIDbit, 5);
    if(HcalNum==HcalOther) setbit(tmpdetIDbit, 6);

    hbherhSubDetID_.push_back(tmpdetIDbit);
    if (HcalNum == HcalBarrel) {
      nHBRH_++;
    } else { // HcalEndcap
      nHERH_++;
    }
    
    nhbheRH_++;
  }//for (hbherechit = HBHERecHits->begin(); hbherechit != HBHERecHits->end(); hbherechit++)


  if(debugRH){

    std::cout<<"Looping HO"<<std::endl;
  }

  
  ///HO rechit collection
  HORecHitCollection::const_iterator horechit;

  for (horechit = HORecHits->begin(); horechit != HORecHits->end(); horechit++) {


    HcalDetId det = horechit->id();
    double Energy = horechit->energy();
    Int_t depth = det.depth();
    Int_t ieta = det.ieta();
    Int_t iphi = det.iphi();

    horhE_.push_back(Energy);
    horhDepth_.push_back(depth);
    horhiEta_.push_back(ieta);
    horhiPhi_.push_back(iphi);

    const GlobalPoint & rechitPoint = geo->getPosition(det);
    horhEta_.push_back(rechitPoint.eta());
    horhPhi_.push_back(rechitPoint.phi());
    horhX_.push_back(rechitPoint.x());
    horhY_.push_back(rechitPoint.y());
    horhZ_.push_back(rechitPoint.z());

    nhoRH_++;
  }


  if(debugRH){

    std::cout<<"Looping HF"<<std::endl;
  }

  //HF rechit collection
  HFRecHitCollection::const_iterator hfrechit;

  for (hfrechit = HFRecHits->begin(); hfrechit != HFRecHits->end(); hfrechit++) {

    HcalDetId det = hfrechit->id();
    double Energy = hfrechit->energy();
    Int_t depth = det.depth();
    Int_t ieta = det.ieta();
    Int_t iphi = det.iphi();

    hfrhE_.push_back(Energy);
    hfrhDepth_.push_back(depth);
    hfrhiEta_.push_back(ieta);
    hfrhiPhi_.push_back(iphi);

    const GlobalPoint & rechitPoint = geo->getPosition(det);
    hfrhEta_.push_back(rechitPoint.eta());
    hfrhPhi_.push_back(rechitPoint.phi());
    hfrhX_.push_back(rechitPoint.x());
    hfrhY_.push_back(rechitPoint.y());
    hfrhZ_.push_back(rechitPoint.z());

    nhfRH_++;

  }


  //////////////////////////////ECAL rechits/////////////////////////////////////////
  if(debugRH){

    std::cout<<"Getting ECAL RH"<<std::endl;
  }

  edm::Handle<EcalRecHitCollection> barrelRecHitsHandle;
  edm::Handle<EcalRecHitCollection> endcapRecHitsHandle;
  edm::Handle<EcalRecHitCollection> esRecHitsHandle;
  e.getByToken(ebReducedRecHitCollection_,barrelRecHitsHandle);
  e.getByToken(eeReducedRecHitCollection_,endcapRecHitsHandle);
  e.getByToken(esReducedRecHitCollection_,esRecHitsHandle);

  const EcalRecHitCollection* EBRecHits = nullptr;
  const EcalRecHitCollection* EERecHits = nullptr;
  const EcalRecHitCollection* ESRecHits = nullptr;
  
  if ( !barrelRecHitsHandle.isValid() ){
    LogDebug("") << "ggNtuplizer_rechits: Error! EB rechits can't get product!" << std::endl;
  } else{
    EBRecHits = barrelRecHitsHandle.product();
  }

  if ( !endcapRecHitsHandle.isValid() ){
    LogDebug("") << "ggNtuplizer_rechits: Error! EE rechits can't get product!" << std::endl;
  } else{
    EERecHits = endcapRecHitsHandle.product();
  }
  
  if ( !esRecHitsHandle.isValid() ){
    LogDebug("") << "ggNtuplizer_rechits: Error! ES rechits can't get product!" << std::endl;
  } else{
    ESRecHits = esRecHitsHandle.product();
  }

  if(debugRH){

    std::cout<<"Looping EB"<<std::endl;
  }

  
  //EB
  EcalRecHitCollection::const_iterator ebrechit;

  if(debugRH){
    std::cout<<"Start looping EB now"<<std::endl;
  }

  if(debugRH) std::cout<<"EB EH size "<<EBRecHits->size()<<std::endl;

  for( ebrechit = EBRecHits->begin(); ebrechit != EBRecHits->end(); ebrechit++ ){

    if(debugRH){
      std::cout<<"1st EB rechit"<<std::endl;
    }

    double Energy = ebrechit->energy();

    if(debugRH){
      std::cout<<"E "<<Energy<<std::endl;
    }

    EBDetId det = ebrechit->id();

    if(debugRH){
      std::cout<<"detID "<<det<<std::endl;
    }

    int ieta = det.ieta();
    if(debugRH){
      std::cout<<"ieta "<<ieta<<std::endl;
    }


    int iphi = det.iphi();
    int zside = det.zside();

    ebrhE_.push_back(Energy);
    ebrhiEta_.push_back(ieta);
    ebrhiPhi_.push_back(iphi);
    ebrhZside_.push_back(zside);

    const GlobalPoint & rechitPoint = geo->getPosition(det);
    ebrhEta_.push_back(rechitPoint.eta());

    if(debugRH){
      std::cout<<"eta "<<rechitPoint.eta()<<std::endl;
    }


    ebrhPhi_.push_back(rechitPoint.phi());
    ebrhX_.push_back(rechitPoint.x());
    ebrhY_.push_back(rechitPoint.y());
    ebrhZ_.push_back(rechitPoint.z());

    nebRH_++;
  }


  if(debugRH){

    std::cout<<"Looping EE"<<std::endl;
  }

  ///EE
  EcalRecHitCollection::const_iterator eerechit;
  for( eerechit = EERecHits->begin(); eerechit != EERecHits->end(); eerechit++ ){

    double Energy = eerechit->energy();
    EEDetId det = eerechit->id();
    int ieta = det.ix();
    int iphi = det.iy();
    int zside = det.zside();

    eerhE_.push_back(Energy);
    eerhiEta_.push_back(ieta);
    eerhiPhi_.push_back(iphi);
    eerhZside_.push_back(zside);

    const GlobalPoint & rechitPoint = geo->getPosition(det);
    eerhEta_.push_back(rechitPoint.eta());
    eerhPhi_.push_back(rechitPoint.phi());
    eerhX_.push_back(rechitPoint.x());
    eerhY_.push_back(rechitPoint.y());
    eerhZ_.push_back(rechitPoint.z());

    neeRH_++;
  }


  ///ES
  EcalRecHitCollection::const_iterator esrechit;
  for( esrechit = ESRecHits->begin(); esrechit != ESRecHits->end(); esrechit++ ){

    double Energy = esrechit->energy();
    ESDetId det = esrechit->id();
    int ieta = det.six();
    int iphi = det.siy();
    int zside = det.zside();
    int strip = det.strip();

    int plane = det.plane();

    esrhE_.push_back(Energy);
    esrhiEta_.push_back(ieta);
    esrhiPhi_.push_back(iphi);
    esrhPlane_.push_back(plane);
    esrhStrip_.push_back(strip);
    esrhZside_.push_back(zside);

    const GlobalPoint & rechitPoint = geo->getPosition(det);
    esrhEta_.push_back(rechitPoint.eta());
    esrhPhi_.push_back(rechitPoint.phi());
    esrhX_.push_back(rechitPoint.x());
    esrhY_.push_back(rechitPoint.y());
    esrhZ_.push_back(rechitPoint.z());

    nesRH_++;
  }


  ///////////////////////////////////MUON rechits

  edm::Handle<CSCSegmentCollection> allCSCSegments;
  e.getByToken(cscSegmentsCollection_, allCSCSegments);
  if(allCSCSegments.isValid()){ 
    if(allCSCSegments->size()>0){
      
      LogDebug("rpcefficiency")<<"CSC \t Number of CSC Segments in this event = "<<allCSCSegments->size();
      
      CSCSegmentCollection::const_iterator segment;
      
      for (segment = allCSCSegments->begin();segment!=allCSCSegments->end(); ++segment){
	CSCDetId CSCId = segment->cscDetId();
	int cscEndCap = CSCId.endcap();
	int cscStation = CSCId.station();
	int cscRing = CSCId.ring();

	LocalPoint segmentPosition= segment->localPosition();
	LocalVector segmentDirection=segment->localDirection();
	float dz=segmentDirection.z();
	float Xo=segmentPosition.x();
	float Yo=segmentPosition.y();
	float Zo=segmentPosition.z();
	float dx=segmentDirection.x();
	float dy=segmentDirection.y();

	float eta = segmentPosition.eta();
	float phi = segmentPosition.phi();
	
	int segDim = segment->dimension();
	int segnRH = segment->nRecHits();

	cscSegEE_.push_back(cscEndCap);
	cscSegStation_.push_back(cscStation);
	cscSegRing_.push_back(cscRing);
	cscSegDx_.push_back(dx);
	cscSegDy_.push_back(dy);
	cscSegDz_.push_back(dz);
	cscSegX_.push_back(Xo);
	cscSegY_.push_back(Yo);
	cscSegZ_.push_back(Zo);
	cscSegDim_.push_back(segDim);
	cscSegnRH_.push_back(segnRH);
	cscSegEta_.push_back(eta);
	cscSegPhi_.push_back(phi);

	nCSCSeg_++;
      }//for (segment = allCSCSegments->begin();segment!=allCSCSegments->end(); ++segment)
      
	///////////////////////////////////////////////
	/*
      std::map<CSCDetId,int> CSCSegmentsCounter;
      CSCSegmentCollection::const_iterator segment;
      int segmentsInThisEventInTheEndcap=0;
      for (segment = allCSCSegments->begin();segment!=allCSCSegments->end(); ++segment){
	CSCSegmentsCounter[segment->cscDetId()]++;
	
      }    


      statistics->Fill(allCSCSegments->size()+18);
      LogDebug("rpcefficiency")<<"CSC \t loop over all the CSCSegments ";
      for (segment = allCSCSegments->begin();segment!=allCSCSegments->end(); ++segment){
	CSCDetId CSCId = segment->cscDetId();
	if(CSCSegmentsCounter[CSCId]==1 && CSCId.ring()!=1 && allCSCSegments->size()>=2){
	  LogDebug("rpcefficiency")<<"CSC \t \t yes";
          int cscEndCap = CSCId.endcap();
	  int cscStation = CSCId.station();
	  int cscRing = CSCId.ring();
	  int rpcRegion = 1; if(cscEndCap==2) rpcRegion= -1;//Relation among the endcaps
	  int rpcRing = cscRing;
	  if(cscRing==4)rpcRing =1;
	  int rpcStation = cscStation;
	  int rpcSegment = CSCId.chamber();
	  LocalPoint segmentPosition= segment->localPosition();
	  LocalVector segmentDirection=segment->localDirection();
	  float dz=segmentDirection.z();
	  LogDebug("rpcefficiency")<<"CSC \t \t Is a good Segment? dim = 4, 4 <= nRecHits <= 10 Incident angle int range 45 < "<<acos(dz)*180/3.1415926<<" < 135? ";
	  if(segment->dimension()==4 && (segment->nRecHits()<=10 && segment->nRecHits()>=4)&& acos(dz)*180/3.1415926 > 45. && acos(dz)*180/3.1415926 < 160. ){ 
	    float Xo=segmentPosition.x();
	    float Yo=segmentPosition.y();
	    float Zo=segmentPosition.z();
	    float dx=segmentDirection.x();
	    float dy=segmentDirection.y();
	    float dz=segmentDirection.z();
	    LogDebug("rpcefficiency")<<"CSC \t \t Getting chamber from Geometry";
	    const CSCChamber* TheChamber=cscGeo->chamber(CSCId); 
	    LogDebug("rpcefficiency")<<"CSC \t \t Getting ID from Chamber";
	    const CSCDetId TheId=TheChamber->id();
	    LogDebug("rpcefficiency")<<"CSC \t \t Printing The Id"<<TheId;
	    std::set<RPCDetId> rollsForThisCSC = rollstoreCSC[CSCStationIndex(rpcRegion,rpcStation,rpcRing,rpcSegment)];
	    LogDebug("rpcefficiency")<<"CSC \t \t Number of rolls for this CSC = "<<rollsForThisCSC.size();
	    if(rpcRing!=1){
	      //Loop over all the rolls
	      for (std::set<RPCDetId>::iterator iteraRoll = rollsForThisCSC.begin();iteraRoll != rollsForThisCSC.end(); iteraRoll++){
		const RPCRoll* rollasociated = rpcGeo->roll(*iteraRoll);
		RPCDetId rpcId = rollasociated->id();
		const BoundPlane & RPCSurface = rollasociated->surface(); 
		GlobalPoint CenterPointRollGlobal = RPCSurface.toGlobal(LocalPoint(0,0,0));
		LocalPoint CenterRollinCSCFrame = TheChamber->toLocal(CenterPointRollGlobal);
		float D=CenterRollinCSCFrame.z();
		float X=Xo+dx*D/dz;
		float Y=Yo+dy*D/dz;
		float Z=D;
		const TrapezoidalStripTopology* top_=dynamic_cast<const TrapezoidalStripTopology*>(&(rollasociated->topology()));
		LocalPoint xmin = top_->localPosition(0.);
		LogDebug("rpcefficiency")<<"CSC \t \t \t xmin of this  Roll "<<xmin<<"cm";
		LocalPoint xmax = top_->localPosition((float)rollasociated->nstrips());
		LogDebug("rpcefficiency")<<"CSC \t \t \t xmax of this  Roll "<<xmax<<"cm";
		float rsize = fabs( xmax.x()-xmin.x() );
		LogDebug("rpcefficiency")<<"CSC \t \t \t Roll Size "<<rsize<<"cm";
		float stripl = top_->stripLength();
		float stripw = top_->pitch();
		float extrapolatedDistance = sqrt((X-Xo)*(X-Xo)+(Y-Yo)*(Y-Yo)+(Z-Zo)*(Z-Zo));
		if(extrapolatedDistance<=MaxD){ 
		  GlobalPoint GlobalPointExtrapolated=TheChamber->toGlobal(LocalPoint(X,Y,Z));
		  LocalPoint PointExtrapolatedRPCFrame = RPCSurface.toLocal(GlobalPointExtrapolated);
		  if(fabs(PointExtrapolatedRPCFrame.z()) < 10. && 
		     fabs(PointExtrapolatedRPCFrame.x()) < rsize*0.5 && 
		     fabs(PointExtrapolatedRPCFrame.y()) < stripl*0.5){ 
		    RPCDetId  rollId = rollasociated->id();
		    RPCGeomServ rpcsrv(rollId);
		    std::string nameRoll = rpcsrv.name();
		    LogDebug("rpcefficiency")<<"CSC \t \t \t \t The RPCName is "<<nameRoll;
		    const float stripPredicted = 
		      rollasociated->strip(LocalPoint(PointExtrapolatedRPCFrame.x(),PointExtrapolatedRPCFrame.y(),0.)); 
		    LogDebug("rpcefficiency")<<"CSC  \t \t \t \t \t Candidate"<<rollId<<" "<<"(from CSC Segment) STRIP---> "<<stripPredicted<< std::endl;
		    //--------- HISTOGRAM STRIP PREDICTED FROM CSC  -------------------
		    std::map<std::string, MonitorElement*> meMap=meCollection[rpcId.rawId()];
		    meIdCSC.str("");
		    meIdCSC<<"ExpectedOccupancyFromCSC_"<<rollId.rawId();
		    meMap[meIdCSC.str()]->Fill(stripPredicted);
		    //--------------------------------------------------------------------
		    
		    
		    //-------RecHitPart Just For Residual--------
		    int cluSize = 0;
		    int countRecHits = 0;
		    float minres = 3000.;
		    
		    LogDebug("rpcefficiency")<<"CSC  \t \t \t \t \t Getting RecHits in Roll Asociated";
		    typedef std::pair<RPCRecHitCollection::const_iterator, RPCRecHitCollection::const_iterator> rangeRecHits;
		    rangeRecHits recHitCollection =  rpcHits->get(rollasociated->id());
		    RPCRecHitCollection::const_iterator recHit;
		    
		    for (recHit = recHitCollection.first; recHit != recHitCollection.second ; recHit++) {
		      
		      countRecHits++;
		      LocalPoint recHitPos=recHit->localPosition();
		      float res=PointExtrapolatedRPCFrame.x()- recHitPos.x();
		      LogDebug("rpcefficiency")<<"CSC  \t \t \t \t \t \t Found Rec Hit at "<<res<<"cm of the prediction.";
		      if(fabs(res)<fabs(minres)){
			minres=res;
			cluSize = recHit->clusterSize();
			LogDebug("rpcefficiency")<<"CSC  \t \t \t \t \t \t \t New Min Res "<<res<<"cm.";
		      }
		    }
		    
		    if(countRecHits==0){
		      LogDebug("rpcefficiency") <<"CSC \t \t \t \t \t THIS ROLL DOESN'T HAVE ANY RECHIT";
		    }else{  
		      assert(minres!=3000); 
		      
		      if(fabs(minres)<=(rangestrips+cluSize*0.5)*stripw){
			LogDebug("rpcefficiency")<<"CSC  \t \t \t \t \t \t True!";
			
			if(rollId.ring()==2&&rollId.roll()==1){
			  if(cluSize==1*dupli) { hGlobalResClu1R2A->Fill(minres);} 
			  else if(cluSize==2*dupli) { hGlobalResClu2R2A->Fill(minres); }
			  else if(cluSize==3*dupli)  {hGlobalResClu3R2A->Fill(minres);}
			}
			else if(rollId.ring()==2&&rollId.roll()==2){
			  if(cluSize==1*dupli) { hGlobalResClu1R2B->Fill(minres);} 
			  else if(cluSize==2*dupli) { hGlobalResClu2R2B->Fill(minres); }
			  else if(cluSize==3*dupli)  {hGlobalResClu3R2B->Fill(minres);}
			}
			else if(rollId.ring()==2&&rollId.roll()==3){
			  if(cluSize==1*dupli) { hGlobalResClu1R2C->Fill(minres);} 
			  else if(cluSize==2*dupli)  {hGlobalResClu2R2C->Fill(minres); }
			  else if(cluSize==3*dupli) { hGlobalResClu3R2C->Fill(minres);}
			}
			else if(rollId.ring()==3&&rollId.roll()==1){
			  if(cluSize==1*dupli)  {hGlobalResClu1R3A->Fill(minres); }
			  else if(cluSize==2*dupli) { hGlobalResClu2R3A->Fill(minres); }
			  else if(cluSize==3*dupli)  {hGlobalResClu3R3A->Fill(minres);}
			}
			else if(rollId.ring()==3&&rollId.roll()==2){
			  if(cluSize==1*dupli) { hGlobalResClu1R3B->Fill(minres); }
			  else if(cluSize==2*dupli) { hGlobalResClu2R3B->Fill(minres); }
			  else if(cluSize==3*dupli) { hGlobalResClu3R3B->Fill(minres);}
			}
			else if(rollId.ring()==3&&rollId.roll()==3){
			  if(cluSize==1*dupli) {hGlobalResClu1R3C->Fill(minres); }
			  else if(cluSize==2*dupli) {hGlobalResClu2R3C->Fill(minres); }
			  else if(cluSize==3*dupli) {hGlobalResClu3R3C->Fill(minres);}
			}
			meIdRPC.str("");
			meIdRPC<<"RPCDataOccupancyFromCSC_"<<rollId.rawId();
			meMap[meIdRPC.str()]->Fill(stripPredicted);
		      }
		    }
		    
		  }else{
		    LogDebug("rpcefficiency")<<"CSC \t \t \t \t No the prediction is outside of this roll";
		  }//Condition for the right match
		}else{//if extrapolation distance D is not too long
		  LogDebug("rpcefficiency")<<"CSC \t \t \t No, Exrtrapolation too long!, canceled";
		}//D so big
	      }//loop over the rolls asociated 
	    }//Condition over the LS1 geometry!!!!
	  }//Is the segment 4D?
	}else{
	  LogDebug("rpcefficiency")<<"CSC \t \t More than one segment in this chamber, or we are in Station Ring 1 or in Station 4";
	}
      }
    }else{
      LogDebug("rpcefficiency")<<"CSC This Event doesn't have any CSCSegment";
    }
      */
    }//if(allCSCSegments->size()>0)
  }//if(allCSCSegments.isValid())





/////////////////////////////////////////////////////////////////////////////////////////////
    /*necalSC_ = 0;
	calSCindex_.clear();
	ecalSCeta_.clear();
	ecalSCphi_.clear();
	ecalSCEn_.clear();
	ecalSCRawEn_.clear();
	ecalSCetaWidth_.clear();
	ecalSCphiWidth_.clear();
	ecalSC_LICTD_.clear();
	ecalSC_nL1Spike_.clear();
	ecalSC_nDiweird_.clear();
	ecalSC_nWeird_.clear();
	ecalSC_nSaturated_.clear();
	ecalSC_nOutOfTime_.clear();
	ecalSC_nXtals_.clear();
	ecalSC_maxEnXtalTime_.clear();
	ecalSC_maxEnXtalSwissCross_.clear();
	ecalSC_maxEnXtalBits_.clear();


	edm::Handle<std::vector<reco::SuperCluster>> ecalSChandle;
	e.getByToken(ecalSCcollection_, ecalSChandle);

	noZS::EcalClusterLazyTools lazyToolnoZS(e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);

	if(ecalSChandle.isValid()){

		UShort_t _scIndex_ = 0;
		for(std::vector<reco::SuperCluster>::const_iterator iSC = ecalSChandle->begin(); iSC != ecalSChandle->end(); iSC++){
			_scIndex_++;

			necalSC_++;
			Float_t maxEnXtalTime = -999.;
			UChar_t tmpnL1Spike	= 0;
			UChar_t tmpnDiweird = 0;
			UChar_t tmpnWeird = 0;
			UChar_t tmpnSaturated = 0;
			UChar_t tmpnOutOfTime = 0;
			UChar_t tmpnXtals = 0;
			Float_t tmpmaxEnXtalSwissCross = -999;
			UChar_t tmpmaxEnXtalBits = 0;
			Float_t tmpscEta = iSC->eta();
			Float_t tmpscPhi = iSC->phi();
			Float_t tmpecalSC_LICTD =  getLICTD(&(*iSC), lazyToolnoZS, maxEnXtalTime, tmpnL1Spike, tmpnDiweird, tmpnWeird, tmpnSaturated, tmpnOutOfTime, tmpnXtals, tmpmaxEnXtalBits, tmpmaxEnXtalSwissCross);

			ecalSCindex_.push_back(_scIndex_-1);
			ecalSCeta_.push_back(tmpscEta);
			ecalSCphi_.push_back(tmpscPhi);
			ecalSCEn_.push_back(iSC->energy());
			ecalSCRawEn_.push_back(iSC->rawEnergy());
			ecalSCetaWidth_.push_back(iSC->etaWidth());
			ecalSCphiWidth_.push_back(iSC->phiWidth());
			ecalSC_LICTD_.push_back(tmpecalSC_LICTD);
			ecalSC_nL1Spike_.push_back(tmpnL1Spike);
			ecalSC_nDiweird_.push_back(tmpnDiweird);
			ecalSC_nWeird_.push_back(tmpnWeird);
			ecalSC_nSaturated_.push_back(tmpnSaturated);
			ecalSC_nOutOfTime_.push_back(tmpnOutOfTime);
			ecalSC_nXtals_.push_back(tmpnXtals);
			ecalSC_maxEnXtalTime_.push_back(maxEnXtalTime);
			ecalSC_maxEnXtalSwissCross_.push_back(tmpmaxEnXtalSwissCross);
			ecalSC_maxEnXtalBits_.push_back(tmpmaxEnXtalBits);
		}
	}
    */
};




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////         OOT ECAL SC     ////////////////////////////////////////////////////////
/*
UShort_t necalootSC_;
std::vector<UShort_t> ecalootSCindex_;
std::vector<Float_t> ecalootSCeta_;
std::vector<Float_t> ecalootSCphi_;
std::vector<Float_t> ecalootSCEn_;
std::vector<Float_t> ecalootSCRawEn_;
std::vector<Float_t> ecalootSCetaWidth_;
std::vector<Float_t> ecalootSCphiWidth_;
std::vector<Float_t> ecalootSC_LICTD_;
std::vector<UChar_t> ecalootSC_nL1Spike_;
std::vector<UChar_t> ecalootSC_nDiweird_;
std::vector<UChar_t> ecalootSC_nWeird_;
std::vector<UChar_t> ecalootSC_nSaturated_;
std::vector<UChar_t> ecalootSC_nOutOfTime_;
std::vector<UChar_t> ecalootSC_nXtals_;
std::vector<Float_t> ecalootSC_maxEnXtalTime_;
std::vector<Float_t> ecalootSC_maxEnXtalSwissCross_;
std::vector<UChar_t> ecalootSC_maxEnXtalBits_;



void ggNtuplizer::branchesECALOOTSC(TTree* tree) {
	tree->Branch("necalootSC",                    		&necalootSC_);
	tree->Branch("ecalootSC_index",                    	&ecalootSCindex_);
	tree->Branch("ecalootSC_eta",                    	&ecalootSCeta_);
	tree->Branch("ecalootSC_phi",                    	&ecalootSCphi_);
	tree->Branch("ecalootSC_En",                    	&ecalootSCEn_);
	tree->Branch("ecalootSC_RawEn",                    	&ecalootSCRawEn_);
	tree->Branch("ecalootSC_etaWidth",                  &ecalootSCetaWidth_);
	tree->Branch("ecalootSC_phiWidth",                  &ecalootSCphiWidth_);
	tree->Branch("ecalootSC_LICTD",                    &ecalootSC_LICTD_);
	tree->Branch("ecalootSC_nL1Spike",                 &ecalootSC_nL1Spike_);
	tree->Branch("ecalootSC_nDiweird",                 &ecalootSC_nDiweird_);
	tree->Branch("ecalootSC_nWeird",					&ecalootSC_nWeird_);
	tree->Branch("ecalootSC_nSaturated",				&ecalootSC_nSaturated_);
	tree->Branch("ecalootSC_nOutOfTime",				&ecalootSC_nOutOfTime_);
	tree->Branch("ecalootSC_nXtals",					&ecalootSC_nXtals_);
	tree->Branch("ecalootSC_maxEnXtalTime",			&ecalootSC_maxEnXtalTime_);
	tree->Branch("ecalootSC_maxEnXtalSwissCross",		&ecalootSC_maxEnXtalSwissCross_);
	tree->Branch("ecalootSC_maxEnXtalBits",			&ecalootSC_maxEnXtalBits_);
};


void ggNtuplizer::fillECALOOTSC(const edm::Event& e, const edm::EventSetup& es){
	necalootSC_ = 0;
	ecalootSCindex_.clear();
	ecalootSCeta_.clear();
	ecalootSCphi_.clear();
	ecalootSCEn_.clear();
	ecalootSCRawEn_.clear();
	ecalootSCetaWidth_.clear();
	ecalootSCphiWidth_.clear();
	ecalootSC_LICTD_.clear();
	ecalootSC_nL1Spike_.clear();
	ecalootSC_nDiweird_.clear();
	ecalootSC_nWeird_.clear();
	ecalootSC_nSaturated_.clear();
	ecalootSC_nOutOfTime_.clear();
	ecalootSC_nXtals_.clear();
	ecalootSC_maxEnXtalTime_.clear();
	ecalootSC_maxEnXtalSwissCross_.clear();
	ecalootSC_maxEnXtalBits_.clear();


	edm::Handle<std::vector<reco::SuperCluster>> ecalootSChandle;
	e.getByToken(ecalSC_OOT_collection_, ecalootSChandle);

	noZS::EcalClusterLazyTools lazyToolnoZS(e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);

	if(ecalootSChandle.isValid()){

		UShort_t _scIndex_ = 0;
		for(std::vector<reco::SuperCluster>::const_iterator iSC = ecalootSChandle->begin(); iSC != ecalootSChandle->end(); iSC++){
			_scIndex_++;

			necalootSC_++;
			Float_t maxEnXtalTime = -999.;
			UChar_t tmpnL1Spike	= 0;
			UChar_t tmpnDiweird = 0;
			UChar_t tmpnWeird = 0;
			UChar_t tmpnSaturated = 0;
			UChar_t tmpnOutOfTime = 0;
			UChar_t tmpnXtals = 0;
			Float_t tmpmaxEnXtalSwissCross = -999;
			UChar_t tmpmaxEnXtalBits = 0;
			Float_t tmpscEta = iSC->eta();
			Float_t tmpscPhi = iSC->phi();
			Float_t tmpecalootSC_LICTD =  getLICTD(&(*iSC), lazyToolnoZS, maxEnXtalTime, tmpnL1Spike, tmpnDiweird, tmpnWeird, tmpnSaturated, tmpnOutOfTime, tmpnXtals, tmpmaxEnXtalBits, tmpmaxEnXtalSwissCross);

			ecalootSCindex_.push_back(_scIndex_-1);
			ecalootSCeta_.push_back(tmpscEta);
			ecalootSCphi_.push_back(tmpscPhi);
			ecalootSCEn_.push_back(iSC->energy());
			ecalootSCRawEn_.push_back(iSC->rawEnergy());
			ecalootSCetaWidth_.push_back(iSC->etaWidth());
			ecalootSCphiWidth_.push_back(iSC->phiWidth());
			ecalootSC_LICTD_.push_back(tmpecalootSC_LICTD);
			ecalootSC_nL1Spike_.push_back(tmpnL1Spike);
			ecalootSC_nDiweird_.push_back(tmpnDiweird);
			ecalootSC_nWeird_.push_back(tmpnWeird);
			ecalootSC_nSaturated_.push_back(tmpnSaturated);
			ecalootSC_nOutOfTime_.push_back(tmpnOutOfTime);
			ecalootSC_nXtals_.push_back(tmpnXtals);
			ecalootSC_maxEnXtalTime_.push_back(maxEnXtalTime);
			ecalootSC_maxEnXtalSwissCross_.push_back(tmpmaxEnXtalSwissCross);
			ecalootSC_maxEnXtalBits_.push_back(tmpmaxEnXtalBits);
		}
	}
};
*/
