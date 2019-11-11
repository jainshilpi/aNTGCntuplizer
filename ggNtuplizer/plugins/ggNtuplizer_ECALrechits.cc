#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

UShort_t necalSC_;
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



void ggNtuplizer::branchesECALSC(TTree* tree) {
  tree->Branch("necalSC",                       &necalSC_);
  // tree->Branch("ecalSCindex",                     &ecalSCindex_);
  tree->Branch("ecalSCeta",                     &ecalSCeta_);
  tree->Branch("ecalSCphi",                     &ecalSCphi_);
  tree->Branch("ecalSCEn",                      &ecalSCEn_);
  tree->Branch("ecalSCRawEn",                     &ecalSCRawEn_);
  tree->Branch("ecalSCetaWidth",                  &ecalSCetaWidth_);
  tree->Branch("ecalSCphiWidth",                  &ecalSCphiWidth_);
  tree->Branch("ecalSC_LICTD",                    &ecalSC_LICTD_);
  tree->Branch("ecalSC_nL1Spike",                 &ecalSC_nL1Spike_);
  tree->Branch("ecalSC_nDiweird",                 &ecalSC_nDiweird_);
  tree->Branch("ecalSC_nWeird",         &ecalSC_nWeird_);
  tree->Branch("ecalSC_nSaturated",       &ecalSC_nSaturated_);
  tree->Branch("ecalSC_nOutOfTime",       &ecalSC_nOutOfTime_);
  tree->Branch("ecalSC_nXtals",         &ecalSC_nXtals_);
  tree->Branch("ecalSC_maxEnXtalTime",      &ecalSC_maxEnXtalTime_);
  tree->Branch("ecalSC_maxEnXtalSwissCross",    &ecalSC_maxEnXtalSwissCross_);
  tree->Branch("ecalSC_maxEnXtalBits",      &ecalSC_maxEnXtalBits_);
};



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

  const EcalRecHit *maxEnergyXtalRecHit_  = nullptr;
  Float_t maxEnergyXtal_          = std::numeric_limits<Float_t>::min();
  Float_t maxTime_            = std::numeric_limits<Float_t>::min();
  Float_t minTime_            = std::numeric_limits<Float_t>::max();
  Bool_t  _timeValid            = 0;

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
  _maxEnXtalTime    = maxEnergyXtalRecHit_->time();
  _maxEnXtalBits    = 0;
  if(maxEnergyXtalRecHit_->checkFlag(EcalRecHit::kDiWeird)) setbit(_maxEnXtalBits, 0);
  if(maxEnergyXtalRecHit_->checkFlag(EcalRecHit::kWeird)) setbit(_maxEnXtalBits, 1);
  if(maxEnergyXtalRecHit_->checkFlag(EcalRecHit::kSaturated)) setbit(_maxEnXtalBits, 2);
  if(maxEnergyXtalRecHit_->checkFlag(EcalRecHit::kOutOfTime)) setbit(_maxEnXtalBits, 3);

  _maxEnXtalSwissCross = swissCross(maxEnergyXtalRecHit_->detid(), ltNoZS);
  if(maxEnergyXtalRecHit_->detid().rawId() != seedDetID.rawId()) setbit(_maxEnXtalBits, 4);

  Float_t LICTD_ = (_timeValid) ? std::max(std::abs(maxTime_ - _maxEnXtalTime), std::abs(_maxEnXtalTime - minTime_)) : -9999999.;

  return LICTD_;
};



void ggNtuplizer::fillECALSC(const edm::Event& e, const edm::EventSetup& es){
  necalSC_ = 0;
  ecalSCindex_.clear();
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
      UChar_t tmpnL1Spike = 0;
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
};




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////         OOT ECAL SC     ////////////////////////////////////////////////////////
UShort_t necalootSC_;
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
  tree->Branch("necalootSC",                        &necalootSC_);
  // tree->Branch("ecalootSC_index",                     &ecalootSCindex_);
  tree->Branch("ecalootSC_eta",                     &ecalootSCeta_);
  tree->Branch("ecalootSC_phi",                     &ecalootSCphi_);
  tree->Branch("ecalootSC_En",                      &ecalootSCEn_);
  tree->Branch("ecalootSC_RawEn",                     &ecalootSCRawEn_);
  tree->Branch("ecalootSC_etaWidth",                  &ecalootSCetaWidth_);
  tree->Branch("ecalootSC_phiWidth",                  &ecalootSCphiWidth_);
  tree->Branch("ecalootSC_LICTD",                    &ecalootSC_LICTD_);
  tree->Branch("ecalootSC_nL1Spike",                 &ecalootSC_nL1Spike_);
  tree->Branch("ecalootSC_nDiweird",                 &ecalootSC_nDiweird_);
  tree->Branch("ecalootSC_nWeird",          &ecalootSC_nWeird_);
  tree->Branch("ecalootSC_nSaturated",        &ecalootSC_nSaturated_);
  tree->Branch("ecalootSC_nOutOfTime",        &ecalootSC_nOutOfTime_);
  tree->Branch("ecalootSC_nXtals",          &ecalootSC_nXtals_);
  tree->Branch("ecalootSC_maxEnXtalTime",     &ecalootSC_maxEnXtalTime_);
  tree->Branch("ecalootSC_maxEnXtalSwissCross",   &ecalootSC_maxEnXtalSwissCross_);
  tree->Branch("ecalootSC_maxEnXtalBits",     &ecalootSC_maxEnXtalBits_);
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
      UChar_t tmpnL1Spike = 0;
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



Short_t ggNtuplizer::findSecondaryIndex(const Short_t & searchIndex, const std::vector<UShort_t> & container){
  Short_t secondaryIndex = -999;
  for(UInt_t i = 0; i < container.size(); i++){
    if( container[i] ==  searchIndex) {
      secondaryIndex = i;
      break;
    }
  }
  return secondaryIndex;
};
