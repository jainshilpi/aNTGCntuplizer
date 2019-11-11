#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"

using namespace std;

UShort_t metFilters_;
Float_t genMET_;
Float_t genMETPhi_;
Float_t pfMET_;
Float_t pfMETPhi_;
Float_t pfMET_T1JERUp_;
Float_t pfMET_T1JERDo_;
Float_t pfMET_T1JESUp_;
Float_t pfMET_T1JESDo_;
Float_t pfMET_T1MESUp_;
Float_t pfMET_T1MESDo_;
Float_t pfMET_T1EESUp_;
Float_t pfMET_T1EESDo_;
Float_t pfMET_T1PESUp_;
Float_t pfMET_T1PESDo_;
Float_t pfMET_T1TESUp_;
Float_t pfMET_T1TESDo_;
Float_t pfMET_T1UESUp_;
Float_t pfMET_T1UESDo_;
Float_t pfMET_T1TxyPhi_;
Float_t pfMET_T1TxyPt_;
Float_t pfMETPhi_T1JESUp_;
Float_t pfMETPhi_T1JESDo_;
Float_t pfMETPhi_T1UESUp_;
Float_t pfMETPhi_T1UESDo_;
// Float_t pfMET_caloMetSig_;
Float_t pfMET_metSig_;
Float_t pfMET_EtSig_;

const std::vector<std::string> filterNamesToCheck = {
  "Flag_goodVertices",
  "Flag_globalSuperTightHalo2016Filter",
  "Flag_HBHENoiseFilter",
  "Flag_HBHENoiseIsoFilter",
  "Flag_EcalDeadCellTriggerPrimitiveFilter",
  "Flag_BadPFMuonFilter",
  "Flag_BadChargedCandidateFilter"
  "Flag_eeBadScFilter",
  "Flag_ecalBadCalibFilter",
};

void ggNtuplizer::branchesMET(TTree* tree) {
  if (doGenParticles_) {
    tree->Branch("genMET",      &genMET_);
    tree->Branch("genMETPhi",   &genMETPhi_);
  }
  tree->Branch("metFilters",       &metFilters_);
  tree->Branch("pfMET",            &pfMET_);
  tree->Branch("pfMETPhi",         &pfMETPhi_);
  tree->Branch("pfMET_T1JERUp",    &pfMET_T1JERUp_);
  tree->Branch("pfMET_T1JERDo",    &pfMET_T1JERDo_);
  tree->Branch("pfMET_T1JESUp",    &pfMET_T1JESUp_);
  tree->Branch("pfMET_T1JESDo",    &pfMET_T1JESDo_);
  // tree->Branch("pfMET_T1MESUp",    &pfMET_T1MESUp_);
  // tree->Branch("pfMET_T1MESDo",    &pfMET_T1MESDo_);
  // tree->Branch("pfMET_T1EESUp",    &pfMET_T1EESUp_);
  // tree->Branch("pfMET_T1EESDo",    &pfMET_T1EESDo_);
  // tree->Branch("pfMET_T1PESUp",    &pfMET_T1PESUp_);
  // tree->Branch("pfMET_T1PESDo",    &pfMET_T1PESDo_);
  // tree->Branch("pfMET_T1TESUp",    &pfMET_T1TESUp_);
  // tree->Branch("pfMET_T1TESDo",    &pfMET_T1TESDo_);
  tree->Branch("pfMET_T1UESUp",    &pfMET_T1UESUp_);
  tree->Branch("pfMET_T1UESDo",    &pfMET_T1UESDo_);
  tree->Branch("pfMETPhi_T1JESUp", &pfMETPhi_T1JESUp_);
  tree->Branch("pfMETPhi_T1JESDo", &pfMETPhi_T1JESDo_);
  tree->Branch("pfMETPhi_T1UESUp", &pfMETPhi_T1UESUp_);
  tree->Branch("pfMETPhi_T1UESDo", &pfMETPhi_T1UESDo_);
  // tree->Branch("pfMET_caloMetSig", &pfMET_caloMetSig_);
  tree->Branch("pfMET_metSig", &pfMET_metSig_);
  tree->Branch("pfMET_EtSig", &pfMET_EtSig_);
}

void ggNtuplizer::fillMET(const edm::Event& e, const edm::EventSetup& es) {
  metFilters_ = 0;
  if (addFilterInfoMINIAOD_) {
    edm::Handle<edm::TriggerResults> patFilterResultsHandle;
    e.getByToken(patTrgResultsLabel_, patFilterResultsHandle);
    edm::TriggerResults const& patFilterResults = *patFilterResultsHandle;

    auto&& filterNames = e.triggerNames(patFilterResults);

    ////=== the following lines allow us to find the filters stored in the event ! ===
    // edm::TriggerNames const& triggerNames = e.triggerNames(patFilterResults);
    // for ( edm::TriggerNames::Strings::const_iterator triggerName = triggerNames.triggerNames().begin(); triggerName != triggerNames.triggerNames().end(); ++triggerName){
    //   int triggerId = triggerNames.triggerIndex(*triggerName);
    //   if ( triggerId >= 0 && triggerId < (int)triggerNames.size() ){
    //     std::string triggerDecision = ( patFilterResultsHandle->accept(triggerId) ) ? "passed" : "failed";
    //     std::cout << " triggerName = " << (*triggerName) << " " << triggerDecision << std::endl;
    //   }
    // }
    for (unsigned iF = 0; iF < filterNamesToCheck.size(); iF++){
      unsigned index = filterNames.triggerIndex(filterNamesToCheck[iF]);
      if( index == filterNames.size() ) LogDebug("METFilters") << filterNamesToCheck[iF] << " is missing, exiting";
      else if (!patFilterResults.accept(index)) metFilters_ += std::pow(2, iF);
    }

    //////https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#Accessing_the_filter_decision_in
    edm::Handle<bool> passecalBadCalibFilterUpdate ;
    e.getByToken(ecalBadCalibFilterUpdateToken_, passecalBadCalibFilterUpdate);
    Bool_t _passecalBadCalibFilterUpdate =  (*passecalBadCalibFilterUpdate );
    if(!_passecalBadCalibFilterUpdate) metFilters_ += std::pow(2, filterNamesToCheck.size());
  }

  edm::Handle<edm::View<pat::MET> > pfMETHandle;
  e.getByToken(pfMETlabel_, pfMETHandle);

  genMET_    = -99;
  genMETPhi_ = -99;
  pfMET_     = -99;
  pfMETPhi_  = -99;

  if (pfMETHandle.isValid()) {
    const pat::MET *pfMET = 0;
    pfMET     = &(pfMETHandle->front());
    pfMET_    = pfMET->et();
    pfMETPhi_ = pfMET->phi();

    // Type1MET uncertainties =======================================
    pfMET_T1JERUp_ = pfMET->shiftedPt(pat::MET::JetResUp);
    pfMET_T1JERDo_ = pfMET->shiftedPt(pat::MET::JetResDown);
    pfMET_T1JESUp_ = pfMET->shiftedPt(pat::MET::JetEnUp);
    pfMET_T1JESDo_ = pfMET->shiftedPt(pat::MET::JetEnDown);

    // pfMET_T1MESUp_ = pfMET->shiftedPt(pat::MET::MuonEnUp);
    // pfMET_T1MESDo_ = pfMET->shiftedPt(pat::MET::MuonEnDown);
    // pfMET_T1EESUp_ = pfMET->shiftedPt(pat::MET::ElectronEnUp);
    // pfMET_T1EESDo_ = pfMET->shiftedPt(pat::MET::ElectronEnDown);
    // pfMET_T1PESUp_ = pfMET->shiftedPt(pat::MET::PhotonEnUp);
    // pfMET_T1PESDo_ = pfMET->shiftedPt(pat::MET::PhotonEnDown);
    // pfMET_T1TESUp_ = pfMET->shiftedPt(pat::MET::TauEnUp);
    // pfMET_T1TESDo_ = pfMET->shiftedPt(pat::MET::TauEnDown);

    pfMET_T1UESUp_ = pfMET->shiftedPt(pat::MET::UnclusteredEnUp);
    pfMET_T1UESDo_ = pfMET->shiftedPt(pat::MET::UnclusteredEnDown);

    pfMETPhi_T1JESUp_ = pfMET->shiftedPhi(pat::MET::JetEnUp);
    pfMETPhi_T1JESDo_ = pfMET->shiftedPhi(pat::MET::JetEnDown);
    pfMETPhi_T1UESUp_ = pfMET->shiftedPhi(pat::MET::UnclusteredEnUp);
    pfMETPhi_T1UESDo_ = pfMET->shiftedPhi(pat::MET::UnclusteredEnDown);

    // pfMET_caloMetSig_ = pfMET->caloMetSignificance();
    pfMET_metSig_ = pfMET->metSignificance();
    pfMET_EtSig_ = pfMET->mEtSig();

    if (!e.isRealData()) {
      genMET_    = pfMET->genMET()->et();
      genMETPhi_ = pfMET->genMET()->phi();
    }
  }
}
