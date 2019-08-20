#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "Math/VectorUtil.h"

using namespace std;

UShort_t          nPho_;
vector<Float_t>  phoE_;
vector<Float_t>  phoSigmaE_;
vector<Float_t>  phoEt_;
vector<Float_t>  phoEta_;
vector<Float_t>  phoPhi_;
vector<Float_t>  phoCalibE_;
vector<Float_t>  phoSigmaCalibE_;
vector<Float_t>  phoCalibEt_;
std::vector<Short_t> phoSCindex_;
vector<Float_t>  phoESEnP1_;
vector<Float_t>  phoESEnP2_;
vector<UChar_t> phoFiducialRegion_;
vector<UChar_t>   phoQualityBits_;
vector<Float_t>  phoR9_;
vector<Float_t>  phoHoverE_;
vector<Float_t>  phoESEffSigmaRR_;
vector<Float_t>  phoSigmaIEtaIEtaFull5x5_;
vector<Float_t>  phoSigmaIEtaIPhiFull5x5_;
vector<Float_t>  phoSigmaIPhiIPhiFull5x5_;
vector<Float_t>  phoE2x2Full5x5_;
vector<Float_t>  phoE5x5Full5x5_;
vector<Float_t>  phoR9Full5x5_;
vector<Float_t>  phoPFChIso_;
vector<Float_t>  phoPFPhoIso_;
vector<Float_t>  phoPFNeuIso_;
vector<Float_t>  phoPFChWorstIso_;
vector<Float_t>  phoSeedBCE_;
vector<Float_t>  phoSeedBCEta_;
vector<Float_t>  phoIDMVA_;
vector<ULong64_t> phoFiredSingleTrgs_;
vector<ULong64_t> phoFiredDoubleTrgs_;
vector<ULong64_t> phoFiredTripleTrgs_;
vector<ULong64_t> phoFiredL1Trgs_;
vector<Float_t>  phoSeedTime_;
vector<Float_t>  phoSeedEnergy_;
vector<Float_t>  phoMIPChi2_;
vector<Float_t>  phoMIPTotEnergy_;
vector<Float_t>  phoMIPSlope_;
vector<Float_t>  phoMIPIntercept_;
vector<Short_t>  phoMIPNhitCone_;
vector<UChar_t> phoIDbit_;
vector<Float_t>    phoScale_stat_up_;
vector<Float_t>    phoScale_stat_dn_;
vector<Float_t>    phoScale_syst_up_;
vector<Float_t>    phoScale_syst_dn_;
vector<Float_t>    phoScale_gain_up_;
vector<Float_t>    phoScale_gain_dn_;
vector<Float_t>    phoResol_rho_up_;
vector<Float_t>    phoResol_rho_dn_;
vector<Float_t>    phoResol_phi_up_;
vector<Float_t>    phoResol_phi_dn_;
vector<Short_t> pho_gen_index_;
vector<float>  phoSCEta_;
vector<float>  phoSCPhi_;
vector<float>  phoSCEtaWidth_;
vector<float>  phoSCPhiWidth_;
vector<float>  phoSCBrem_;
vector<float>  phoSCE_;
vector<float>  phoSCRawE_;

//Necessary for the Photon Footprint removal
template <class T, class U>
bool isInFootprint(const T& thefootprint, const U& theCandidate) {
  for ( auto itr = thefootprint.begin(); itr != thefootprint.end(); ++itr ) {

    if( itr.key() == theCandidate.key() ) return true;

  }
  return false;
};


void ggNtuplizer::branchesPhotons(TTree* tree) {
  tree->Branch("nPho",                    &nPho_);
  tree->Branch("phoE",                    &phoE_);
  tree->Branch("phoSigmaE",               &phoSigmaE_);
  tree->Branch("phoEt",                   &phoEt_);
  tree->Branch("phoEta",                  &phoEta_);
  tree->Branch("phoPhi",                  &phoPhi_);
  tree->Branch("phoCalibE",               &phoCalibE_);
  tree->Branch("phoSigmaCalibE",          &phoSigmaCalibE_);
  tree->Branch("phoCalibEt",              &phoCalibEt_);
  tree->Branch("phoSCindex",              &phoSCindex_);
  tree->Branch("phoESEnP1",               &phoESEnP1_);
  tree->Branch("phoESEnP2",               &phoESEnP2_);
  tree->Branch("phoFiducialRegion",       &phoFiducialRegion_);
  tree->Branch("phoQualityBits",              &phoQualityBits_);
  tree->Branch("phoR9",                   &phoR9_);
  tree->Branch("phoHoverE",               &phoHoverE_);
  tree->Branch("phoESEffSigmaRR",         &phoESEffSigmaRR_);
  tree->Branch("phoSigmaIEtaIEtaFull5x5", &phoSigmaIEtaIEtaFull5x5_);
  tree->Branch("phoSigmaIEtaIPhiFull5x5", &phoSigmaIEtaIPhiFull5x5_);
  tree->Branch("phoSigmaIPhiIPhiFull5x5", &phoSigmaIPhiIPhiFull5x5_);
  tree->Branch("phoE2x2Full5x5",          &phoE2x2Full5x5_);
  tree->Branch("phoE5x5Full5x5",          &phoE5x5Full5x5_);
  tree->Branch("phoR9Full5x5",            &phoR9Full5x5_);
  tree->Branch("phoSeedBCE",              &phoSeedBCE_);
  tree->Branch("phoSeedBCEta",            &phoSeedBCEta_);
  tree->Branch("phoPFChIso",              &phoPFChIso_);
  tree->Branch("phoPFPhoIso",             &phoPFPhoIso_);
  tree->Branch("phoPFNeuIso",             &phoPFNeuIso_);
  tree->Branch("phoPFChWorstIso",         &phoPFChWorstIso_);
  tree->Branch("phoIDMVA",                &phoIDMVA_);
  tree->Branch("phoFiredSingleTrgs",      &phoFiredSingleTrgs_);
  tree->Branch("phoFiredDoubleTrgs",      &phoFiredDoubleTrgs_);
  tree->Branch("phoFiredTripleTrgs",      &phoFiredTripleTrgs_);
  tree->Branch("phoFiredL1Trgs",          &phoFiredL1Trgs_);
  tree->Branch("phoSeedTime",             &phoSeedTime_);
  tree->Branch("phoSeedEnergy",           &phoSeedEnergy_);
  // tree->Branch("phoSeedTimeFull5x5",              &phoSeedTimeFull5x5_);
  tree->Branch("phoMIPChi2",                      &phoMIPChi2_);
  tree->Branch("phoMIPTotEnergy",                 &phoMIPTotEnergy_);
  tree->Branch("phoMIPSlope",                     &phoMIPSlope_);
  tree->Branch("phoMIPIntercept",                 &phoMIPIntercept_);
  tree->Branch("phoMIPNhitCone",                  &phoMIPNhitCone_);
  // tree->Branch("phoxtalBits",      &phoxtalBits_);
  tree->Branch("phoIDbit",         &phoIDbit_);
  tree->Branch("phoScale_stat_up", &phoScale_stat_up_);
  tree->Branch("phoScale_stat_dn", &phoScale_stat_dn_);
  tree->Branch("phoScale_syst_up", &phoScale_syst_up_);
  tree->Branch("phoScale_syst_dn", &phoScale_syst_dn_);
  tree->Branch("phoScale_gain_up", &phoScale_gain_up_);
  tree->Branch("phoScale_gain_dn", &phoScale_gain_dn_);
  tree->Branch("phoResol_rho_up",  &phoResol_rho_up_);
  tree->Branch("phoResol_rho_dn",  &phoResol_rho_dn_);
  tree->Branch("phoResol_phi_up",  &phoResol_phi_up_);
  tree->Branch("phoResol_phi_dn",  &phoResol_phi_dn_);

  tree->Branch("phoSCEta",                &phoSCEta_);
  tree->Branch("phoSCPhi",                &phoSCPhi_);
  tree->Branch("phoSCE",                  &phoSCE_);
  tree->Branch("phoSCEtaWidth",           &phoSCEtaWidth_);
  tree->Branch("phoSCPhiWidth",           &phoSCPhiWidth_);
  tree->Branch("phoSCBrem",               &phoSCBrem_);
  tree->Branch("phoSCRawE",               &phoSCRawE_);


  if(doGenParticles_){
    tree->Branch("pho_gen_index",  &pho_gen_index_);
  }
}

void ggNtuplizer::fillPhotons(const edm::Event& e, const edm::EventSetup& es) {

  // cleanup from previous execution
  phoE_                   .clear();
  phoSigmaE_              .clear();
  phoEt_                  .clear();
  phoEta_                 .clear();
  phoPhi_                 .clear();
  phoCalibE_              .clear();
  phoSigmaCalibE_         .clear();
  phoCalibEt_             .clear();
  phoSCindex_             .clear();
  phoESEnP1_              .clear();
  phoESEnP2_              .clear();
  phoFiducialRegion_      .clear();
  phoQualityBits_.clear();
  phoR9_                  .clear();
  phoHoverE_              .clear();
  phoESEffSigmaRR_        .clear();
  phoSigmaIEtaIEtaFull5x5_.clear();
  phoSigmaIEtaIPhiFull5x5_.clear();
  phoSigmaIPhiIPhiFull5x5_.clear();
  phoE2x2Full5x5_         .clear();
  phoE5x5Full5x5_         .clear();
  phoR9Full5x5_           .clear();
  phoPFChIso_             .clear();
  phoPFPhoIso_            .clear();
  phoPFNeuIso_            .clear();
  phoPFChWorstIso_        .clear();
  phoSeedBCE_           .clear();
  phoSeedBCEta_         .clear();
  phoIDMVA_               .clear();
  phoFiredSingleTrgs_     .clear();
  phoFiredDoubleTrgs_     .clear();
  phoFiredTripleTrgs_     .clear();
  phoFiredL1Trgs_         .clear();
  phoSeedTime_            .clear();
  phoSeedEnergy_          .clear();
  phoMIPChi2_           .clear();
  phoMIPTotEnergy_      .clear();
  phoMIPSlope_          .clear();
  phoMIPIntercept_      .clear();
  phoMIPNhitCone_       .clear();
  phoIDbit_        .clear();
  phoScale_stat_up_.clear();
  phoScale_stat_dn_.clear();
  phoScale_syst_up_.clear();
  phoScale_syst_dn_.clear();
  phoScale_gain_up_.clear();
  phoScale_gain_dn_.clear();
  phoResol_rho_up_ .clear();
  phoResol_rho_dn_ .clear();
  phoResol_phi_up_ .clear();
  phoResol_phi_dn_ .clear();
  pho_gen_index_.clear();

  phoSCE_               .clear();
  phoSCEta_             .clear();
  phoSCPhi_             .clear();
  phoSCEtaWidth_        .clear();
  phoSCPhiWidth_        .clear();
  phoSCBrem_            .clear();
  phoSCRawE_            .clear();

  nPho_ = 0;

  edm::Handle<edm::View<pat::Photon> > photonHandle;
  e.getByToken(photonCollection_, photonHandle);

  if (!photonHandle.isValid()) {
    edm::LogWarning("ggNtuplizer") << "no pat::Photons in event";
    return;
  }

  edm::Handle<vector<reco::GenParticle>> genParticlesHandle;
  if(doGenParticles_) e.getByToken(genParticlesCollection_, genParticlesHandle);

  edm::Handle<std::vector<reco::SuperCluster>> ecalSChandle;
  e.getByToken(ecalSCcollection_, ecalSChandle);

  ///specific for AOD////
  edm::Handle<edm::ValueMap<bool> >  loose_id_decisions;
  edm::Handle<edm::ValueMap<bool> >  medium_id_decisions;
  edm::Handle<edm::ValueMap<bool> >  tight_id_decisions;
  edm::Handle<edm::ValueMap<float> > mvaValues;
  edm::Handle<edm::ValueMap<float> > phoChargedIsolationMap;
  edm::Handle<edm::ValueMap<float> > phoNeutralHadronIsolationMap;
  edm::Handle<edm::ValueMap<float> > phoPhotonIsolationMap;
  edm::Handle<edm::ValueMap<float> > phoWorstChargedIsolationMap;

  e.getByToken(phoLooseIdMapToken_ ,  loose_id_decisions);
  e.getByToken(phoMediumIdMapToken_,  medium_id_decisions);
  e.getByToken(phoTightIdMapToken_ ,  tight_id_decisions);
  e.getByToken(phoMVAValuesMapToken_, mvaValues);

  e.getByToken(phoChargedIsolationToken_,       phoChargedIsolationMap);
  e.getByToken(phoNeutralHadronIsolationToken_, phoNeutralHadronIsolationMap);
  e.getByToken(phoPhotonIsolationToken_,        phoPhotonIsolationMap);
  e.getByToken(phoWorstChargedIsolationToken_,  phoWorstChargedIsolationMap);
  

  ///END of specific for AOD///

  /*
  edm::Handle<std::vector<reco::SuperCluster>> ecalSChandleEB;
  e.getByToken(ecalSCcollectionEB_, ecalSChandleEB);

  edm::Handle<std::vector<reco::SuperCluster>> ecalSChandleEE;
  e.getByToken(ecalSCcollectionEE_, ecalSChandleEE);
  */

  EcalClusterLazyTools       lazyTool    (e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);
  noZS::EcalClusterLazyTools lazyToolnoZS(e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);

  for (edm::View<pat::Photon>::const_iterator iPho = photonHandle->begin(); iPho != photonHandle->end(); ++iPho){

    // if(iPho->et() < 15.) continue;

    phoE_             .push_back(iPho->energy());
    //phoCalibE_        .push_back(iPho->userFloat("ecalEnergyPostCorr"));
    phoEt_            .push_back(iPho->et());
    //phoCalibEt_       .push_back(iPho->et()*iPho->userFloat("ecalEnergyPostCorr")/iPho->energy());
    //phoSigmaE_        .push_back(iPho->userFloat("ecalEnergyErrPreCorr"));
    //phoSigmaCalibE_   .push_back(iPho->userFloat("ecalEnergyErrPostCorr"));
    phoEta_           .push_back(iPho->eta());
    phoPhi_           .push_back(iPho->phi());
    phoESEnP1_        .push_back(iPho->superCluster()->preshowerEnergyPlane1());
    phoESEnP2_        .push_back(iPho->superCluster()->preshowerEnergyPlane2());

    phoSCE_           .push_back(iPho->superCluster()->energy());
    phoSCRawE_        .push_back(iPho->superCluster()->rawEnergy());
    phoSCEta_         .push_back(iPho->superCluster()->eta());
    phoSCPhi_         .push_back(iPho->superCluster()->phi());
    phoSCEtaWidth_    .push_back(iPho->superCluster()->etaWidth());
    phoSCPhiWidth_    .push_back(iPho->superCluster()->phiWidth());
    phoSCBrem_        .push_back(iPho->superCluster()->phiWidth()/iPho->superCluster()->etaWidth());

    UChar_t _phoQualityBits = 0;
    if(iPho->hasPixelSeed()) setbit(_phoQualityBits, 0);
    if(iPho->passElectronVeto()) setbit(_phoQualityBits, 1);
    phoQualityBits_.push_back(_phoQualityBits);

    phoR9_            .push_back(iPho->r9());
    phoHoverE_        .push_back(iPho->hadTowOverEm());
    phoESEffSigmaRR_  .push_back(lazyTool.eseffsirir(*(iPho->superCluster())));

    /*phoPFChIso_       .push_back(iPho->userFloat("phoChargedIsolation"));
    phoPFPhoIso_      .push_back(iPho->userFloat("phoPhotonIsolation"));
    phoPFNeuIso_      .push_back(iPho->userFloat("phoNeutralHadronIsolation"));
    phoPFChWorstIso_  .push_back(iPho->userFloat("phoWorstChargedIsolation"));
    phoIDMVA_         .push_back(iPho->userFloat("PhotonMVAEstimatorRunIIFall17v2Values"));
    */


    const auto pho = photonHandle->ptrAt(nPho_);
    phoPFChIso_              .push_back((*phoChargedIsolationMap)[pho->originalObjectRef()]);
    phoPFPhoIso_             .push_back((*phoPhotonIsolationMap)[pho->originalObjectRef()]);
    phoPFNeuIso_             .push_back((*phoNeutralHadronIsolationMap)[pho->originalObjectRef()]);
    phoPFChWorstIso_         .push_back((*phoWorstChargedIsolationMap)[pho->originalObjectRef()]);

    
    // VID decisions
    UShort_t tmpphoIDbit = 0;
    //cout<<"Photons "<<endl;
    bool isPassLoose  = (*loose_id_decisions)[pho->originalObjectRef()];
    if(isPassLoose) setbit(tmpphoIDbit, 0);
    //cout<<"isPassLoose "<<isPassLoose<<endl;

    bool isPassMedium = (*medium_id_decisions)[pho->originalObjectRef()];
    if(isPassMedium) setbit(tmpphoIDbit, 1);
    //cout<<"isPassMedium "<<isPassMedium<<endl;

    bool isPassTight  = (*tight_id_decisions)[pho->originalObjectRef()];
    if(isPassTight) setbit(tmpphoIDbit, 2);
    //cout<<"isPassTight "<<isPassTight<<endl;

    if (!(iPho->mipIsHalo()))  setbit(tmpphoIDbit, 7);
    
    phoIDbit_.push_back(tmpphoIDbit);    
    
    phoIDMVA_.push_back((*mvaValues)[pho->originalObjectRef()]);


    phoSeedBCE_        .push_back(iPho->superCluster()->seed()->energy());
    phoSeedBCEta_      .push_back(iPho->superCluster()->seed()->eta());
    // phoSeedTimeFull5x5_.push_back(lazyToolnoZS.SuperClusterSeedTime(*(iPho->superCluster())));
    phoMIPChi2_        .push_back(iPho->mipChi2());
    phoMIPTotEnergy_   .push_back(iPho->mipTotEnergy());
    phoMIPSlope_       .push_back(iPho->mipSlope());
    phoMIPIntercept_   .push_back(iPho->mipIntercept());
    phoMIPNhitCone_    .push_back(iPho->mipNhitCone());


    
    // get photon supercluster index (for looking up from the SC branches)
    /*
    if(ecalSChandle.isValid()){
      const reco::SuperCluster * _tmpPhoSC = (iPho->superCluster().isAvailable()) ? iPho->superCluster().get() : nullptr;
      Short_t tmpPhoSCindex = (_tmpPhoSC == nullptr) ? -999 : std::distance(ecalSChandle->begin(), (std::vector<reco::SuperCluster>::const_iterator) _tmpPhoSC);
      phoSCindex_.push_back(tmpPhoSCindex);
      //check
      // if(tmpPhoSCindex>=0){
      //   std::cout<<_tmpPhoSC->eta()<<" "<<(ecalSChandle->begin()+tmpPhoSCindex)->eta()<<" "<<_tmpPhoSC->phi()<<" "<<(ecalSChandle->begin()+tmpPhoSCindex)->phi()<<std::endl;
      // }
    }
    */
    
    UChar_t tmpphoFiducialRegion = 0;
    if(iPho->isEB()) setbit(tmpphoFiducialRegion, 0);
    if(iPho->isEE()) setbit(tmpphoFiducialRegion, 1);
    if(iPho->isEBEEGap()) setbit(tmpphoFiducialRegion, 2);
    if(iPho->isEBEtaGap ()) setbit(tmpphoFiducialRegion, 3);
    if(iPho->isEBPhiGap ()) setbit(tmpphoFiducialRegion, 4);
    if(iPho->isEEDeeGap ()) setbit(tmpphoFiducialRegion, 5);
    if(iPho->isEERingGap()) setbit(tmpphoFiducialRegion, 6);
    phoFiducialRegion_  .push_back(tmpphoFiducialRegion);

    /*
    // VID decisions
    UShort_t tmpphoIDbit = 0;
    // if(year_ == 2017){
    //// https://twiki.cern.ch/twiki/bin/view/CMS/EgammaRunIIRecommendations?rev=9#Fall17v2_AN1
    bool isPassLoose  = iPho->photonID("cutBasedPhotonID-Fall17-94X-V2-loose");
    if (isPassLoose)  setbit(tmpphoIDbit, 0);
    bool isPassMedium = iPho->photonID("cutBasedPhotonID-Fall17-94X-V2-medium");
    if (isPassMedium) setbit(tmpphoIDbit, 1);
    bool isPassTight  = iPho->photonID("cutBasedPhotonID-Fall17-94X-V2-tight");
    if (isPassTight)  setbit(tmpphoIDbit, 2);
    bool isPassMVAv2wp80  = iPho->photonID("mvaPhoID-RunIIFall17-v2-wp80");
    if (isPassMVAv2wp80)  setbit(tmpphoIDbit, 3);
    bool isPassMVAv2wp90  = iPho->photonID("mvaPhoID-RunIIFall17-v2-wp90");
    if (isPassMVAv2wp90)  setbit(tmpphoIDbit, 4);
    // bool isPassMVAisov2wp80  = iPho->photonID("mvaPhoID-Fall17-iso-V2-wp80");
    // if (isPassMVAisov2wp80)  setbit(tmpphoIDbit, 5);
    // bool isPassMVAisov2wp90  = iPho->photonID("mvaPhoID-Fall17-iso-V2-wp90");
    // if (isPassMVAisov2wp90)  setbit(tmpphoIDbit, 6);

    if (!(iPho->mipIsHalo()))  setbit(tmpphoIDbit, 7);
    phoIDbit_.push_back(tmpphoIDbit);



    // systematics for energy scale and resolution
    phoScale_stat_up_.push_back(iPho->userFloat("energyScaleStatUp"));
    phoScale_stat_dn_.push_back(iPho->userFloat("energyScaleStatDown"));
    phoScale_syst_up_.push_back(iPho->userFloat("energyScaleSystUp"));
    phoScale_syst_dn_.push_back(iPho->userFloat("energyScaleSystDown"));
    phoScale_gain_up_.push_back(iPho->userFloat("energyScaleGainUp"));
    phoScale_gain_dn_.push_back(iPho->userFloat("energyScaleGainDown"));
    phoResol_rho_up_ .push_back(iPho->userFloat("energySigmaRhoUp"));
    phoResol_rho_dn_ .push_back(iPho->userFloat("energySigmaRhoDown"));
    phoResol_phi_up_ .push_back(iPho->userFloat("energySigmaPhiUp"));
    phoResol_phi_dn_ .push_back(iPho->userFloat("energySigmaPhiDown"));
    */

    ///////////////////////////////SATURATED/UNSATURATED ///from ggFlash////
    DetId seed = (iPho->superCluster()->seed()->hitsAndFractions())[0].first;
    bool isBarrel = (seed.subdetId() == EcalBarrel);
    const EcalRecHitCollection * rechits = (isBarrel ? lazyTool.getEcalEBRecHitCollection() : lazyTool.getEcalEERecHitCollection());

    EcalRecHitCollection::const_iterator theSeedHit = rechits->find(seed);
    if (theSeedHit != rechits->end()) {
      phoSeedTime_.push_back(theSeedHit->time());
      phoSeedEnergy_.push_back(theSeedHit->energy());
    } else{
      phoSeedTime_  .push_back(-99.);
      phoSeedEnergy_.push_back(-99.);
    }

    phoFiredSingleTrgs_     .push_back(matchSinglePhotonTriggerFilters(iPho->et(), iPho->eta(), iPho->phi()));
    phoFiredDoubleTrgs_     .push_back(matchDoublePhotonTriggerFilters(iPho->et(), iPho->eta(), iPho->phi()));
    phoFiredTripleTrgs_     .push_back(matchTriplePhotonTriggerFilters(iPho->et(), iPho->eta(), iPho->phi()));
    phoFiredL1Trgs_         .push_back(matchL1TriggerFilters(iPho->et(), iPho->eta(), iPho->phi()));

    std::vector<Float_t> vCov = lazyToolnoZS.localCovariances( *(iPho->superCluster()->seed()) );
    const Float_t spp = (isnan(vCov[2]) ? 0. : sqrt(vCov[2]));
    const Float_t sep = vCov[1];

    phoSigmaIEtaIEtaFull5x5_ .push_back(iPho->full5x5_sigmaIetaIeta());
    phoSigmaIEtaIPhiFull5x5_ .push_back(sep);
    phoSigmaIPhiIPhiFull5x5_ .push_back(spp);
    phoE2x2Full5x5_          .push_back(lazyToolnoZS.e2x2(*(iPho->superCluster()->seed())));
    phoE5x5Full5x5_          .push_back(iPho->full5x5_e5x5());
    phoR9Full5x5_            .push_back(iPho->full5x5_r9());


    if(doGenParticles_){
      const reco::GenParticle * phoGen_ = iPho->genParticle(); // I don't know what matching algoritm is used - https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2017#MC_Truth
      Short_t phoGenPos_ = (phoGen_ == nullptr) ? -999 : std::distance(genParticlesHandle->begin(), (std::vector<reco::GenParticle>::const_iterator) phoGen_);
      // if(phoGenPos_>=0) std::cout<<"pho->genParticle->pdgID "<<phoGen_->pdgId()<<" prunedGenParticle("<<phoGenPos_<<")->pdgid() "<< (&*genParticlesHandle->begin()+phoGenPos_)->pdgId()<<std::endl;
      pho_gen_index_.push_back(phoGenPos_);
    }

    nPho_++;
  }
}






/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////         OOT Photons     ////////////////////////////////////////////////////////
UShort_t          nootPho_;
vector<Float_t>  ootPho_E_;
vector<Float_t>  ootPho_Et_;
vector<Float_t>  ootPho_Eta_;
vector<Float_t>  ootPho_Phi_;
std::vector<Short_t> ootPho_SCindex_;
vector<UChar_t> ootPho_FiducialRegion_;
vector<UChar_t>   ootPho_QualityBits_;
vector<Float_t>  ootPho_R9_;
vector<Float_t>  ootPho_HoverE_;
vector<Float_t>  ootPho_ESEffSigmaRR_;
vector<Float_t>  ootPho_SigmaIEtaIEtaFull5x5_;
vector<Float_t>  ootPho_SigmaIEtaIPhiFull5x5_;
vector<Float_t>  ootPho_SigmaIPhiIPhiFull5x5_;
vector<Float_t>  ootPho_R9Full5x5_;
// vector<Float_t>  ootPho_R9noZS_;
// vector<Float_t>  ootPho_sigmaIetaIeta_NoZS_;
// vector<Float_t>  ootPho_sigmaIetaIphi_NoZS_;
// vector<Float_t>  ootPho_sigmaIphiIphi_NoZS_;
vector<ULong64_t> ootPho_FiredSingleTrgs_;
vector<ULong64_t> ootPho_FiredDoubleTrgs_;
vector<ULong64_t> ootPho_FiredTripleTrgs_;
vector<ULong64_t> ootPho_FiredL1Trgs_;
vector<Float_t>  ootPho_SeedTime_;
vector<Float_t>  ootPho_SeedEnergy_;
vector<Float_t>  ootPho_MIPChi2_;
vector<Float_t>  ootPho_MIPTotEnergy_;
vector<Float_t>  ootPho_MIPSlope_;
vector<Float_t>  ootPho_MIPIntercept_;
vector<Short_t>  ootPho_MIPNhitCone_;
vector<UChar_t>   ootPho_IDbit_;


void ggNtuplizer::branchesPhotonsOOT(TTree* tree) {
  tree->Branch("nootPho",                     &nootPho_);
  tree->Branch("ootPho_E",                    &ootPho_E_);
  tree->Branch("ootPho_Et",                   &ootPho_Et_);
  tree->Branch("ootPho_Eta",                  &ootPho_Eta_);
  tree->Branch("ootPho_Phi",                  &ootPho_Phi_);
  tree->Branch("ootPho_SCindex",              &ootPho_SCindex_);
  tree->Branch("ootPho_FiducialRegion",       &ootPho_FiducialRegion_);
  tree->Branch("ootPho_QualityBits",          &ootPho_QualityBits_);
  tree->Branch("ootPho_R9",                   &ootPho_R9_);
  tree->Branch("ootPho_HoverE",               &ootPho_HoverE_);
  tree->Branch("ootPho_ESEffSigmaRR",         &ootPho_ESEffSigmaRR_);
  tree->Branch("ootPho_SigmaIEtaIEtaFull5x5", &ootPho_SigmaIEtaIEtaFull5x5_);
  tree->Branch("ootPho_SigmaIEtaIPhiFull5x5", &ootPho_SigmaIEtaIPhiFull5x5_);
  tree->Branch("ootPho_SigmaIPhiIPhiFull5x5", &ootPho_SigmaIPhiIPhiFull5x5_);
  tree->Branch("ootPho_R9Full5x5",            &ootPho_R9Full5x5_);
  // tree->Branch("ootPho_R9noZS",             &ootPho_R9noZS_);
  // tree->Branch("ootPho_sigmaIetaIeta_NoZS",              &ootPho_sigmaIetaIeta_NoZS_);
  // tree->Branch("ootPho_sigmaIetaIphi_NoZS",             &ootPho_sigmaIetaIphi_NoZS_);
  // tree->Branch("ootPho_sigmaIphiIphi_NoZS",         &ootPho_sigmaIphiIphi_NoZS_);
  tree->Branch("ootPho_FiredSingleTrgs",      &ootPho_FiredSingleTrgs_);
  tree->Branch("ootPho_FiredDoubleTrgs",      &ootPho_FiredDoubleTrgs_);
  tree->Branch("ootPho_FiredTripleTrgs",      &ootPho_FiredTripleTrgs_);
  tree->Branch("ootPho_FiredL1Trgs",          &ootPho_FiredL1Trgs_);
  tree->Branch("ootPho_SeedTime",             &ootPho_SeedTime_);
  tree->Branch("ootPho_SeedEnergy",           &ootPho_SeedEnergy_);
  tree->Branch("ootPho_MIPChi2",              &ootPho_MIPChi2_);
  tree->Branch("ootPho_MIPTotEnergy",         &ootPho_MIPTotEnergy_);
  tree->Branch("ootPho_MIPSlope",             &ootPho_MIPSlope_);
  tree->Branch("ootPho_MIPIntercept",         &ootPho_MIPIntercept_);
  tree->Branch("ootPho_MIPNhitCone",          &ootPho_MIPNhitCone_);
  tree->Branch("ootPho_IDbit",                &ootPho_IDbit_);
}


void ggNtuplizer::fillPhotonsOOT(const edm::Event& e, const edm::EventSetup& es) {

  // cleanup from previous execution
  ootPho_E_                   .clear();
  ootPho_Et_                  .clear();
  ootPho_Eta_                 .clear();
  ootPho_Phi_                 .clear();
  ootPho_SCindex_             .clear();
  ootPho_FiducialRegion_      .clear();
  ootPho_QualityBits_.clear();
  ootPho_R9_                  .clear();
  ootPho_HoverE_              .clear();
  ootPho_ESEffSigmaRR_        .clear();
  ootPho_SigmaIEtaIEtaFull5x5_.clear();
  ootPho_SigmaIEtaIPhiFull5x5_.clear();
  ootPho_SigmaIPhiIPhiFull5x5_.clear();
  ootPho_R9Full5x5_           .clear();
  // ootPho_R9noZS_             .clear();
  // ootPho_sigmaIetaIeta_NoZS_            .clear();
  // ootPho_sigmaIetaIphi_NoZS_            .clear();
  // ootPho_sigmaIphiIphi_NoZS_        .clear();
  ootPho_FiredSingleTrgs_     .clear();
  ootPho_FiredDoubleTrgs_     .clear();
  ootPho_FiredTripleTrgs_     .clear();
  ootPho_FiredL1Trgs_         .clear();
  ootPho_SeedTime_            .clear();
  ootPho_SeedEnergy_          .clear();
  ootPho_MIPChi2_           .clear();
  ootPho_MIPTotEnergy_      .clear();
  ootPho_MIPSlope_          .clear();
  ootPho_MIPIntercept_      .clear();
  ootPho_MIPNhitCone_       .clear();
  ootPho_IDbit_        .clear();
  nootPho_ = 0;

  edm::Handle<edm::View<pat::Photon> > photonOOT_Handle;
  e.getByToken(photonOOTCollection_, photonOOT_Handle);

  if (!photonOOT_Handle.isValid()) {
    edm::LogWarning("ggNtuplizer") << "no OOT pat::Photons in event";
    return;
  }

  edm::Handle<std::vector<reco::SuperCluster>> ecalSC_OOT_handle;
  e.getByToken(ecalSC_OOT_collection_, ecalSC_OOT_handle);

  EcalClusterLazyTools       lazyTool    (e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);
  noZS::EcalClusterLazyTools lazyToolnoZS(e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);

  for (edm::View<pat::Photon>::const_iterator iootPho_ = photonOOT_Handle->begin(); iootPho_ != photonOOT_Handle->end(); ++iootPho_){

    // if(iootPho_->et() < 15.) continue;

    ootPho_E_             .push_back(iootPho_->energy());
    ootPho_Et_            .push_back(iootPho_->et());
    ootPho_Eta_           .push_back(iootPho_->eta());
    ootPho_Phi_           .push_back(iootPho_->phi());

    UChar_t _ootPho_QualityBits = 0;
    if(iootPho_->hasPixelSeed()) setbit(_ootPho_QualityBits, 0);
    if(iootPho_->passElectronVeto()) setbit(_ootPho_QualityBits, 1);
    ootPho_QualityBits_.push_back(_ootPho_QualityBits);

    ootPho_R9_            .push_back(iootPho_->r9());
    ootPho_HoverE_        .push_back(iootPho_->hadTowOverEm());
    ootPho_ESEffSigmaRR_  .push_back(lazyTool.eseffsirir(*(iootPho_->superCluster())));
    // ootPho_R9noZS_       .push_back(iootPho_->userFloat("r9_NoZS"));
    // ootPho_sigmaIetaIeta_NoZS_      .push_back(iootPho_->userFloat("sigmaIetaIeta_NoZS"));
    // ootPho_sigmaIetaIphi_NoZS_      .push_back(iootPho_->userFloat("sigmaIetaIphi_NoZS"));
    // ootPho_sigmaIphiIphi_NoZS_  .push_back(iootPho_->userFloat("sigmaIphiIphi_NoZS"));
    ootPho_MIPChi2_        .push_back(iootPho_->mipChi2());
    ootPho_MIPTotEnergy_   .push_back(iootPho_->mipTotEnergy());
    ootPho_MIPSlope_       .push_back(iootPho_->mipSlope());
    ootPho_MIPIntercept_   .push_back(iootPho_->mipIntercept());
    ootPho_MIPNhitCone_    .push_back(iootPho_->mipNhitCone());


    // get ootPho_ton supercluster index (for looking up from the SC branches)
    if(ecalSC_OOT_handle.isValid()){
      const reco::SuperCluster * _tmpootPho_SC = (iootPho_->superCluster().isAvailable()) ? iootPho_->superCluster().get() : nullptr;
      Short_t tmpootPho_SCindex = (_tmpootPho_SC == nullptr) ? -999 : std::distance(ecalSC_OOT_handle->begin(), (std::vector<reco::SuperCluster>::const_iterator) _tmpootPho_SC);
      ootPho_SCindex_.push_back(tmpootPho_SCindex);
      //check
      // if(tmpootPho_SCindex>=0){
      //   std::cout<<_tmpootPho_SC->eta()<<" "<<(ecalSChandle->begin()+tmpootPho_SCindex)->eta()<<" "<<_tmpootPho_SC->phi()<<" "<<(ecalSChandle->begin()+tmpootPho_SCindex)->phi()<<std::endl;
      // }
    }

    UChar_t tmpootPho_FiducialRegion = 0;
    if(iootPho_->isEB()) setbit(tmpootPho_FiducialRegion, 0);
    if(iootPho_->isEE()) setbit(tmpootPho_FiducialRegion, 1);
    if(iootPho_->isEBEEGap()) setbit(tmpootPho_FiducialRegion, 2);
    if(iootPho_->isEBEtaGap ()) setbit(tmpootPho_FiducialRegion, 3);
    if(iootPho_->isEBPhiGap ()) setbit(tmpootPho_FiducialRegion, 4);
    if(iootPho_->isEEDeeGap ()) setbit(tmpootPho_FiducialRegion, 5);
    if(iootPho_->isEERingGap()) setbit(tmpootPho_FiducialRegion, 6);
    ootPho_FiducialRegion_  .push_back(tmpootPho_FiducialRegion);


    // VID decisions
    UShort_t tmpootPho_IDbit = 0;
    // if(year_ == 2017){
    //// https://twiki.cern.ch/twiki/bin/view/CMS/EgammaRunIIRecommendations?rev=9#Fall17v2_AN1
    // bool isPassLoose  = iootPho_->photonID("cutBasedphotonID-Fall17-94X-V2-loose");
    // if (isPassLoose)  setbit(tmpootPho_IDbit, 0);
    // bool isPassMedium = iootPho_->photonID("cutBasedphotonID-Fall17-94X-V2-medium");
    // if (isPassMedium) setbit(tmpootPho_IDbit, 1);
    // bool isPassTight  = iootPho_->photonID("cutBasedphotonID-Fall17-94X-V2-tight");
    // if (isPassTight)  setbit(tmpootPho_IDbit, 2);
    // bool isPassMVAv2wp80  = iootPho_->photonID("mvaootPho_ID-RunIIFall17-v2-wp80");
    // if (isPassMVAv2wp80)  setbit(tmpootPho_IDbit, 3);
    // bool isPassMVAv2wp90  = iootPho_->photonID("mvaootPho_ID-RunIIFall17-v2-wp90");
    // if (isPassMVAv2wp90)  setbit(tmpootPho_IDbit, 4);
    // bool isPassMVAisov2wp80  = iootPho_->photonID("mvaootPho_ID-Fall17-iso-V2-wp80");
    // if (isPassMVAisov2wp80)  setbit(tmpootPho_IDbit, 5);
    // bool isPassMVAisov2wp90  = iootPho_->photonID("mvaootPho_ID-Fall17-iso-V2-wp90");
    // if (isPassMVAisov2wp90)  setbit(tmpootPho_IDbit, 6);
    if (!(iootPho_->mipIsHalo()))  setbit(tmpootPho_IDbit, 7);

    ootPho_IDbit_.push_back(tmpootPho_IDbit);


    ///////////////////////////////SATURATED/UNSATURATED ///from ggFlash////
    DetId seed = (iootPho_->superCluster()->seed()->hitsAndFractions())[0].first;
    bool isBarrel = (seed.subdetId() == EcalBarrel);
    const EcalRecHitCollection * rechits = isBarrel ? lazyTool.getEcalEBRecHitCollection() : lazyTool.getEcalEERecHitCollection();

    EcalRecHitCollection::const_iterator theSeedHit = rechits->find(seed);
    if (theSeedHit != rechits->end()) {
      ootPho_SeedTime_.push_back(theSeedHit->time());
      ootPho_SeedEnergy_.push_back(theSeedHit->energy());
    } else{
      ootPho_SeedTime_  .push_back(-99.);
      ootPho_SeedEnergy_.push_back(-99.);
    }

    ootPho_FiredSingleTrgs_     .push_back(matchSinglePhotonTriggerFilters(iootPho_->et(), iootPho_->eta(), iootPho_->phi()));
    ootPho_FiredDoubleTrgs_     .push_back(matchDoublePhotonTriggerFilters(iootPho_->et(), iootPho_->eta(), iootPho_->phi()));
    ootPho_FiredTripleTrgs_     .push_back(matchTriplePhotonTriggerFilters(iootPho_->et(), iootPho_->eta(), iootPho_->phi()));
    ootPho_FiredL1Trgs_         .push_back(matchL1TriggerFilters(iootPho_->et(), iootPho_->eta(), iootPho_->phi()));

    std::vector<Float_t> vCov = lazyToolnoZS.localCovariances( *(iootPho_->superCluster()->seed()) );
    const Float_t spp = (isnan(vCov[2]) ? 0. : sqrt(vCov[2]));
    const Float_t sep = vCov[1];

    ootPho_SigmaIEtaIEtaFull5x5_ .push_back(iootPho_->full5x5_sigmaIetaIeta());
    ootPho_SigmaIEtaIPhiFull5x5_ .push_back(sep);
    ootPho_SigmaIPhiIPhiFull5x5_ .push_back(spp);
    ootPho_R9Full5x5_            .push_back(iootPho_->full5x5_r9());

    nootPho_++;
  }
}
