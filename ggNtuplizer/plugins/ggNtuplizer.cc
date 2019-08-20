#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"

using namespace std;
using namespace edm;

void ggNtuplizer::endJob() {
  tree_->BuildIndex("run", "event");
  std::cout<<"\tTTree "<<tree_->GetName()<<" indexed on run and event."<<std::endl<<" Complete!"<<std::endl;
};


ggNtuplizer::ggNtuplizer(const edm::ParameterSet& ps) :
hltPrescaleProvider_(ps, consumesCollector(), *this)
{

  getECALprefiringWeights_      = ps.getParameter<bool>("getECALprefiringWeights");
  if(getECALprefiringWeights_){
    prefweight_token              = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProb"));
    prefweightup_token            = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbUp"));
    prefweightdown_token          = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbDown"));
  }


  development_                = ps.getParameter<bool>("development");
  addFilterInfoAOD_           = ps.getParameter<bool>("addFilterInfoAOD");
  addFilterInfoMINIAOD_       = ps.getParameter<bool>("addFilterInfoMINIAOD");
  doNoHFMET_                  = ps.getParameter<bool>("doNoHFMET");


  doOOTphotons_               = ps.getParameter<bool>("doOOTphotons");
  doGenParticles_             = ps.getParameter<bool>("doGenParticles");
  runOnParticleGun_           = ps.getParameter<bool>("runOnParticleGun");
  runOnSherpa_                = ps.getParameter<bool>("runOnSherpa");
  dumpPFPhotons_              = ps.getParameter<bool>("dumpPFPhotons");
  dumpJets_                   = ps.getParameter<bool>("dumpJets");
  dumpAK8Jets_                = ps.getParameter<bool>("dumpAK8Jets");
  dumpSoftDrop_               = ps.getParameter<bool>("dumpSoftDrop");
  dumpTaus_                   = ps.getParameter<bool>("dumpTaus");
  dumpPDFSystWeight_          = ps.getParameter<bool>("dumpPDFSystWeight");
  dumpHFElectrons_            = ps.getParameter<bool>("dumpHFElectrons");
  doRecHits_                  = ps.getParameter<bool>("doRecHits");

  year_                       = ps.getParameter<int>("year");

  trgFilterDeltaPtCut_        = ps.getParameter<double>("trgFilterDeltaPtCut");
  trgFilterDeltaRCut_         = ps.getParameter<double>("trgFilterDeltaRCut");

  beamHaloSummaryToken_       = consumes<reco::BeamHaloSummary>         (ps.getParameter<InputTag>("beamHaloSummary"));
  vtxLabel_                   = consumes<reco::VertexCollection>        (ps.getParameter<InputTag>("VtxLabel"));
  vtxBSLabel_                 = consumes<reco::VertexCollection>        (ps.getParameter<InputTag>("VtxBSLabel"));
  rhoLabel_                   = consumes<double>                        (ps.getParameter<InputTag>("rhoLabel"));
  rhoCentralLabel_            = consumes<double>                        (ps.getParameter<InputTag>("rhoCentralLabel"));
  // trgEventLabel_              = consumes<trigger::TriggerEvent>         (ps.getParameter<InputTag>("triggerEvent"));
  triggerObjectsLabel_        = consumes<pat::TriggerObjectStandAloneCollection>(ps.getParameter<edm::InputTag>("triggerEvent"));
  trgResultsLabel_            = consumes<edm::TriggerResults>           (ps.getParameter<InputTag>("triggerResults"));
  patTrgResultsLabel_         = consumes<edm::TriggerResults>           (ps.getParameter<InputTag>("patTriggerResults"));
  trgResultsProcess_          =                                          ps.getParameter<InputTag>("triggerResults").process();
  generatorLabel_             = consumes<GenEventInfoProduct>           (ps.getParameter<InputTag>("generatorLabel"));
  lheEventLabel_              = consumes<LHEEventProduct>               (ps.getParameter<InputTag>("LHEEventLabel"));
  puCollection_               = consumes<vector<PileupSummaryInfo> >    (ps.getParameter<InputTag>("pileupCollection"));
  genParticlesCollection_     = consumes<vector<reco::GenParticle> >    (ps.getParameter<InputTag>("genParticleSrc"));
  pfMETlabel_                 = consumes<View<pat::MET> >               (ps.getParameter<InputTag>("pfMETLabel"));
  electronCollection_         = consumes<View<pat::Electron> >          (ps.getParameter<InputTag>("electronSrc"));
  gsfTracks_                  = consumes<View<reco::GsfTrack>>          (ps.getParameter<InputTag>("gsfTrackSrc"));
  hbheRecHitCollection_                 = consumes<HBHERecHitCollection >               (ps.getParameter<InputTag>("hbheRecHitCollection"));
  hfRecHitCollection_                 = consumes<HFRecHitCollection >               (ps.getParameter<InputTag>("hfRecHitCollection"));
  hoRecHitCollection_                 = consumes<HORecHitCollection >               (ps.getParameter<InputTag>("hoRecHitCollection"));

  cscSegmentsCollection_                = consumes<CSCSegmentCollection>(ps.getParameter<edm::InputTag>("cscSegmentsCollection"));

  //BadChCandFilterToken_       = consumes<bool>                          (ps.getParameter<InputTag>("BadChargedCandidateFilter"));
  //BadPFMuonFilterToken_       = consumes<bool>                          (ps.getParameter<edm::InputTag>("BadPFMuonFilter"));

  photonCollection_           = consumes<View<pat::Photon> >            (ps.getParameter<InputTag>("photonSrc"));
  photonOOTCollection_        = consumes<View<pat::Photon> >            (ps.getParameter<InputTag>("photonOOTSrc"));
  muonCollection_             = consumes<View<pat::Muon> >              (ps.getParameter<InputTag>("muonSrc"));
  ebReducedRecHitCollection_  = consumes<EcalRecHitCollection>          (ps.getParameter<InputTag>("ebReducedRecHitCollection"));
  eeReducedRecHitCollection_  = consumes<EcalRecHitCollection>          (ps.getParameter<InputTag>("eeReducedRecHitCollection"));
  esReducedRecHitCollection_  = consumes<EcalRecHitCollection>          (ps.getParameter<InputTag>("esReducedRecHitCollection"));

  ecalSC_OOT_collection_      = consumes<std::vector<reco::SuperCluster>>(ps.getParameter<InputTag>("ecalSCOOTcollection"));
  ecalSCcollection_           = consumes<std::vector<reco::SuperCluster>>(ps.getParameter<InputTag>("ecalSCcollection"));

  //ecalSCcollectionEB_           = consumes<std::vector<reco::SuperCluster>>(ps.getParameter<InputTag>("ecalSCcollectionEB"));
  //ecalSCcollectionEE_           = consumes<std::vector<reco::SuperCluster>>(ps.getParameter<InputTag>("ecalSCcollectionEE"));
  //ecalSC_OOT_collectionEB_      = consumes<std::vector<reco::SuperCluster>>(ps.getParameter<InputTag>("ecalSCOOTcollectionEB"));
  //ecalSC_OOT_collectionEE_      = consumes<std::vector<reco::SuperCluster>>(ps.getParameter<InputTag>("ecalSCOOTcollectionEE"));

  recophotonCollection_       = consumes<reco::PhotonCollection>        (ps.getParameter<InputTag>("recoPhotonSrc"));
  tracklabel_                 = consumes<reco::TrackCollection>         (ps.getParameter<InputTag>("TrackLabel"));
  gsfElectronlabel_           = consumes<reco::GsfElectronCollection>   (ps.getParameter<InputTag>("gsfElectronLabel"));
  tauCollection_              = consumes<vector<pat::Tau> >             (ps.getParameter<InputTag>("tauSrc"));
  pfAllParticles_             = consumes<reco::PFCandidateCollection>   (ps.getParameter<InputTag>("PFAllCandidates"));
  pckPFCandidateCollection_   = consumes<pat::PackedCandidateCollection>(ps.getParameter<InputTag>("packedPFCands"));
  pckPFCdsLabel_              = consumes<vector<pat::PackedCandidate>>  (ps.getParameter<InputTag>("packedPFCands"));
  recoCdsLabel_               = consumes<View<reco::Candidate>>         (ps.getParameter<InputTag>("packedPFCands"));

  ak4PFJetsCHSLabel_          = consumes<View<pat::Jet> >               (ps.getParameter<InputTag>("ak4PFJetsCHSSrc"));
  ak4PFJetsPUPPILabel_        = consumes<View<pat::Jet> >               (ps.getParameter<InputTag>("ak4PFJetsPUPPISrc"));
  ak8JetsPUPPILabel_          = consumes<View<pat::Jet> >               (ps.getParameter<InputTag>("ak8JetsPUPPISrc"));

  ak4PFJetsCHSGenJetLabel_    = consumes<std::vector<reco::GenJet> >    (ps.getParameter<InputTag>("ak4PFJetsCHSGenJetLabel"));
  ak8GenJetLabel_             = consumes<std::vector<reco::GenJet> >    (ps.getParameter<InputTag>("ak8GenJetLabel"));
  newparticles_               =                                          ps.getParameter< vector<int > >("newParticles");
  //ecalBadCalibFilterUpdateToken_ = consumes< Bool_t >(ps.getParameter<InputTag>("ecalBadCalibFilter"));

  /////////PHOTON AND ELECTRON ID
  // electron ID 
  eleVetoIdMapToken_    = consumes<edm::ValueMap<bool> >(ps.getParameter<edm::InputTag>("eleVetoIdMap"));
  eleLooseIdMapToken_   = consumes<edm::ValueMap<bool> >(ps.getParameter<edm::InputTag>("eleLooseIdMap"));
  eleMediumIdMapToken_  = consumes<edm::ValueMap<bool> >(ps.getParameter<edm::InputTag>("eleMediumIdMap"));
  eleTightIdMapToken_   = consumes<edm::ValueMap<bool> >(ps.getParameter<edm::InputTag>("eleTightIdMap"));
  eleHEEPIdMapToken_    = consumes<edm::ValueMap<bool> >(ps.getParameter<edm::InputTag>("eleHEEPIdMap"));
  eleMVAValuesMapToken_ = consumes<edm::ValueMap<float> >(ps.getParameter<edm::InputTag>("eleMVAValuesMap"));

  // Photon ID in VID framwork 
  phoLooseIdMapToken_             = consumes<edm::ValueMap<bool> >(ps.getParameter<edm::InputTag>("phoLooseIdMap"));
  phoMediumIdMapToken_            = consumes<edm::ValueMap<bool> >(ps.getParameter<edm::InputTag>("phoMediumIdMap"));
  phoTightIdMapToken_             = consumes<edm::ValueMap<bool> >(ps.getParameter<edm::InputTag>("phoTightIdMap"));
  phoMVAValuesMapToken_           = consumes<edm::ValueMap<float> >(ps.getParameter<edm::InputTag>("phoMVAValuesMap")); 
  phoChargedIsolationToken_       = consumes <edm::ValueMap<float> >(ps.getParameter<edm::InputTag>("phoChargedIsolation"));
  phoNeutralHadronIsolationToken_ = consumes <edm::ValueMap<float> >(ps.getParameter<edm::InputTag>("phoNeutralHadronIsolation"));
  phoPhotonIsolationToken_        = consumes <edm::ValueMap<float> >(ps.getParameter<edm::InputTag>("phoPhotonIsolation"));
  phoWorstChargedIsolationToken_  = consumes <edm::ValueMap<float> >(ps.getParameter<edm::InputTag>("phoWorstChargedIsolation"));
  /////////END of PHOTON AND ELECTRON ID


  Service<TFileService> fs;
  tree_    = fs->make<TTree>("EventTree", "Event data");
  hEvents_ = fs->make<TH1F>("hEvents",    "total processed and skimmed events",   2,  0,   1);

  branchesGlobalEvent(tree_);

  if (doGenParticles_) {
    branchesGenInfo(tree_, fs);
    branchesGenPart(tree_);
  }
  if(dumpJets_ && doGenParticles_) branchesGenAK4JetPart(tree_);
  if(dumpAK8Jets_ && doGenParticles_) branchesGenAK8JetPart(tree_);
  branchesMET(tree_);
  branchesPhotons(tree_);
  //branchesECALSC(tree_);
  if(doOOTphotons_) {
    branchesPhotonsOOT(tree_);
    branchesECALOOTSC(tree_);
  }
  branchesElectrons(tree_);
  branchesMuons(tree_);
  if (dumpJets_)        branchesAK4CHSJets(tree_);
  if (dumpAK8Jets_)     branchesAK8PUPPIJets(tree_);

  if(doRecHits_) branchesRecHit(tree_);

  //debug = true;
  debug = false;
}

ggNtuplizer::~ggNtuplizer() {
}

void ggNtuplizer::analyze(const edm::Event& e, const edm::EventSetup& es) {



  hEvents_->Fill(0.2);

  if(debug) std::cout<<"taking vtx collection"<<std::endl;
  edm::Handle<reco::VertexCollection> vtxHandle;
  e.getByToken(vtxLabel_, vtxHandle);

  if(debug) std::cout<<"GOT vtx collection"<<std::endl;

  reco::Vertex vtx;

  // best-known primary vertex coordinates
  math::XYZPoint pv(0, 0, 0);
  for (vector<reco::Vertex>::const_iterator v = vtxHandle->begin(); v != vtxHandle->end(); ++v) {
    // replace isFake() for miniAOD since it requires tracks while miniAOD vertices don't have tracks:
    // Vertex.h: bool isFake() const {return (chi2_==0 && ndof_==0 && tracks_.empty());}
    bool isFake = (v->chi2() == 0 && v->ndof() == 0);

    if (!isFake) {
      pv.SetXYZ(v->x(), v->y(), v->z());
      vtx = *v;
      break;
    }
  }

  if(debug) std::cout<<"Doing GlobalEvent"<<std::endl;

  //initTriggerFilters(e);
  fillGlobalEvent(e, es);

  if(debug) std::cout<<"Done GlobalEvent"<<std::endl;

  if (!e.isRealData() && doGenParticles_) {
    fillGenInfo(e);
    if (doGenParticles_)
      fillGenPart(e);
  }

  if(debug) std::cout<<"Doing MET"<<std::endl;
  fillMET(e, es);

  if(debug) std::cout<<"Doing Electrons"<<std::endl;
  fillElectrons(e, es, pv);

  if(debug) std::cout<<"Doing Muons"<<std::endl;
  fillMuons(e, pv, vtx);

  if(debug) std::cout<<"Doing Photons"<<std::endl;
  fillPhotons(e, es);

  if(debug) std::cout<<"Doing Rechits"<<std::endl;
  if(doRecHits_) fillRecHits(e, es);


  if(doOOTphotons_) {
    if(debug) std::cout<<"Doing OOTPhotons"<<std::endl;
    fillPhotonsOOT(e, es);
    //fillECALOOTSC(e, es);
  }


  if (dumpJets_){
    if(debug) std::cout<<"Doing Jets"<<std::endl;
    fillAK4CHSJets(e,es);
    if(debug) std::cout<<"Done Jets"<<std::endl;
  }

  if (dumpAK8Jets_){
    if(debug) std::cout<<"Dumping AK8"<<std::endl;
    fillAK8PUPPIJets(e,es);
    if(debug) std::cout<<"Done AK8"<<std::endl;
  }

  if(dumpJets_ && doGenParticles_){

    if(debug) std::cout<<"Doing AK4 and gen"<<std::endl;
    fillGenAK4JetInfo(e, vtx.z());
    if(debug) std::cout<<"Done AK4 and gen"<<std::endl;
  }

  if(dumpAK8Jets_ && doGenParticles_){
    if(debug) std::cout<<"Doing AK8 and gen"<<std::endl;
    fillGenAK8JetInfo(e, vtx.z());
    if(debug) std::cout<<"Done AK8 and gen"<<std::endl;
  }

  //fillECALSC(e, es);

  tree_->Fill();
  hEvents_->Fill(0.8);

}

// void ggNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
// {
//   //The following says we do not know what parameters are allowed so do no validation
//   // Please change this to state exactly what you do use, even if it is no parameters
//   edm::ParameterSetDescription desc;
//   desc.setUnknown();
//   descriptions.addDefault(desc);
// }

DEFINE_FWK_MODULE(ggNtuplizer);
