#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"
#include "ggAnalysis/ggNtuplizer/interface/GenParticleParentage.h"
using namespace std;

// (local) variables associated with tree branches
UShort_t        nGenAK4Jets_;
vector<float>  GenAK4JetEn_;
vector<float>  GenAK4JetPt_;
vector<float>  GenAK4JetEta_;
vector<float>  GenAK4JetDetEta_;
vector<float>  GenAK4JetPhi_;
vector<float>  GenAK4JetEMenergy_;
vector<float>  GenAK4JetHADenergy_;
vector<float>  GenAK4JetINVenergy_;


void ggNtuplizer::branchesGenAK4JetPart(TTree* tree) {

  tree->Branch("nGenAK4Jets",     &nGenAK4Jets_);
  tree->Branch("GenAK4JetEn",       &GenAK4JetEn_);
  tree->Branch("GenAK4JetPt",       &GenAK4JetPt_);
  tree->Branch("GenAK4JetEta",      &GenAK4JetEta_);
  tree->Branch("GenAK4JetDetEta",      &GenAK4JetDetEta_);
  tree->Branch("GenAK4JetPhi",      &GenAK4JetPhi_);
  tree->Branch("GenAK4JetEMenergy",      &GenAK4JetEMenergy_);
  tree->Branch("GenAK4JetHADenergy",      &GenAK4JetHADenergy_);
}

void ggNtuplizer::fillGenAK4JetInfo(const edm::Event& e, Float_t zPV) {

  nGenAK4Jets_ = -99;
  GenAK4JetEn_.clear();
  GenAK4JetPt_.clear();
  GenAK4JetEta_.clear();
  GenAK4JetDetEta_.clear();
  GenAK4JetPhi_.clear();
  GenAK4JetEMenergy_.clear();
  GenAK4JetHADenergy_.clear();
  GenAK4JetINVenergy_.clear();

    // Get GenJets
  edm::Handle<std::vector<reco::GenJet>> ak4PFJetsCHSgenjets;
  e.getByToken(ak4PFJetsCHSGenJetLabel_,ak4PFJetsCHSgenjets);
  if(!ak4PFJetsCHSgenjets.isValid()) {
    edm::LogWarning("ggNtuplizer") << "Could not find Ak4GenJet vector";
    return;
  }

  nGenAK4Jets_= ak4PFJetsCHSgenjets->size();

  for (vector<reco::GenJet>::const_iterator ip = ak4PFJetsCHSgenjets->begin(); ip != ak4PFJetsCHSgenjets->end(); ++ip) {
    GenAK4JetEn_.push_back(ip->energy());
    GenAK4JetPt_.push_back(ip->pt());
    GenAK4JetEta_.push_back(ip->eta());
    GenAK4JetDetEta_.push_back(ip->detectorEta(zPV));
    GenAK4JetPhi_.push_back(ip->phi());
    GenAK4JetEMenergy_.push_back(ip->emEnergy () );
    GenAK4JetHADenergy_.push_back(ip->hadEnergy ());
    GenAK4JetINVenergy_.push_back(ip->invisibleEnergy ());
  }
}


