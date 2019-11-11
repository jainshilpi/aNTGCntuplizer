#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"
#include "ggAnalysis/ggNtuplizer/interface/GenParticleParentage.h"
using namespace std;

// (local) variables associated with tree branches
UShort_t        nGenAK8Jets_;
vector<float>  GenAK8JetEn_;
vector<float>  GenAK8JetPt_;
vector<float>  GenAK8JetEta_;
vector<float>  GenAK8JetDetEta_;
vector<float>  GenAK8JetPhi_;
vector<float>  GenAK8JetEMenergy_;
vector<float>  GenAK8JetHADenergy_;
vector<float>  GenAK8JetINVenergy_;


void ggNtuplizer::branchesGenAK8JetPart(TTree* tree) {

  tree->Branch("nGenAK8Jets",     &nGenAK8Jets_);
  tree->Branch("GenAK8JetEn",       &GenAK8JetEn_);
  tree->Branch("GenAK8JetPt",       &GenAK8JetPt_);
  tree->Branch("GenAK8JetEta",      &GenAK8JetEta_);
  tree->Branch("GenAK8JetDetEta",      &GenAK8JetDetEta_);
  tree->Branch("GenAK8JetPhi",      &GenAK8JetPhi_);
  tree->Branch("GenAK8JetEMenergy",      &GenAK8JetEMenergy_);
  tree->Branch("GenAK8JetHADenergy",      &GenAK8JetHADenergy_);
}

void ggNtuplizer::fillGenAK8JetInfo(const edm::Event& e, Float_t zPV) {

  nGenAK8Jets_ = -99;
  GenAK8JetEn_.clear();
  GenAK8JetPt_.clear();
  GenAK8JetEta_.clear();
  GenAK8JetDetEta_.clear();
  GenAK8JetPhi_.clear();
  GenAK8JetEMenergy_.clear();
  GenAK8JetHADenergy_.clear();
  GenAK8JetINVenergy_.clear();

    // Get GenJets
  edm::Handle<std::vector<reco::GenJet>> ak8genjetsHandle;
  e.getByToken(ak8GenJetLabel_,ak8genjetsHandle);
  if(!ak8genjetsHandle.isValid()) {
    edm::LogWarning("ggNtuplizer") << "Could not find Ak8GenJet vector";
    return;
  }

  nGenAK8Jets_= ak8genjetsHandle->size();

  for (vector<reco::GenJet>::const_iterator ip = ak8genjetsHandle->begin(); ip != ak8genjetsHandle->end(); ++ip) {
    GenAK8JetEn_.push_back(ip->energy());
    GenAK8JetPt_.push_back(ip->pt());
    GenAK8JetEta_.push_back(ip->eta());
    GenAK8JetDetEta_.push_back(ip->detectorEta(zPV));
    GenAK8JetPhi_.push_back(ip->phi());
    GenAK8JetEMenergy_.push_back(ip->emEnergy());
    GenAK8JetHADenergy_.push_back(ip->hadEnergy());
    GenAK8JetINVenergy_.push_back(ip->invisibleEnergy());
  }
}


