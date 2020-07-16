#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"
#include "DataFormats/EgammaReco/interface/HFEMClusterShape.h"
#include "RecoEgamma/EgammaHFProducers/plugins/HFEGammaSLCorrector.cc"

using namespace std;
using namespace reco;

// (local) variables associated with tree branches
UShort_t          nhfEle_;
vector<float>  hfeleEn_;
vector<float>  hfeleSCEn_;
vector<float>  hfelePt_;
vector<float>  hfeleEta_;
vector<float>  hfelePhi_;
vector<float>  hfeleSCEta_;
vector<float>  hfeleSCPhi_;

vector<float>  hfeleELong3x3_;
vector<float>  hfeleELong5x5_;
vector<float>  hfeleELong1x1_;

vector<float>  hfeleEShort3x3_;
vector<float>  hfeleEShort5x5_;
vector<float>  hfeleEShort1x1_;

vector<float>  hfeleESeL_;
vector<float>  hfeleECOREe9_;
vector<float>  hfeleE9E25_;

vector<float>  hfeleECore_;
vector<float>  hfeleE9E25L_;
vector<float>  hfeleE1E9_;
vector<float>  hfeleVar2D_;


void ggNtuplizer::branchesHFElectrons(TTree* tree) {

  tree->Branch("nhfEle",                    &nhfEle_);
  // tree->Branch("eleChargeConsistent",     &eleChargeConsistent_);
  tree->Branch("hfeleEn",                   &hfeleEn_);
  tree->Branch("hfeleSCEn",                 &hfeleSCEn_);
  tree->Branch("hfelePt",                   &hfelePt_);
  tree->Branch("hfeleEta",                  &hfeleEta_);
  tree->Branch("hfelePhi",                  &hfelePhi_);
  tree->Branch("hfeleSCEta",                &hfeleSCEta_);
  tree->Branch("hfeleSCPhi",                &hfeleSCPhi_);

  tree->Branch("hfeleELong3x3",                &hfeleELong3x3_);
  tree->Branch("hfeleELong5x5",                &hfeleELong5x5_);
  tree->Branch("hfeleELong1x1",                &hfeleELong1x1_);

  tree->Branch("hfeleEShort3x3",                &hfeleEShort3x3_);
  tree->Branch("hfeleEShort5x5",                &hfeleEShort5x5_);
  tree->Branch("hfeleEShort1x1",                &hfeleEShort1x1_);

  tree->Branch("hfeleESeL",                &hfeleESeL_);
  tree->Branch("hfeleECOREe9",                &hfeleECOREe9_);
  tree->Branch("hfeleE9E25",                &hfeleE9E25_);

  tree->Branch("hfeleECore",                &hfeleECore_);
  tree->Branch("hfeleE9E25L",                &hfeleE9E25L_);
  tree->Branch("hfeleE1E9",                &hfeleE1E9_);

  tree->Branch("hfeleVar2D",                &hfeleVar2D_);


}

void ggNtuplizer::fillHFElectrons(const edm::Event &e, const edm::EventSetup &es, math::XYZPoint &pv) {

  hfeleEn_                      .clear();
  hfeleSCEn_                    .clear();
  hfelePt_                      .clear();
  hfeleEta_                     .clear();
  hfelePhi_                     .clear();
  hfeleSCEta_                   .clear();
  hfeleSCPhi_                   .clear();

  nhfEle_ = 0;

  //vector<reco::RecoEcalCandidate>       "hfRecoEcalCandidate"       ""                "RECO"
  edm::Handle<reco::RecoEcalCandidateCollection> electronHandle;
  e.getByToken(hfelectronCollection_, electronHandle);

  edm::Handle<reco::HFEMClusterShapeAssociationCollection> hf_assoc;
  e.getByToken(hfclustersHFEM_, hf_assoc);

  if (!electronHandle.isValid()) {
    edm::LogWarning("ggNtuplizer") << "no pat::Electrons in event";
    return;
  }

  for (reco::RecoEcalCandidateCollection::const_iterator iEle = electronHandle->begin(); iEle != electronHandle->end(); ++iEle) {

    hfeleEn_              .push_back(iEle->energy());
    hfelePt_              .push_back(iEle->pt());
    hfeleEta_             .push_back(iEle->eta());
    hfelePhi_             .push_back(iEle->phi());
    hfeleSCEn_            .push_back(iEle->superCluster()->energy());
    hfeleSCEta_             .push_back(iEle->superCluster()->eta());
    hfeleSCPhi_             .push_back(iEle->superCluster()->phi());
    
    const HFEMClusterShapeRef clusShapeRef = hf_assoc->find(iEle->superCluster())->val;
    const HFEMClusterShape& clusShape = *clusShapeRef;

    double m_intercept2DSlope = 0.475 ;
    double e9e25 = clusShape.eLong3x3() / clusShape.eLong5x5();
    double e1e9 = clusShape.eLong1x1() / clusShape.eLong3x3();
    double eSeL = hf_egamma::eSeLCorrected(clusShape.eShort3x3(), clusShape.eLong3x3(), 4);
    double var2d = (clusShape.eCOREe9() - (eSeL * m_intercept2DSlope));

    hfeleELong3x3_    .push_back(clusShape.eLong3x3());
    hfeleELong5x5_    .push_back(clusShape.eLong5x5());
    hfeleELong1x1_    .push_back(clusShape.eLong1x1());


    hfeleEShort3x3_    .push_back(clusShape.eShort3x3());
    hfeleEShort5x5_    .push_back(clusShape.eShort5x5());
    hfeleEShort1x1_    .push_back(clusShape.eShort1x1());

    hfeleESeL_         .push_back(clusShape.eSeL());
    hfeleECOREe9_         .push_back(clusShape.eCOREe9());
    hfeleE9E25_         .push_back(clusShape.e9e25());
    hfeleECore_         .push_back(clusShape.eCore());

    hfeleE9E25L_        .push_back(e9e25); 
    hfeleE1E9_        .push_back(e1e9); 
    
    hfeleVar2D_       .push_back(var2d);

    nhfEle_++;
  }
};
