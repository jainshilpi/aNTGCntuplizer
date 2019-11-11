#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
// #include "FWCore/Utilities/interface/RandomNumberGenerator.h"
// #include "CLHEP/Random/RandomEngine.h"
// #include "CLHEP/Random/RandGauss.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetResolutionObject.h"

using namespace std;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

UShort_t			nAK4CHSJet_;
vector<Char_t>  AK4CHSJet_Charge_;
vector<float>  AK4CHSJet_Pt_;
vector<float>  AK4CHSJet_En_;
vector<float>  AK4CHSJet_Eta_;
vector<float>  AK4CHSJet_Phi_;
vector<float>  AK4CHSJet_RawPt_;
vector<float>  AK4CHSJet_RawEn_;
vector<float>  AK4CHSJet_Mt_;
vector<float>  AK4CHSJet_Area_;
vector<float>  AK4CHSJet_LeadTrackPt_;
vector<float>  AK4CHSJet_LeadTrackEta_;
vector<float>  AK4CHSJet_LeadTrackPhi_;
vector<int>    AK4CHSJet_LepTrackPID_;
vector<float>  AK4CHSJet_LepTrackPt_;
vector<float>  AK4CHSJet_LepTrackEta_;
vector<float>  AK4CHSJet_LepTrackPhi_;
vector<float> AK4CHSJet_ConsituentEtaPhiSpread_;
vector<float> AK4CHSJet_MaxConstituentDistance_;
vector<float> AK4CHSJet_constituentPtDistribution_;
vector<float>  AK4CHSJet_CHF_;
vector<float>  AK4CHSJet_NHF_;
vector<float>  AK4CHSJet_CEF_;
vector<float>  AK4CHSJet_NEF_;
vector<UShort_t>    AK4CHSJet_NCH_;
vector<UShort_t>    AK4CHSJet_NNP_;
vector<float>  AK4CHSJet_MUF_;
vector<float>  AK4CHSJet_HFHAE_;
vector<float>  AK4CHSJet_HFEME_;
vector<UShort_t>    AK4CHSJet_NConstituents_;
vector<float>  AK4CHSJet_CSV2BJetTags_;
vector<float>  AK4CHSJet_DeepCSVTags_b_;
vector<float>  AK4CHSJet_DeepCSVTags_bb_;
vector<float>  AK4CHSJet_DeepCSVTags_c_;
vector<float>  AK4CHSJet_DeepCSVTags_udsg_;
vector<float>  AK4CHSJet_DeepFlavTags_bb_;
vector<float>  AK4CHSJet_DeepFlavTags_bbb_;
vector<float>  AK4CHSJet_DeepFlavTags_lepb_;
vector<float>  AK4CHSJet_DeepFlavTags_c_;
vector<float>  AK4CHSJet_DeepFlavTags_uds_;
vector<float>  AK4CHSJet_DeepFlavTags_g_;
vector<float>  AK4CHSJet_CombMVA2Tags_;
vector<int>    AK4CHSJet_PartonFlavour_;
vector<int>    AK4CHSJet_HadronFlavour_;
vector<Char_t>	AK4CHSJet_ID_;
vector<float>  AK4CHSJet_PUID_;
vector<Char_t>    AK4CHSJet_PUFullID_;
vector<Float_t>		AK4CHSJet_qgLikelihood_;
// vector<Float_t>    AK4CHSJet_JECfactor_;
vector<float>  AK4CHSJet_JECUnc_;
vector<Float_t> AK4CHSJet_JER_ptRes_;
vector<Float_t> AK4CHSJet_JER_phiRes_;
vector<Float_t> AK4CHSJet_JER_SF_;
vector<Float_t> AK4CHSJet_JER_SFup_;
vector<Float_t> AK4CHSJet_JER_SFdown_;
vector<ULong64_t> AK4CHSJet_FiredTrgs_;
//gen-info
vector<Short_t>  AK4CHSJet_GenJetIndex_;
vector<Short_t>    AK4CHSJet_GenPartonIndex_;


void ggNtuplizer::branchesAK4CHSJets(TTree* tree) {
	tree->Branch("nAK4CHSJet",                &nAK4CHSJet_);
	tree->Branch("AK4CHSJet_Charge",               &AK4CHSJet_Charge_);
	tree->Branch("AK4CHSJet_Pt",               &AK4CHSJet_Pt_);
	tree->Branch("AK4CHSJet_En",               &AK4CHSJet_En_);
	tree->Branch("AK4CHSJet_Eta",              &AK4CHSJet_Eta_);
	tree->Branch("AK4CHSJet_Phi",              &AK4CHSJet_Phi_);
	tree->Branch("AK4CHSJet_RawPt",            &AK4CHSJet_RawPt_);
	tree->Branch("AK4CHSJet_RawEn",            &AK4CHSJet_RawEn_);
	tree->Branch("AK4CHSJet_Mt",               &AK4CHSJet_Mt_);
	tree->Branch("AK4CHSJet_Area",             &AK4CHSJet_Area_);
	tree->Branch("AK4CHSJet_LeadTrackPt",      &AK4CHSJet_LeadTrackPt_);
	tree->Branch("AK4CHSJet_LeadTrackEta",     &AK4CHSJet_LeadTrackEta_);
	tree->Branch("AK4CHSJet_LeadTrackPhi",     &AK4CHSJet_LeadTrackPhi_);
	tree->Branch("AK4CHSJet_LepTrackPID",      &AK4CHSJet_LepTrackPID_);
	tree->Branch("AK4CHSJet_LepTrackPt",       &AK4CHSJet_LepTrackPt_);
	tree->Branch("AK4CHSJet_LepTrackEta",      &AK4CHSJet_LepTrackEta_);
	tree->Branch("AK4CHSJet_LepTrackPhi",      &AK4CHSJet_LepTrackPhi_);
	tree->Branch("AK4CHSJet_ConsituentEtaPhiSpread",      &AK4CHSJet_ConsituentEtaPhiSpread_);
	tree->Branch("AK4CHSJet_constituentPtDistribution",      &AK4CHSJet_constituentPtDistribution_);
	tree->Branch("AK4CHSJet_MaxConstituentDistance",      &AK4CHSJet_MaxConstituentDistance_);
	tree->Branch("AK4CHSJet_CSV2BJetTags",     &AK4CHSJet_CSV2BJetTags_);
	tree->Branch("AK4CHSJet_DeepCSVTags_b",    &AK4CHSJet_DeepCSVTags_b_);
	tree->Branch("AK4CHSJet_DeepCSVTags_bb",   &AK4CHSJet_DeepCSVTags_bb_);
	tree->Branch("AK4CHSJet_DeepCSVTags_c",    &AK4CHSJet_DeepCSVTags_c_);
	tree->Branch("AK4CHSJet_DeepCSVTags_udsg", &AK4CHSJet_DeepCSVTags_udsg_);
	tree->Branch("AK4CHSJet_DeepFlavTags_bb", &AK4CHSJet_DeepFlavTags_bb_);
	tree->Branch("AK4CHSJet_DeepFlavTags_bbb", &AK4CHSJet_DeepFlavTags_bbb_);
	tree->Branch("AK4CHSJet_DeepFlavTags_lepb", &AK4CHSJet_DeepFlavTags_lepb_);
	tree->Branch("AK4CHSJet_DeepFlavTags_c", &AK4CHSJet_DeepFlavTags_c_);
	tree->Branch("AK4CHSJet_DeepFlavTags_uds", &AK4CHSJet_DeepFlavTags_uds_);
	tree->Branch("AK4CHSJet_DeepFlavTags_g", &AK4CHSJet_DeepFlavTags_g_);
	tree->Branch("AK4CHSJet_CombMVA2Tags", &AK4CHSJet_CombMVA2Tags_);


	if (doGenParticles_){
		tree->Branch("AK4CHSJet_PartonFlavour",       &AK4CHSJet_PartonFlavour_);
		tree->Branch("AK4CHSJet_HadronFlavour",        &AK4CHSJet_HadronFlavour_);
		tree->Branch("AK4CHSJet_GenJetIndex",        &AK4CHSJet_GenJetIndex_);
		tree->Branch("AK4CHSJet_GenPartonIndex",    &AK4CHSJet_GenPartonIndex_);
	}
	tree->Branch("AK4CHSJet_ID",        &AK4CHSJet_ID_);
	tree->Branch("AK4CHSJet_PUID",      &AK4CHSJet_PUID_);
	tree->Branch("AK4CHSJet_PUFullID",  &AK4CHSJet_PUFullID_);
	tree->Branch("AK4CHSJet_qgLikelihood",  &AK4CHSJet_qgLikelihood_);

	// tree->Branch("AK4CHSJet_JECfactor",  &AK4CHSJet_JECfactor_);
	tree->Branch("AK4CHSJet_JECUnc",    &AK4CHSJet_JECUnc_);
	tree->Branch("AK4CHSJet_JER_ptRes",    &AK4CHSJet_JER_ptRes_);
	tree->Branch("AK4CHSJet_JER_phiRes",    &AK4CHSJet_JER_phiRes_);
	tree->Branch("AK4CHSJet_JER_SF",    &AK4CHSJet_JER_SF_);
	tree->Branch("AK4CHSJet_JER_SFup",    &AK4CHSJet_JER_SFup_);
	tree->Branch("AK4CHSJet_JER_SFdown",    &AK4CHSJet_JER_SFdown_);
	tree->Branch("AK4CHSJet_FiredTrgs", &AK4CHSJet_FiredTrgs_);
	tree->Branch("AK4CHSJet_CHF",       &AK4CHSJet_CHF_);
	tree->Branch("AK4CHSJet_NHF",       &AK4CHSJet_NHF_);
	tree->Branch("AK4CHSJet_CEF",       &AK4CHSJet_CEF_);
	tree->Branch("AK4CHSJet_NEF",       &AK4CHSJet_NEF_);
	tree->Branch("AK4CHSJet_NCH",       &AK4CHSJet_NCH_);
	tree->Branch("AK4CHSJet_NNP",       &AK4CHSJet_NNP_);
	tree->Branch("AK4CHSJet_MUF",       &AK4CHSJet_MUF_);
	if (development_) {
		tree->Branch("AK4CHSJet_HFHAE",         &AK4CHSJet_HFHAE_);
		tree->Branch("AK4CHSJet_HFEME",         &AK4CHSJet_HFEME_);
		tree->Branch("AK4CHSJet_NConstituents", &AK4CHSJet_NConstituents_);
	}
}


void ggNtuplizer::fillAK4CHSJets(const edm::Event& e, const edm::EventSetup& es) {
	AK4CHSJet_Charge_. clear();
	AK4CHSJet_Pt_                                  .clear();
	AK4CHSJet_En_                                  .clear();
	AK4CHSJet_Eta_                                 .clear();
	AK4CHSJet_Phi_                                 .clear();
	AK4CHSJet_RawPt_                               .clear();
	AK4CHSJet_RawEn_                               .clear();
	AK4CHSJet_Mt_                                  .clear();
	AK4CHSJet_Area_                                .clear();
	AK4CHSJet_LeadTrackPt_                         .clear();
	AK4CHSJet_LeadTrackEta_                        .clear();
	AK4CHSJet_LeadTrackPhi_                        .clear();
	AK4CHSJet_LepTrackPt_                          .clear();
	AK4CHSJet_LepTrackPID_                         .clear();
	AK4CHSJet_LepTrackEta_                         .clear();
	AK4CHSJet_LepTrackPhi_                         .clear();
	AK4CHSJet_ConsituentEtaPhiSpread_				.clear();
	AK4CHSJet_MaxConstituentDistance_				.clear();
	AK4CHSJet_constituentPtDistribution_			.clear();
	AK4CHSJet_CSV2BJetTags_                        .clear();
	AK4CHSJet_DeepCSVTags_b_                       .clear();
	AK4CHSJet_DeepCSVTags_bb_                      .clear();
	AK4CHSJet_DeepCSVTags_c_                       .clear();
	AK4CHSJet_DeepCSVTags_udsg_                    .clear();
	AK4CHSJet_DeepFlavTags_bb_ .clear();
	AK4CHSJet_DeepFlavTags_bbb_ .clear();
	AK4CHSJet_DeepFlavTags_lepb_ .clear();
	AK4CHSJet_DeepFlavTags_c_ .clear();
	AK4CHSJet_DeepFlavTags_uds_ .clear();
	AK4CHSJet_DeepFlavTags_g_ .clear();
	AK4CHSJet_CombMVA2Tags_ .clear();

	AK4CHSJet_PartonFlavour_                            .clear();
	AK4CHSJet_HadronFlavour_                             .clear();
	AK4CHSJet_ID_                                  .clear();
	AK4CHSJet_PUID_                                .clear();
	AK4CHSJet_PUFullID_                            .clear();
	AK4CHSJet_qgLikelihood_                        .clear();
	// AK4CHSJet_JECfactor_							.clear();
	AK4CHSJet_JECUnc_                              .clear();
	AK4CHSJet_JER_ptRes_                              .clear();
	AK4CHSJet_JER_phiRes_                             .clear();
	AK4CHSJet_JER_SF_                             .clear();
	AK4CHSJet_JER_SFup_                             .clear();
	AK4CHSJet_JER_SFdown_                             .clear();
	AK4CHSJet_FiredTrgs_                           .clear();
	AK4CHSJet_CHF_                                 .clear();
	AK4CHSJet_NHF_                                 .clear();
	AK4CHSJet_CEF_                                 .clear();
	AK4CHSJet_NEF_                                 .clear();
	AK4CHSJet_NCH_                                 .clear();
	AK4CHSJet_NNP_                                 .clear();
	AK4CHSJet_MUF_                                 .clear();
	if (development_) {
		AK4CHSJet_HFHAE_                               .clear();
		AK4CHSJet_HFEME_                               .clear();
		AK4CHSJet_NConstituents_                       .clear();
	}
	AK4CHSJet_GenJetIndex_.clear();
	AK4CHSJet_GenPartonIndex_.clear();

	nAK4CHSJet_ = 0;


	edm::Handle<edm::View<pat::Jet> > ak4PFJetsCHSHandle;
	e.getByToken(ak4PFJetsCHSLabel_, ak4PFJetsCHSHandle);

	edm::Handle<std::vector<reco::GenJet>> genJetsHandle;
	if (doGenParticles_) e.getByToken(ak4PFJetsCHSGenJetLabel_,genJetsHandle);

	edm::Handle<double> rhoHandle;
	e.getByToken(rhoLabel_, rhoHandle);
	float rho_ = *(rhoHandle.product());

	if (!ak4PFJetsCHSHandle.isValid()) {
		edm::LogWarning("ggNtuplizer") << "no pat::Jets (AK4) in event";
		return;
	}

	edm::Handle<vector<reco::GenParticle> > genParticlesHandle;
	if(doGenParticles_)e.getByToken(genParticlesCollection_, genParticlesHandle);

	edm::Handle<reco::VertexCollection> vtxHandle;
	e.getByToken(vtxLabel_, vtxHandle);
	if (!vtxHandle.isValid()) edm::LogWarning("ggNtuplizer") << "Primary vertices info not unavailable";

  	// Accessing the JEC uncertainties
	edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
	es.get<JetCorrectionsRecord>().get("AK4PFchs",JetCorParColl);
	JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
	JetCorrectionUncertainty *jecUnc=0;
	jecUnc = new JetCorrectionUncertainty(JetCorPar);

	// Access JER
	JME::JetResolution jet_pt_resolution_ = JME::JetResolution::get(es, "AK4PFchs_pt");
	JME::JetResolution jet_phi_resolution_ = JME::JetResolution::get(es, "AK4PFchs_phi");
	JME::JetResolutionScaleFactor resolution_sf_ = JME::JetResolutionScaleFactor::get(es, "AK4PFchs");
	JME::JetParameters jet_parameters;

	for (edm::View<pat::Jet>::const_iterator iJet = ak4PFJetsCHSHandle->begin(); iJet != ak4PFJetsCHSHandle->end(); ++iJet) {

		// if (iJet->pt() < 30.) continue;

		// if((!iJet->isPFJet()) && (!iJet->isJPTJet())) continue;

		AK4CHSJet_Charge_.push_back(    iJet->charge());
		AK4CHSJet_Pt_.push_back(    iJet->pt());
		AK4CHSJet_En_.push_back(    iJet->energy());
		AK4CHSJet_Eta_.push_back(   iJet->eta());
		AK4CHSJet_Phi_.push_back(   iJet->phi());
		AK4CHSJet_RawPt_.push_back( iJet->correctedJet("Uncorrected").pt());
		AK4CHSJet_RawEn_.push_back( iJet->correctedJet("Uncorrected").energy());
		AK4CHSJet_Mt_.push_back(    iJet->mt());
		AK4CHSJet_Area_.push_back(  iJet->jetArea());
		AK4CHSJet_CEF_.push_back(   iJet->chargedEmEnergyFraction());
		AK4CHSJet_NEF_.push_back(   iJet->neutralEmEnergyFraction());
		AK4CHSJet_CHF_.push_back(   iJet->chargedHadronEnergyFraction());
		AK4CHSJet_NHF_.push_back(   iJet->neutralHadronEnergyFraction());
		AK4CHSJet_NNP_.push_back(   iJet->neutralMultiplicity());
		AK4CHSJet_MUF_.push_back(   iJet->muonEnergyFraction());
		if (development_) {
			AK4CHSJet_HFHAE_.push_back( iJet->HFHadronEnergy());
			AK4CHSJet_HFEME_.push_back( iJet->HFEMEnergy());
			AK4CHSJet_NConstituents_.push_back(iJet->numberOfDaughters());
		}

		// JEC https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections
		if (std::abs(iJet->eta()) < 5.2){
			jecUnc->setJetEta(iJet->eta());
      		jecUnc->setJetPt(iJet->pt()); // here you must use the CORRECTED jet pt
      		AK4CHSJet_JECUnc_.push_back(jecUnc->getUncertainty(true));
      	} else {
      		AK4CHSJet_JECUnc_.push_back(-1.);
      	}
      	// AK4CHSJet_JECfactor_.push_back(iJet->jecFactor());

      	AK4CHSJet_FiredTrgs_.push_back(matchJetTriggerFilters(iJet->pt(), iJet->eta(), iJet->phi()));


      	if(doGenParticles_){
    		//Searching for leading track and lepton
      		float leadTrkPt  = -99;
      		float leadTrkEta = -99;
      		float leadTrkPhi = -99;
      		int   lepTrkPID  = -99;
      		float lepTrkPt   = -99;
      		float lepTrkEta  = -99;
      		float lepTrkPhi  = -99;

      		for (unsigned id = 0; id < iJet->getJetConstituents().size(); id++) {

      			const edm::Ptr<reco::Candidate> daughter = iJet->getJetConstituents().at(id);

      			if (daughter.isNonnull() && daughter.isAvailable()) {
      				if (daughter->charge() != 0 && daughter->pt() > leadTrkPt) {
      					leadTrkPt  = daughter->pt();
      					leadTrkEta = daughter->eta();
      					leadTrkPhi = daughter->phi();
      				}

      				if (abs(daughter->pdgId()) == 11 || abs(daughter->pdgId()) == 13) {
      					if (daughter->pt() > lepTrkPt) {
      						lepTrkPID = daughter->pdgId();
      						lepTrkPt  = daughter->pt();
      						lepTrkEta = daughter->eta();
      						lepTrkPhi = daughter->phi();
      					}
      				}
      			}
      		}

      		AK4CHSJet_LeadTrackPt_ .push_back(leadTrkPt);
      		AK4CHSJet_LeadTrackEta_.push_back(leadTrkEta);
      		AK4CHSJet_LeadTrackPhi_.push_back(leadTrkPhi);
      		AK4CHSJet_LepTrackPID_ .push_back(lepTrkPID);
      		AK4CHSJet_LepTrackPt_  .push_back(lepTrkPt);
      		AK4CHSJet_LepTrackEta_ .push_back(lepTrkEta);
      		AK4CHSJet_LepTrackPhi_ .push_back(lepTrkPhi);
      	}

      	AK4CHSJet_ConsituentEtaPhiSpread_.push_back(iJet->constituentEtaPhiSpread());
      	AK4CHSJet_MaxConstituentDistance_.push_back(iJet->maxDistance());
      	AK4CHSJet_constituentPtDistribution_.push_back(iJet->constituentPtDistribution());

    	//b/c-tagging
      	AK4CHSJet_CSV2BJetTags_    .push_back(iJet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
      	AK4CHSJet_DeepCSVTags_b_   .push_back(iJet->bDiscriminator("pfDeepCSVJetTags:probb"));
      	AK4CHSJet_DeepCSVTags_bb_  .push_back(iJet->bDiscriminator("pfDeepCSVJetTags:probbb"));
      	AK4CHSJet_DeepCSVTags_c_   .push_back(iJet->bDiscriminator("pfDeepCSVJetTags:probc"));
      	AK4CHSJet_DeepCSVTags_udsg_.push_back(iJet->bDiscriminator("pfDeepCSVJetTags:probudsg"));
      	AK4CHSJet_DeepFlavTags_bb_ .push_back(iJet->bDiscriminator("pfDeepFlavourJetTags:probb"));
      	AK4CHSJet_DeepFlavTags_bbb_ .push_back(iJet->bDiscriminator("pfDeepFlavourJetTags:probbb"));
      	AK4CHSJet_DeepFlavTags_lepb_ .push_back(iJet->bDiscriminator("pfDeepFlavourJetTags:problepb"));
      	AK4CHSJet_DeepFlavTags_c_ .push_back(iJet->bDiscriminator("pfDeepFlavourJetTags:probc"));
      	AK4CHSJet_DeepFlavTags_uds_ .push_back(iJet->bDiscriminator("pfDeepFlavourJetTags:probuds"));
      	AK4CHSJet_DeepFlavTags_g_ .push_back(iJet->bDiscriminator("pfDeepFlavourJetTags:probg"));
      	AK4CHSJet_CombMVA2Tags_ .push_back(iJet->bDiscriminator("pfCombinedMVAV2BJetTags"));

    	//parton id
      	AK4CHSJet_PartonFlavour_.push_back(iJet->partonFlavour());
      	AK4CHSJet_HadronFlavour_.push_back(iJet->hadronFlavour());

    	//jet PF Loose ID

      	double NHF      = iJet->neutralHadronEnergyFraction();
      	double NEMF     = iJet->neutralEmEnergyFraction();
      	double CHF      = iJet->chargedHadronEnergyFraction();
      	Int_t CHM      = iJet->chargedMultiplicity();
      	Int_t NNP      = iJet->neutralMultiplicity();
      	Int_t NumConst = CHM + NNP;
      	double CEMF     = iJet->chargedEmEnergyFraction();
      	double MUF      = iJet->muonEnergyFraction();

      	AK4CHSJet_NCH_.push_back( CHM);

      	//https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2017?rev=7
      	bool looseJetID = false;
      	bool tightJetID = false;
      	Bool_t tightLeptVetoID = false;
      	if (fabs(iJet->eta()) <= 2.7) {
      		looseJetID = (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((fabs(iJet->eta())<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || fabs(iJet->eta())>2.4);
      		tightJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((fabs(iJet->eta())<=2.4 && CHF>0 && CHM>0) || fabs(iJet->eta())>2.4);
      		tightLeptVetoID = (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF < 0.8) && ((fabs(iJet->eta())<=2.4 && CHF>0 && CHM>0 && CEMF<0.8) || fabs(iJet->eta())>2.4);
      	} else if (fabs(iJet->eta()) <= 3.0) {
      		looseJetID = (NEMF>0.01 && NHF<0.98 && NNP>2);
      		tightJetID = (NEMF>0.02 && NEMF <0.99 && NNP>2);
      		tightLeptVetoID = tightJetID;
      	} else {
      		looseJetID = (NEMF<0.90 && NNP>10);
      		tightJetID = (NEMF<0.90 && NNP>10 && NHF > 0.02);
      		tightLeptVetoID = tightJetID;
      	}
      	Char_t jetIDdecision = 0;
      	if(looseJetID) setbit(jetIDdecision, 0);
      	if(tightJetID) setbit(jetIDdecision, 1);
      	if(tightLeptVetoID) setbit(jetIDdecision, 2);
      	AK4CHSJet_ID_.push_back(jetIDdecision);

    	// PUJet ID from slimmedJets - not available for PUPPI

      	Float_t tmpPUID = -999;
      	Char_t tmpPUfullID = -99;

      	if(iJet->hasUserFloat("AK4PFCHSupdatedpileupJetIdEvaluator.:fullDiscriminant")) tmpPUID = iJet->userFloat("AK4PFCHSupdatedpileupJetIdEvaluator:fullDiscriminant");
      	else if(iJet->hasUserFloat("pileupJetId:fullDiscriminant")) tmpPUID = iJet->userFloat("pileupJetId:fullDiscriminant");

      	if(iJet->hasUserInt("AK4PFCHSupdatedpileupJetIdEvaluator:fullId")) tmpPUfullID = iJet->userInt("AK4PFCHSupdatedpileupJetIdEvaluator:fullId");
      	else if(iJet->hasUserInt("pileupJetId:fullId")) tmpPUfullID = iJet->userInt("pileupJetId:fullId");

      	AK4CHSJet_PUID_.push_back(tmpPUID); // not for PUPPI https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2017#Jets
      	AK4CHSJet_PUFullID_.emplace_back(tmpPUfullID);

      	Float_t tmp_qgLikelihood = -999.;
      	if(iJet->hasUserFloat("QGTagger:qgLikelihood")) tmp_qgLikelihood = iJet->userFloat("QGTagger:qgLikelihood");
      	else if(iJet->hasUserFloat("QGTaggerAK4PFCHS:qgLikelihood")) tmp_qgLikelihood = iJet->userFloat("QGTaggerAK4PFCHS:qgLikelihood");
      	else if(iJet->hasUserFloat("QGTaggerAK4PFCHSupdated:qgLikelihood")) tmp_qgLikelihood = iJet->userFloat("QGTaggerAK4PFCHSupdated:qgLikelihood");
      	AK4CHSJet_qgLikelihood_.emplace_back(tmp_qgLikelihood);
    	// gen jet and parton
      	if (doGenParticles_) {
      		int jetGenPartonIndex = -99;
      		if (iJet->genParton() && genParticlesHandle.isValid()) {
      			jetGenPartonIndex = std::distance(genParticlesHandle->begin(), (vector<reco::GenParticle>::const_iterator) iJet->genParton());
      		}
      		AK4CHSJet_GenPartonIndex_.push_back(jetGenPartonIndex);
      		Int_t _genJetIndex = -99;
      		if (iJet->genJet() && genJetsHandle.isValid()) {
      			_genJetIndex = std::distance(genJetsHandle->begin(), (vector<reco::GenJet>::const_iterator) iJet->genJet());
      			// cout<<"  "<<iJet->genJet()->pt()<<"  "<<(genJetsHandle->begin()+_genJetIndex)->pt()<<endl;
      		}
      		AK4CHSJet_GenJetIndex_.push_back(_genJetIndex);
      	}

      	// JER https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyResolution
      	jet_parameters.setJetPt(iJet->pt()).setJetEta(iJet->eta()).setRho(rho_);

      	Float_t _jer_ptRes = jet_pt_resolution_.getResolution(jet_parameters);
      	Float_t _jer_phiRes = jet_phi_resolution_.getResolution(jet_parameters);
      	Float_t _jer_sf = resolution_sf_.getScaleFactor(jet_parameters);
      	Float_t _jersf_up = resolution_sf_.getScaleFactor(jet_parameters, Variation::UP);
      	Float_t _jersf_down = resolution_sf_.getScaleFactor(jet_parameters, Variation::DOWN);

      	AK4CHSJet_JER_ptRes_.push_back(_jer_ptRes);
      	AK4CHSJet_JER_phiRes_.push_back(_jer_phiRes);
      	AK4CHSJet_JER_SF_.push_back(_jer_sf);
      	AK4CHSJet_JER_SFup_.push_back(_jersf_up);
      	AK4CHSJet_JER_SFdown_.push_back(_jersf_down);

      	nAK4CHSJet_++;
      }

      delete jecUnc;
  }