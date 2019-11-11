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
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfo.h"

////////////////////////////////////////////////////////////////////////
/////Generate ROOT dictionary for std::vector<std::vector<Char_t> >
// #include <vector>
// #ifdef __MAKECINT__
// #pragma link C++ class std::vector<std::vector<Char_t> >+;
// #endif
////////////////////////////////////////////////////////////////////////

using namespace std;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

UShort_t         nAK8PuppiJet_;
vector<Char_t> AK8PuppiJet_Charge_;
vector<float> AK8PuppiJet_Pt_;
vector<float> AK8PuppiJet_En_;
vector<float> AK8PuppiJet_RawPt_;
vector<float> AK8PuppiJet_RawEn_;
vector<float> AK8PuppiJet_Eta_;
vector<float> AK8PuppiJet_Phi_;
vector<float> AK8PuppiJet_Mass_;
vector<float> AK8PuppiJet_tau1_;
vector<float> AK8PuppiJet_tau2_;
vector<float> AK8PuppiJet_tau3_;
vector<float> AK8PuppiJet_tau4_;
vector<float> AK8PuppiJet_CHF_;
vector<float> AK8PuppiJet_NHF_;
vector<float> AK8PuppiJet_CEF_;
vector<float> AK8PuppiJet_NEF_;
vector<UChar_t>   AK8PuppiJet_NCH_;
vector<UChar_t>   AK8PuppiJet_NNP_;
vector<float> AK8PuppiJet_MUF_;
// vector<float> AK8PuppiJet_puppiMultiplicity_;
// vector<float> AK8PuppiJet_neutral_puppiMultiplicity_;
vector<Char_t> AK8PuppiJet_PFid_;
vector<float> AK8PuppiJet_PUPPISoftDropMass_;
vector<float> AK8PuppiJet_CHSPrunedMass_;
vector<float> AK8PuppiJet_CHSSoftDropMass_;
vector<float> AK8PuppiJet_pfBoostedDSVBTag_;

// vector<float> AK8PuppiJet_TrimmedMass_;
// vector<float> AK8PuppiJet_FilteredMass_;
// vector<float> AK8PuppiJet_DSVnewV4_;
vector<float> AK8PuppiJet_CSV_;
vector<float> AK8PuppiJet_DDBvLTags_qcd_;
vector<float> AK8PuppiJet_DDBvLTags_Hbb_;
vector<float> AK8PuppiJet_DDCvLTags_qcd_;
vector<float> AK8PuppiJet_DDCvLTags_Hcc_;
vector<float> AK8PuppiJet_DDCvBTags_Hbb_;
vector<float> AK8PuppiJet_DDCvBTags_Hcc_;
vector<float>  AK8PuppiJet_DDBDTags_bbcsLight_;
vector<float>  AK8PuppiJet_DDBDTags_ccvsLight_;
vector<float>  AK8PuppiJet_DDBDTags_TvsQCD_;
vector<float>  AK8PuppiJet_DDBDTags_ZHccvsQCD_;
vector<float>  AK8PuppiJet_DDBDTags_WvsQCD_;
vector<float>  AK8PuppiJet_DDBDTags_ZHbbvsQCD_;
vector<float> AK8PuppiJet_JECUnc_;
// vector<float> AK8PuppiJet_L2L3corr_;
vector<float> AK8PuppiJet_JER_sf_;
vector<float> AK8PuppiJet_JER_sf_up_;
vector<float> AK8PuppiJet_JER_sf_down_;
vector<Float_t> AK8PuppiJet_JER_ptRes_;
vector<Float_t> AK8PuppiJet_JER_phiRes_;
vector<float>  AK8PuppiJet_nb1ecf2_;
vector<float>  AK8PuppiJet_nb1ecf3_;
vector<float>  AK8PuppiJet_nb2ecf2_;
vector<float>  AK8PuppiJet_nb2ecf3_;
vector<int>   AK8PuppiJet_PartonFlavour_;
vector<int>   AK8PuppiJet_HadronFlavour_;
// vector<float> AK8PuppiJet_GenJetEn_;
// vector<float> AK8PuppiJet_GenJetPt_;
// vector<float> AK8PuppiJet_GenJetEta_;
// vector<float> AK8PuppiJet_GenJetPhi_;
vector<Short_t>   AK8PuppiJet_GenJetIndex_;
// vector<vector<Short_t>> AK8PuppiJet_GenPartonPID_;
vector<UShort_t>             AK8PuppiJet_nSDSJpuppi_ ;
vector< vector<float> > AK8PuppiJet_SDSJpuppiPt_ ;
vector< vector<float> > AK8PuppiJet_SDSJpuppiEta_ ;
vector< vector<float> > AK8PuppiJet_SDSJpuppiMass_ ;
vector< vector<float> > AK8PuppiJet_SDSJpuppiPhi_ ;
vector< vector<float> > AK8PuppiJet_SDSJpuppiE_ ;
vector< vector<Char_t > >  AK8PuppiJet_SDSJpuppiCharge_ ;
vector< vector<int > >  AK8PuppiJet_SDSJpuppiFlavour_;
vector< vector<float> > AK8PuppiJet_SDSJpuppiCSV_ ;
vector< vector<float> > AK8PuppiJet_SDSJpuppinb1ecf2_;
vector< vector<float> > AK8PuppiJet_SDSJpuppinb1ecf3_;
vector< vector<float> > AK8PuppiJet_SDSJpuppinb2ecf2_;
vector< vector<float> > AK8PuppiJet_SDSJpuppinb2ecf3_;


void ggNtuplizer::branchesAK8PUPPIJets(TTree* tree) {
	tree->Branch("nAK8PUPPIJet",                  &nAK8PuppiJet_);
	tree->Branch("AK8PuppiJet_Charge",                 &AK8PuppiJet_Charge_);
	tree->Branch("AK8PuppiJet_Pt",                 &AK8PuppiJet_Pt_);
	tree->Branch("AK8PuppiJet_En",                 &AK8PuppiJet_En_);
	tree->Branch("AK8PuppiJet_RawPt",              &AK8PuppiJet_RawPt_);
	tree->Branch("AK8PuppiJet_RawEn",              &AK8PuppiJet_RawEn_);
	tree->Branch("AK8PuppiJet_Eta",                &AK8PuppiJet_Eta_);
	tree->Branch("AK8PuppiJet_Phi",                &AK8PuppiJet_Phi_);
	tree->Branch("AK8PuppiJet_Mass",               &AK8PuppiJet_Mass_);
	tree->Branch("AK8PuppiJet_tau1",              &AK8PuppiJet_tau1_);
	tree->Branch("AK8PuppiJet_tau2",              &AK8PuppiJet_tau2_);
	tree->Branch("AK8PuppiJet_tau3",              &AK8PuppiJet_tau3_);
	tree->Branch("AK8PuppiJet_tau4",              &AK8PuppiJet_tau3_);
	tree->Branch("AK8PuppiJet_CHF",                &AK8PuppiJet_CHF_);
	tree->Branch("AK8PuppiJet_NHF",                &AK8PuppiJet_NHF_);
	tree->Branch("AK8PuppiJet_CEF",                &AK8PuppiJet_CEF_);
	tree->Branch("AK8PuppiJet_NEF",                &AK8PuppiJet_NEF_);
	tree->Branch("AK8PuppiJet_NCH",                &AK8PuppiJet_NCH_);
	tree->Branch("AK8PuppiJet_NNP",                &AK8PuppiJet_NNP_);
	tree->Branch("AK8PuppiJet_MUF",                &AK8PuppiJet_MUF_);
	// tree->Branch("AK8PuppiJet_puppiMultiplicity",                &AK8PuppiJet_puppiMultiplicity_);
	// tree->Branch("AK8PuppiJet_neutral_puppiMultiplicity",                &AK8PuppiJet_neutral_puppiMultiplicity_);
	tree->Branch("AK8PuppiJet_PFid",          &AK8PuppiJet_PFid_);
	tree->Branch("AK8PuppiJet_PUPPISoftDropMass",       &AK8PuppiJet_PUPPISoftDropMass_);
	// tree->Branch("AK8PuppiJet_TrimmedMass",       &AK8PuppiJet_TrimmedMass_);
	// tree->Branch("AK8PuppiJet_FilteredMass",       &AK8PuppiJet_FilteredMass_);
	tree->Branch("AK8PuppiJet_CHSPrunedMass",         &AK8PuppiJet_CHSPrunedMass_);
	tree->Branch("AK8PuppiJet_CHSSoftDropMass",         &AK8PuppiJet_CHSSoftDropMass_);
	tree->Branch("AK8PuppiJet_pfBoostedDSVBTag",   &AK8PuppiJet_pfBoostedDSVBTag_);
	// tree->Branch("AK8PuppiJet_DSVnewV4",           &AK8PuppiJet_DSVnewV4_);
	tree->Branch("AK8PuppiJet_CSV",                &AK8PuppiJet_CSV_);
	tree->Branch("AK8PuppiJet_DDBvLTags_qcd",                &AK8PuppiJet_DDBvLTags_qcd_);
	tree->Branch("AK8PuppiJet_DDBvLTags_Hbb",                &AK8PuppiJet_DDBvLTags_Hbb_);
	tree->Branch("AK8PuppiJet_DDCvLTags_qcd",                &AK8PuppiJet_DDCvLTags_qcd_);
	tree->Branch("AK8PuppiJet_DDCvLTags_Hcc",                &AK8PuppiJet_DDCvLTags_Hcc_);
	tree->Branch("AK8PuppiJet_DDCvBTags_Hbb",                &AK8PuppiJet_DDCvBTags_Hbb_);
	tree->Branch("AK8PuppiJet_DDCvBTags_Hcc",                &AK8PuppiJet_DDCvBTags_Hcc_);
	tree->Branch("AK8PuppiJet_DDBDTags_bbcsLight",                &AK8PuppiJet_DDBDTags_bbcsLight_);
	tree->Branch("AK8PuppiJet_DDBDTags_ccvsLight",                &AK8PuppiJet_DDBDTags_ccvsLight_);
	tree->Branch("AK8PuppiJet_DDBDTags_TvsQCD",                &AK8PuppiJet_DDBDTags_TvsQCD_);
	tree->Branch("AK8PuppiJet_DDBDTags_ZHccvsQCD",                &AK8PuppiJet_DDBDTags_ZHccvsQCD_);
	tree->Branch("AK8PuppiJet_DDBDTags_WvsQCD",                &AK8PuppiJet_DDBDTags_WvsQCD_);
	tree->Branch("AK8PuppiJet_DDBDTags_ZHbbvsQCD",                &AK8PuppiJet_DDBDTags_ZHbbvsQCD_);
	tree->Branch("AK8PuppiJet_JECUnc",             &AK8PuppiJet_JECUnc_);
	// tree->Branch("AK8PuppiJet_L2L3corr",           &AK8PuppiJet_L2L3corr_);
	tree->Branch("AK8PuppiJet_JER_sf",           &AK8PuppiJet_JER_sf_);
	tree->Branch("AK8PuppiJet_JER_sf_up",           &AK8PuppiJet_JER_sf_up_);
	tree->Branch("AK8PuppiJet_JER_sf_down",           &AK8PuppiJet_JER_sf_down_);
	tree->Branch("AK8PuppiJet_JER_ptRes",           &AK8PuppiJet_JER_ptRes_);
	tree->Branch("AK8PuppiJet_JER_phiRes",           &AK8PuppiJet_JER_phiRes_);
	tree->Branch("AK8PuppiJet_nb1ecf2",           &AK8PuppiJet_nb1ecf2_);
	tree->Branch("AK8PuppiJet_nb1ecf3",           &AK8PuppiJet_nb1ecf3_);
	tree->Branch("AK8PuppiJet_nb2ecf2",           &AK8PuppiJet_nb2ecf2_);
	tree->Branch("AK8PuppiJet_nb2ecf3",           &AK8PuppiJet_nb2ecf3_);
	if (doGenParticles_) {
		tree->Branch("AK8PuppiJet_PartonFlavour",       &AK8PuppiJet_PartonFlavour_);
		tree->Branch("AK8PuppiJet_HadronFlavour",        &AK8PuppiJet_HadronFlavour_);
		tree->Branch("AK8PuppiJet_GenJetIndex",    &AK8PuppiJet_GenJetIndex_);
		// tree->Branch("AK8PuppiJet_GenPartonPID",    &AK8PuppiJet_GenPartonPID_);
		// tree->Branch("AK8PuppiJet_GenJetEn",       &AK8PuppiJet_GenJetEn_);
		// tree->Branch("AK8PuppiJet_GenJetPt",       &AK8PuppiJet_GenJetPt_);
		// tree->Branch("AK8PuppiJet_GenJetEta",      &AK8PuppiJet_GenJetEta_);
		// tree->Branch("AK8PuppiJet_GenJetPhi",      &AK8PuppiJet_GenJetPhi_);
	}
	tree->Branch("AK8PuppiJet_nSDSJpuppi",            &AK8PuppiJet_nSDSJpuppi_);
	tree->Branch("AK8PuppiJet_SDSJpuppiPt",           &AK8PuppiJet_SDSJpuppiPt_);
	tree->Branch("AK8PuppiJet_SDSJpuppiEta",          &AK8PuppiJet_SDSJpuppiEta_);
	tree->Branch("AK8PuppiJet_SDSJpuppiMass",         &AK8PuppiJet_SDSJpuppiMass_);
	tree->Branch("AK8PuppiJet_SDSJpuppiPhi",          &AK8PuppiJet_SDSJpuppiPhi_);
	tree->Branch("AK8PuppiJet_SDSJpuppiE",            &AK8PuppiJet_SDSJpuppiE_);
	tree->Branch("AK8PuppiJet_SDSJpuppiCharge",       &AK8PuppiJet_SDSJpuppiCharge_);
	tree->Branch("AK8PuppiJet_SDSJpuppiFlavour",      &AK8PuppiJet_SDSJpuppiFlavour_);
	tree->Branch("AK8PuppiJet_SDSJpuppiCSV",          &AK8PuppiJet_SDSJpuppiCSV_);
	tree->Branch("AK8PuppiJet_SDSJpuppinb1ecf2",          &AK8PuppiJet_SDSJpuppinb1ecf2_);
	tree->Branch("AK8PuppiJet_SDSJpuppinb1ecf3",          &AK8PuppiJet_SDSJpuppinb1ecf3_);
	tree->Branch("AK8PuppiJet_SDSJpuppinb2ecf2",          &AK8PuppiJet_SDSJpuppinb2ecf2_);
	tree->Branch("AK8PuppiJet_SDSJpuppinb2ecf3",          &AK8PuppiJet_SDSJpuppinb2ecf3_);
}

void ggNtuplizer::fillAK8PUPPIJets(const edm::Event& e, const edm::EventSetup& es){
	AK8PuppiJet_Charge_.clear();
	AK8PuppiJet_Pt_.clear();
	AK8PuppiJet_En_.clear();
	AK8PuppiJet_RawPt_.clear();
	AK8PuppiJet_RawEn_.clear();
	AK8PuppiJet_Eta_.clear();
	AK8PuppiJet_Phi_.clear();
	AK8PuppiJet_Mass_.clear();
	AK8PuppiJet_tau1_.clear();
	AK8PuppiJet_tau2_.clear();
	AK8PuppiJet_tau3_.clear();
	AK8PuppiJet_CHF_.clear();
	AK8PuppiJet_NHF_.clear();
	AK8PuppiJet_CEF_.clear();
	AK8PuppiJet_NEF_.clear();
	AK8PuppiJet_NCH_.clear();
	AK8PuppiJet_NNP_.clear();
	AK8PuppiJet_MUF_.clear();
	// AK8PuppiJet_puppiMultiplicity_.clear();
	// AK8PuppiJet_neutral_puppiMultiplicity_.clear();
	AK8PuppiJet_PFid_.clear();
	AK8PuppiJet_PUPPISoftDropMass_.clear();
	// AK8PuppiJet_TrimmedMass_.clear();
	// AK8PuppiJet_FilteredMass_.clear();
	AK8PuppiJet_CHSPrunedMass_.clear();
	AK8PuppiJet_CHSSoftDropMass_.clear();
	AK8PuppiJet_pfBoostedDSVBTag_.clear();
	// AK8PuppiJet_DSVnewV4_.clear();
	AK8PuppiJet_CSV_.clear();
	AK8PuppiJet_DDBvLTags_qcd_.clear();
	AK8PuppiJet_DDBvLTags_Hbb_.clear();
	AK8PuppiJet_DDCvLTags_qcd_.clear();
	AK8PuppiJet_DDCvLTags_Hcc_.clear();
	AK8PuppiJet_DDCvBTags_Hbb_.clear();
	AK8PuppiJet_DDCvBTags_Hcc_.clear();
	AK8PuppiJet_DDBDTags_bbcsLight_.clear();
	AK8PuppiJet_DDBDTags_ccvsLight_.clear();
	AK8PuppiJet_DDBDTags_TvsQCD_.clear();
	AK8PuppiJet_DDBDTags_ZHccvsQCD_.clear();
	AK8PuppiJet_DDBDTags_WvsQCD_.clear();
	AK8PuppiJet_DDBDTags_ZHbbvsQCD_.clear();
	AK8PuppiJet_JECUnc_.clear();
	// AK8PuppiJet_L2L3corr_.clear();
	AK8PuppiJet_JER_sf_.clear();
	AK8PuppiJet_JER_sf_up_.clear();
	AK8PuppiJet_JER_sf_down_.clear();
	AK8PuppiJet_JER_ptRes_.clear();
	AK8PuppiJet_JER_phiRes_.clear();
	AK8PuppiJet_nb1ecf2_.clear();
	AK8PuppiJet_nb1ecf3_.clear();
	AK8PuppiJet_nb2ecf2_.clear();
	AK8PuppiJet_nb2ecf3_.clear();
	AK8PuppiJet_PartonFlavour_.clear();
	AK8PuppiJet_HadronFlavour_.clear();
	AK8PuppiJet_GenJetIndex_.clear();
	// AK8PuppiJet_GenPartonPID_.clear();
	// AK8PuppiJet_GenJetEn_.clear();
	// AK8PuppiJet_GenJetPt_.clear();
	// AK8PuppiJet_GenJetEta_.clear();
	// AK8PuppiJet_GenJetPhi_.clear();
	AK8PuppiJet_nSDSJpuppi_ .clear();
	AK8PuppiJet_SDSJpuppiPt_ .clear();
	AK8PuppiJet_SDSJpuppiEta_ .clear();
	AK8PuppiJet_SDSJpuppiMass_ .clear();
	AK8PuppiJet_SDSJpuppiPhi_ .clear();
	AK8PuppiJet_SDSJpuppiE_ .clear();
	AK8PuppiJet_SDSJpuppiCharge_ .clear();
	AK8PuppiJet_SDSJpuppiFlavour_.clear();
	AK8PuppiJet_SDSJpuppiCSV_ .clear();
	AK8PuppiJet_SDSJpuppinb1ecf2_.clear();
	AK8PuppiJet_SDSJpuppinb1ecf3_.clear();
	AK8PuppiJet_SDSJpuppinb2ecf2_.clear();
	AK8PuppiJet_SDSJpuppinb2ecf3_.clear();
	nAK8PuppiJet_ = 0;

	// edm::Handle<vector<reco::GenParticle> > genParticlesHandle;
	// if(doGenParticles_)e.getByToken(genParticlesCollection_, genParticlesHandle);

	edm::Handle<std::vector<reco::GenJet>> ak8genjetsHandle;
	if(doGenParticles_)e.getByToken(ak8GenJetLabel_,ak8genjetsHandle);

	edm::Handle<double> rhoHandle;
	e.getByToken(rhoLabel_, rhoHandle);
	float rho_ = *(rhoHandle.product());

	edm::ESHandle<JetCorrectorParametersCollection> AK8JetCorParColl;
	es.get<JetCorrectionsRecord>().get("AK8PFPuppi",AK8JetCorParColl);
	JetCorrectorParameters const & AK8JetCorPar = (*AK8JetCorParColl)["Uncertainty"];
	JetCorrectionUncertainty *AK8jecUnc = new JetCorrectionUncertainty(AK8JetCorPar);

	edm::Handle<edm::View<pat::Jet> > jetsAK8PUPPIhandle;
	e.getByToken(ak8JetsPUPPILabel_, jetsAK8PUPPIhandle);


	// Access JER
	JME::JetResolution jet_pt_resolution_ = JME::JetResolution::get(es, "AK8PFPuppi_pt");
	JME::JetResolution jet_phi_resolution_ = JME::JetResolution::get(es, "AK8PFPuppi_phi");
	JME::JetResolutionScaleFactor resolution_sf_ = JME::JetResolutionScaleFactor::get(es, "AK8PFPuppi");
	JME::JetParameters jet_parameters;

	for (edm::View<pat::Jet>::const_iterator ijetAK8 = jetsAK8PUPPIhandle->begin(); ijetAK8 != jetsAK8PUPPIhandle->end(); ++ijetAK8) {
		if( ijetAK8->pt() < 170.0 ) continue;

		if((!ijetAK8->isPFJet()) && (!ijetAK8->isJPTJet())) continue;

		nAK8PuppiJet_++;
		AK8PuppiJet_Charge_.push_back(ijetAK8->charge());
		AK8PuppiJet_Pt_.push_back(ijetAK8->pt());
		AK8PuppiJet_En_.push_back(ijetAK8->energy());
		AK8PuppiJet_RawPt_.push_back(ijetAK8->correctedJet("Uncorrected").pt());
		AK8PuppiJet_RawEn_.push_back(ijetAK8->correctedJet("Uncorrected").energy());
		AK8PuppiJet_Eta_.push_back(ijetAK8->eta());
		AK8PuppiJet_Phi_.push_back(ijetAK8->phi());
		AK8PuppiJet_Mass_.push_back(ijetAK8->mass());
		AK8PuppiJet_tau1_.push_back(ijetAK8->userFloat("NjettinessAK8Puppi:tau1"));
		AK8PuppiJet_tau2_.push_back(ijetAK8->userFloat("NjettinessAK8Puppi:tau2"));
		AK8PuppiJet_tau3_.push_back(ijetAK8->userFloat("NjettinessAK8Puppi:tau3"));
		AK8PuppiJet_tau4_.push_back(ijetAK8->userFloat("NjettinessAK8Puppi:tau4"));

		//jet PF Loose ID
		double NHF      = ijetAK8->neutralHadronEnergyFraction();
		double NEMF     = ijetAK8->neutralEmEnergyFraction();
		double CHF      = ijetAK8->chargedHadronEnergyFraction();
		Int_t CHM      = ijetAK8->chargedMultiplicity();
		Int_t NNP      = ijetAK8->neutralMultiplicity();
		Int_t NumConst = CHM + NNP;
		double CEMF     = ijetAK8->chargedEmEnergyFraction();
		double MUF      = ijetAK8->muonEnergyFraction();

		AK8PuppiJet_NEF_.push_back(NEMF);
		AK8PuppiJet_NCH_.push_back(CHM);
		AK8PuppiJet_NNP_.push_back(NNP);
		AK8PuppiJet_MUF_.push_back(MUF);
		// AK8PuppiJet_puppiMultiplicity_.push_back(ijetAK8->userFloat("patPuppiJetSpecificProducer:puppiMultiplicity"));
		// AK8PuppiJet_neutral_puppiMultiplicity_.push_back(ijetAK8->userFloat("patPuppiJetSpecificProducer:neutralPuppiMultiplicity"));

		AK8PuppiJet_CHF_.push_back(CHF);
		AK8PuppiJet_NHF_.push_back(NHF);
		AK8PuppiJet_CEF_.push_back(CEMF);

		//https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2017?rev=7
		bool looseJetID = false;
		bool tightJetID = false;
		Bool_t tightLeptVetoID = false;
		if (fabs(ijetAK8->eta()) <= 2.7) {
			looseJetID = (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((fabs(ijetAK8->eta())<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || fabs(ijetAK8->eta())>2.4);
			tightJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((fabs(ijetAK8->eta())<=2.4 && CHF>0 && CHM>0) || fabs(ijetAK8->eta())>2.4);
			tightLeptVetoID = (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF < 0.8) && ((fabs(ijetAK8->eta())<=2.4 && CHF>0 && CHM>0 && CEMF<0.8) || fabs(ijetAK8->eta())>2.4);
		} else if (fabs(ijetAK8->eta()) <= 3.0) {
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
		AK8PuppiJet_PFid_.push_back(jetIDdecision);

		AK8PuppiJet_PUPPISoftDropMass_.push_back(ijetAK8->userFloat("ak8PFJetsPuppiSoftDropMass"));
		// AK8PuppiJet_TrimmedMass_.push_back(ijetAK8->userFloat("ak8PFJetsCHSTrimmedLinks"));;
		// AK8PuppiJet_FilteredMass_.push_back(ijetAK8->userFloat("ak8PFJetsCHSFilteredLinks"));;
		AK8PuppiJet_CHSPrunedMass_.push_back(ijetAK8->userFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSPrunedMass"));
		AK8PuppiJet_CHSSoftDropMass_.push_back(ijetAK8->userFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSSoftDropMass"));
		AK8PuppiJet_pfBoostedDSVBTag_.push_back(ijetAK8->bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags"));
		// AK8PuppiJet_DSVnewV4_.push_back(ijetAK8->bDiscriminator("newV4pfBoostedDoubleSecondaryVertexAK8BJetTags"));
		AK8PuppiJet_CSV_.push_back(ijetAK8->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
		AK8PuppiJet_DDBvLTags_qcd_.push_back(ijetAK8->bDiscriminator("pfMassIndependentDeepDoubleBvLJetTags:probQCD"));
		AK8PuppiJet_DDBvLTags_Hbb_.push_back(ijetAK8->bDiscriminator("pfMassIndependentDeepDoubleBvLJetTags:probHbb"));
		AK8PuppiJet_DDCvLTags_qcd_.push_back(ijetAK8->bDiscriminator("pfMassIndependentDeepDoubleCvLJetTags:probQCD"));
		AK8PuppiJet_DDCvLTags_Hcc_.push_back(ijetAK8->bDiscriminator("pfMassIndependentDeepDoubleCvLJetTags:probHcc"));
		AK8PuppiJet_DDCvBTags_Hbb_.push_back(ijetAK8->bDiscriminator("pfMassIndependentDeepDoubleCvBJetTags:probHbb"));
		AK8PuppiJet_DDCvBTags_Hcc_.push_back(ijetAK8->bDiscriminator("pfMassIndependentDeepDoubleCvBJetTags:probHcc"));
		AK8PuppiJet_DDBDTags_bbcsLight_.push_back(ijetAK8->bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:bbvsLight"));
		AK8PuppiJet_DDBDTags_ccvsLight_.push_back(ijetAK8->bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ccvsLight"));
		AK8PuppiJet_DDBDTags_TvsQCD_.push_back(ijetAK8->bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:TvsQCD"));
		AK8PuppiJet_DDBDTags_ZHccvsQCD_.push_back(ijetAK8->bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZHccvsQCD"));
		AK8PuppiJet_DDBDTags_WvsQCD_.push_back(ijetAK8->bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:WvsQCD"));
		AK8PuppiJet_DDBDTags_ZHbbvsQCD_.push_back(ijetAK8->bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZHbbvsQCD"));


		// JEC https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections
		if (std::abs(ijetAK8->eta()) < 5.2){
			AK8jecUnc->setJetEta(ijetAK8->eta());
			// here you must use the CORRECTED jet pt
			AK8jecUnc->setJetPt(ijetAK8->pt());
			AK8PuppiJet_JECUnc_.push_back(AK8jecUnc->getUncertainty(true));
		} else {
			AK8PuppiJet_JECUnc_.push_back(-1.);
		}

		jet_parameters.setJetPt(ijetAK8->pt()).setJetEta(ijetAK8->eta()).setRho(rho_);
		Float_t _jer_ptRes = jet_pt_resolution_.getResolution(jet_parameters);
		Float_t _jer_phiRes = jet_phi_resolution_.getResolution(jet_parameters);
		Float_t _jer_sf = resolution_sf_.getScaleFactor(jet_parameters);
		Float_t _jersf_up = resolution_sf_.getScaleFactor(jet_parameters, Variation::UP);
		Float_t _jersf_down = resolution_sf_.getScaleFactor(jet_parameters, Variation::DOWN);

		AK8PuppiJet_JER_ptRes_.push_back(_jer_ptRes);
		AK8PuppiJet_JER_phiRes_.push_back(_jer_phiRes);
		AK8PuppiJet_JER_sf_.push_back(_jer_sf);
		AK8PuppiJet_JER_sf_up_.push_back(_jersf_up);
		AK8PuppiJet_JER_sf_down_.push_back(_jersf_down);

      	// AK8PuppiJet_L2L3corr_.push_back(ijetAK8->);

		AK8PuppiJet_nb1ecf2_.push_back(ijetAK8->userFloat("ak8PFJetsPuppiSoftDropValueMap:nb1AK8PuppiSoftDropN2"));
		AK8PuppiJet_nb1ecf3_.push_back(ijetAK8->userFloat("ak8PFJetsPuppiSoftDropValueMap:nb1AK8PuppiSoftDropN3"));
		AK8PuppiJet_nb2ecf2_.push_back(ijetAK8->userFloat("ak8PFJetsPuppiSoftDropValueMap:nb2AK8PuppiSoftDropN2"));
		AK8PuppiJet_nb2ecf3_.push_back(ijetAK8->userFloat("ak8PFJetsPuppiSoftDropValueMap:nb2AK8PuppiSoftDropN3"));

		AK8PuppiJet_PartonFlavour_.push_back(ijetAK8->partonFlavour());
		AK8PuppiJet_HadronFlavour_.push_back(ijetAK8->hadronFlavour());

		if(doGenParticles_){
			Int_t _genJetIndex = -99;
			if (ijetAK8->genJet() && ak8genjetsHandle.isValid()) {
				_genJetIndex = std::distance(ak8genjetsHandle->begin(), (vector<reco::GenJet>::const_iterator) ijetAK8->genJet());
			}
			AK8PuppiJet_GenJetIndex_.push_back(_genJetIndex);

			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// Exception Message:
			// RefCore: A request to resolve a reference to a product of type 'std::vector<reco::GenParticle>' with ProductID '3:2895'
			// can not be satisfied because the product cannot be found.
			// Probably the branch containing the product is not stored in the input file.
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// std::vector<Short_t> __genPartonPIDs;
			// const reco::GenParticleRefVector & partonsRefVector = ijetAK8->jetFlavourInfo().getPartons();
			// if(partonsRefVector.empty () == 0){

			// 	for(UInt_t i = 0; i < partonsRefVector.size(); i++){
			// 		const reco::GenParticle* part = partonsRefVector.at(i).get();
			// 		cout<<part->pdgId()<<endl;
			// 	}
				// cout<<partonsRefVector.size()<<endl;
				// for(reco::GenParticleRefVector::const_iterator iParton = partonsRefVector.begin(); iParton != partonsRefVector.end(); ++iParton){

				// 	const reco::GenParticle& part = **iParton;
				// 	cout<<part.pdgId()<<endl;

				// 	const reco::GenParticle * _part = &(**iParton);

				// 	Short_t partGenPos_ = (_part == nullptr) ? -999 : std::distance(genParticlesHandle->begin(), (vector<reco::GenParticle>::const_iterator) _part);

				// 	// __genPartonPIDs.push_back((*iParton)->pdgId());
				// 	cout<<" ak8jet parton pdgID "<<_part->pdgId()<<"  index  "<< partGenPos_<<endl;
				// }
			// }
			// AK8PuppiJet_GenPartonPID_.push_back(__genPartonPIDs);
		}

		if (dumpSoftDrop_) {
			int nsubjets = 0;
			std::vector<float> vecSDSJpuppicsv ;
			std::vector<float> vecSDSJpuppipt ;
			std::vector<float> vecSDSJpuppieta ;
			std::vector<float> vecSDSJpuppimass ;
			std::vector<float> vecSDSJpuppiphi ;
			std::vector<float> vecSDSJpuppie ;
			std::vector<Char_t> vecSDSJpuppicharge ;
			std::vector<int > vecSDSJpuppiflavour;
			std::vector<float>  vecSDSJpuppinb1ecf2_;
			std::vector<float>  vecSDSJpuppinb1ecf3_;
			std::vector<float>  vecSDSJpuppinb2ecf2_;
			std::vector<float>  vecSDSJpuppinb2ecf3_;
			auto const & sdSubjets = ijetAK8->subjets("SoftDropPuppi");
			for ( auto const & SDSJpuppi : sdSubjets ) {
				nsubjets++;
				vecSDSJpuppipt.push_back(SDSJpuppi->pt());
				vecSDSJpuppieta.push_back(SDSJpuppi->eta());
				vecSDSJpuppimass.push_back(SDSJpuppi->mass());
				vecSDSJpuppiphi.push_back(SDSJpuppi->phi());
				vecSDSJpuppie.push_back(SDSJpuppi->energy());
				vecSDSJpuppiflavour.push_back(SDSJpuppi->partonFlavour());
				vecSDSJpuppicharge.push_back(SDSJpuppi->charge());
				vecSDSJpuppicsv.push_back(SDSJpuppi->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") );
				vecSDSJpuppinb1ecf2_.push_back(SDSJpuppi->userFloat("nb1AK8PuppiSoftDropSubjets:ecfN2"));
				vecSDSJpuppinb1ecf3_.push_back(SDSJpuppi->userFloat("nb1AK8PuppiSoftDropSubjets:ecfN3"));
				vecSDSJpuppinb2ecf2_.push_back(SDSJpuppi->userFloat("nb2AK8PuppiSoftDropSubjets:ecfN2"));
				vecSDSJpuppinb2ecf3_.push_back(SDSJpuppi->userFloat("nb2AK8PuppiSoftDropSubjets:ecfN3"));
			}

			AK8PuppiJet_nSDSJpuppi_.push_back(nsubjets);
			AK8PuppiJet_SDSJpuppiPt_.push_back(vecSDSJpuppipt);
			AK8PuppiJet_SDSJpuppiEta_.push_back(vecSDSJpuppieta);
			AK8PuppiJet_SDSJpuppiPhi_.push_back(vecSDSJpuppiphi);
			AK8PuppiJet_SDSJpuppiMass_.push_back(vecSDSJpuppimass);
			AK8PuppiJet_SDSJpuppiE_.push_back(vecSDSJpuppie);
			AK8PuppiJet_SDSJpuppiCharge_.push_back(vecSDSJpuppicharge);
			AK8PuppiJet_SDSJpuppiFlavour_.push_back(vecSDSJpuppiflavour);
			AK8PuppiJet_SDSJpuppiCSV_.push_back(vecSDSJpuppicsv);
			AK8PuppiJet_SDSJpuppinb1ecf2_.push_back(vecSDSJpuppinb1ecf2_);
			AK8PuppiJet_SDSJpuppinb1ecf3_.push_back(vecSDSJpuppinb1ecf3_);
			AK8PuppiJet_SDSJpuppinb2ecf2_.push_back(vecSDSJpuppinb2ecf2_);
			AK8PuppiJet_SDSJpuppinb2ecf3_.push_back(vecSDSJpuppinb2ecf3_);
		}
	}
	delete AK8jecUnc;
};
