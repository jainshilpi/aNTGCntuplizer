import FWCore.ParameterSet.Config as cms

process = cms.Process('ggKit')


##########################################################################################################################
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_v17')
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 200
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
    'file:/afs/cern.ch/work/m/mwadud/private/naTGC/beamHalo/CMSSW_9_4_13/src/GeneratorInterface/BeamHaloGenerator/test/beamHaloMINIAODSIM.root'
    ))

process.load( "PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff" )
process.load( "PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff" )

process.TFileService = cms.Service("TFileService", fileName = cms.string('anTGCtree_mc.root'))
##########################################################################################################################



##########################################################################################################################
### ECAL prefiring correction
from PhysicsTools.PatUtils.l1ECALPrefiringWeightProducer_cfi import l1ECALPrefiringWeightProducer
process.prefiringweight = l1ECALPrefiringWeightProducer.clone(
    DataEra = cms.string("2017BtoF"),
    UseJetEMPt = cms.bool(False),
    PrefiringRateSystematicUncty = cms.double(0.2),
    SkipWarnings = False)
##########################################################################################################################



##########################################################################################################################
### fix a bug in the ECAL-Tracker momentum combination when applying the scale and smearing
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
     runVID=True,
     era='2017-Nov17ReReco',
     eleIDModules=['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V2_cff',
     'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff',
     'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V2_cff',
     'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V2_cff'],
     phoIDModules=['RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V2_cff',
     'RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Fall17_94X_V2_cff']
 )
##########################################################################################################################



##########################################################################################################################
### reduce effect of high eta EE noise on the PF MET measurement
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD (process,
        isData = False, # false for MC
        fixEE2017 = True,
        fixEE2017Params = {'userawPt': True, 'ptThreshold':50.0, 'minEtaThreshold':2.65, 'maxEtaThreshold': 3.139} ,
        postfix = "ModifiedMET"
        )
##########################################################################################################################


##########################################################################################################################
### include jetToolbox to add various jets
from JMEAnalysis.JetToolbox.jetToolbox_cff import jetToolbox

jetToolbox( process, 'ak4', 'ak4JetSubs', 'noOutput',
  updateCollection='slimmedJets',
  JETCorrPayload='AK4PFchs',
  JETCorrLevels = ['L1FastJet', 'L2Relative', 'L3Absolute'],
  postFix='updated'
  )

jetToolbox( process, 'ak8', 'ak8PUPPIJetToolbox', 'noOutput',
  PUMethod='PUPPI',
  updateCollection='slimmedJetsAK8',
  updateCollectionSubjets='slimmedJetsAK8PFPuppiSoftDropPacked',
  JETCorrPayload = 'AK8PFPuppi'
            # subJETCorrPayload='AK4PFPuppi'
            # addPruning=True,
            # addSoftDrop=True
            # addTrimming=True,
            # addFiltering=True,
            # addNsub=True,
            # maxTau=3,
            # addQGTagger=True,
            # addPUJetID=True,
            # addEnergyCorrFunc=True
            )
##########################################################################################################################


##########################################################################################################################
process.load("ggAnalysis.ggNtuplizer.ggNtuplizer_miniAOD_cfi")
process.ggNtuplizer.year=cms.int32(2017)
process.ggNtuplizer.doGenParticles=cms.bool(True)
process.ggNtuplizer.dumpPFPhotons=cms.bool(False)
process.ggNtuplizer.dumpHFElectrons=cms.bool(False)
process.ggNtuplizer.dumpJets=cms.bool(True)
process.ggNtuplizer.dumpAK8Jets=cms.bool(True)
process.ggNtuplizer.dumpSoftDrop= cms.bool(True)
process.ggNtuplizer.dumpTaus=cms.bool(False)
process.ggNtuplizer.triggerEvent=cms.InputTag("slimmedPatTrigger")
process.ggNtuplizer.ak4PFJetsCHSSrc=cms.InputTag("selectedPatJetsAK4PFCHSupdated")
process.ggNtuplizer.ak4PFJetsCHSGenJetLabel      = cms.InputTag("selectedPatJetsAK4PFCHSupdated","genJets","ggKit")
process.ggNtuplizer.ak8GenJetLabel      = cms.InputTag("selectedPatJetsAK8PFPUPPI","genJets","ggKit")
process.ggNtuplizer.ak8JetsPUPPISrc=cms.InputTag("selectedPatJetsAK8PFPUPPI")
process.ggNtuplizer.runOnSherpa = cms.bool(True)
# process.ggNtuplizer.patTriggerResults = cms.InputTag("TriggerResults", "", "PAT")
process.ggNtuplizer.patTriggerResults = cms.InputTag("TriggerResults", "", "RECO")
# process.ggNtuplizer.triggerEvent=cms.InputTag("slimmedPatTrigger", "", "RECO")
process.ggNtuplizer.getECALprefiringWeights = cms.bool(True)
##########################################################################################################################


##########################################################################################################################
process.load('RecoMET.METFilters.ecalBadCalibFilter_cfi')

baddetEcallist = cms.vuint32(
  [872439604,872422825,872420274,872423218,
  872423215,872416066,872435036,872439336,
  872420273,872436907,872420147,872439731,
  872436657,872420397,872439732,872439339,
  872439603,872422436,872439861,872437051,
  872437052,872420649,872422436,872421950,
  872437185,872422564,872421566,872421695,
  872421955,872421567,872437184,872421951,
  872421694,872437056,872437057,872437313])

process.ecalBadCalibReducedMINIAODFilter = cms.EDFilter(
  "EcalBadCalibFilter",
  EcalRecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
  ecalMinEt        = cms.double(50.),
  baddetEcal    = baddetEcallist,
  taggingMode = cms.bool(True),
  debug = cms.bool(False)
  )

process.ggNtuplizer.ecalBadCalibFilter = cms.InputTag("ecalBadCalibReducedMINIAODFilter")
##########################################################################################################################

##########################################################################################################################
process.p = cms.Path(
  process.prefiringweight *
  process.ecalBadCalibReducedMINIAODFilter *
  process.fullPatMetSequenceModifiedMET *
  process.egammaPostRecoSeq *
  process.ggNtuplizer
  )

# print process.dumpPython()
##########################################################################################################################
