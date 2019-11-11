import FWCore.ParameterSet.Config as cms

process = cms.Process('ggKit')

##########################################################################################################################
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff" )
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_v10')
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('file:2016SinglePhotonMINIAOD.root')
                            )
process.load( "PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff" )
process.load( "PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff" )
process.load( "PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff" )
process.TFileService = cms.Service("TFileService", fileName = cms.string('ggtree_2016data.root'))
##########################################################################################################################


##########################################################################################################################
### fix a bug in the ECAL-Tracker momentum combination when applying the scale and smearing
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runVID=True,
                       era='2016-Legacy',
                       )
##########################################################################################################################


##########################################################################################################################
#from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *
runOnData( process,  names=['Photons', 'Electrons','Muons','Taus','Jets'], outputModules = [] )
##########################################################################################################################


##########################################################################################################################
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD (
        process,
        isData = True,
        postfix = "ModifiedMET"
)
##########################################################################################################################


##########################################################################################################################
### include jetToolbox to add various jets
from JMEAnalysis.JetToolbox.jetToolbox_cff import jetToolbox
jetToolbox( process, 'ak4', 'ak4JetSubs', 'noOutput',
        runOnMC=False,
        updateCollection='slimmedJets',
        JETCorrPayload='AK4PFchs',
        postFix='updated'
        )

### ak 0.8 PUPPI jets
jetToolbox( process, 'ak8', 'ak8PUPPIJetToolbox', 'noOutput',
            runOnMC=False,
            PUMethod='PUPPI',
            updateCollection='slimmedJetsAK8',
            updateCollectionSubjets='slimmedJetsAK8PFPuppiSoftDropPacked',
            JETCorrPayload = 'AK8PFPuppi'
            )
##########################################################################################################################


##########################################################################################################################
process.load("ggAnalysis.ggNtuplizer.ggNtuplizer_miniAOD_cfi")
process.ggNtuplizer.year=cms.int32(2016)
process.ggNtuplizer.doGenParticles=cms.bool(False)
process.ggNtuplizer.dumpPFPhotons=cms.bool(False)
process.ggNtuplizer.dumpHFElectrons=cms.bool(False)
process.ggNtuplizer.dumpJets=cms.bool(True)
process.ggNtuplizer.dumpAK8Jets=cms.bool(True)
process.ggNtuplizer.dumpSoftDrop= cms.bool(True)
process.ggNtuplizer.dumpTaus=cms.bool(False)
process.ggNtuplizer.pfMETLabel=cms.InputTag("slimmedMETsModifiedMET")
process.ggNtuplizer.ak4PFJetsCHSSrc=cms.InputTag("selectedPatJetsAK4PFCHSupdated")
process.ggNtuplizer.ak8JetsPUPPISrc=cms.InputTag("selectedPatJetsAK8PFPUPPI")
process.ggNtuplizer.patTriggerResults = cms.InputTag("TriggerResults", "", "RECO")
process.ggNtuplizer.triggerEvent=cms.InputTag("slimmedPatTrigger")
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
    process.ecalBadCalibReducedMINIAODFilter *
    process.fullPatMetSequenceModifiedMET *
    process.egammaPostRecoSeq *
    process.ggNtuplizer
    )

#print process.dumpPython()
##########################################################################################################################