import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process('ggKit',eras.Run2_2017)

##########################################################################################################################
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff" )
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_v11', '')
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.source = cms.Source("PoolSource",
                            #fileNames = cms.untracked.vstring('file:F8DDFDC7-8AD6-E711-BCA2-4C79BA1811CB.root')
                            #fileNames = cms.untracked.vstring('root://xrootd-cms.infn.it//store/data/Run2017B/SinglePhoton/MINIAOD/17Nov2017-v1/60000/3011B1EA-0BE7-E711-8D8B-3417EBE61338.root')
                            
                            #fileNames = cms.untracked.vstring('root://xrootd-cms.infn.it//store/data/Run2017B/SinglePhoton/AOD/17Nov2017-v1/60000/CCC5B4B1-DBE7-E711-9EFD-7845C4FC3602.root')
                            #fileNames = cms.untracked.vstring('root://xrootd-cms.infn.it//store/data/Run2017B/SinglePhoton/AOD/17Nov2017-v1/60000/CCC5B4B1-DBE7-E711-9EFD-7845C4FC3602.root')
                            fileNames = cms.untracked.vstring('file:CCC5B4B1-DBE7-E711-9EFD-7845C4FC3602.root')

                            #fileNames = cms.untracked.vstring('root://xrootd-cms.infn.it//store/data/Run2017B/SinglePhoton/AOD/23Jun2017-v1/10000/3AB2A9ED-5259-E711-9E5D-FA163E65C977.root')
)
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.load( "PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff" )
process.load( "PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff" )
process.load( "PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff" )
process.TFileService = cms.Service("TFileService", fileName = cms.string('anTGCtree_data.root'))
##########################################################################################################################

#process.Tracer = cms.Service("Tracer")

##########################################################################################################################
### fix a bug in the ECAL-Tracker momentum combination when applying the scale and smearing
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runEnergyCorrections=True,
                       runVID=True,
                       isMiniAOD=False,
                       era='2017-Nov17ReReco',
                       eleIDModules=['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V2_cff',
                                     'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff',
                                     'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V2_cff',
                                     'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V2_cff'],
                       phoIDModules=['RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V2_cff',
                       #phoIDModules=['RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V1p1_cff',
                                     'RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Fall17_94X_V2_cff']
)
##########################################################################################################################

#from JMEAnalysis.JetToolbox.jetToolbox_cff import *
#jetToolbox( process, 'ak4', 'ak4PFJetsCHS', 'out', miniAOD= False, addSoftDrop=True, addSoftDropSubjets=True, addNsub=True, addPUJetID=True, JETCorrPayload='AK4PFchs', JETCorrLevels=['L1FastJet','L2Relative', 'L3Absolute','L2L3Residual'] )
#jetToolbox( process, 'ak4', 'ak4PFJetsCHS', 'noOutput', miniAOD= False, runOnMC=False, addSoftDrop=True, addSoftDropSubjets=False, addNsub=False, addPUJetID=True, JETCorrPayload='AK4PFchs' )
#jetToolbox( process, 'ak8', 'ak8PFJetsCHS', 'out', miniAOD= False, addSoftDrop=True, addSoftDropSubjets=True, addNsub=True )
#process.ggNtuplizer.dumpSoftDrop= cms.bool(False)



##########################################################################################################################
from PhysicsTools.PatAlgos.tools.coreTools import *
runOnData( process,  names=['Photons', 'Electrons','Muons','Taus','Jets'], outputModules = [] )
removeMCMatching(process, names=['All'], outputModules=[])

##########################################################################################################################


##########################################################################################################################
### reduce effect of high eta EE noise on the PF MET measurement
#from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
#runMetCorAndUncFromMiniAOD (
#        process,
#        isData = True, # false for MC
#        fixEE2017 = True,
#        fixEE2017Params = {'userawPt': True, 'ptThreshold':50.0, 'minEtaThreshold':2.65, 'maxEtaThreshold': 3.139} ,
#        postfix = "ModifiedMET"
#)
##########################################################################################################################


##########################################################################################################################
### include jetToolbox to add various jets
#from JMEAnalysis.JetToolbox.jetToolbox_cff import jetToolbox
#jetToolbox( process, 'ak4', 'ak4JetSubs', 'noOutput',
#        runOnMC=False,
#        updateCollection='slimmedJets',
#        JETCorrPayload='AK4PFchs',
#        postFix='updated'
#        )

### ak 0.8 PUPPI jets
#jetToolbox( process, 'ak8', 'ak8PUPPIJetToolbox', 'noOutput',
#            runOnMC=False,
#            PUMethod='PUPPI',
#            updateCollection='slimmedJetsAK8',
#            updateCollectionSubjets='slimmedJetsAK8PFPuppiSoftDropPacked',
#            JETCorrPayload = 'AK8PFPuppi'
#            )
##########################################################################################################################


##########################################################################################################################
process.load("ggAnalysis.ggNtuplizer.ggNtuplizer_AOD_cfi")
process.ggNtuplizer.year=cms.int32(2017)
process.ggNtuplizer.doGenParticles=cms.bool(False)
process.ggNtuplizer.dumpPFPhotons=cms.bool(False)
process.ggNtuplizer.dumpHFElectrons=cms.bool(False)
process.ggNtuplizer.dumpJets=cms.bool(True)
process.ggNtuplizer.dumpAK8Jets=cms.bool(False)
process.ggNtuplizer.dumpSoftDrop= cms.bool(False)
process.ggNtuplizer.dumpTaus=cms.bool(False)
#process.ggNtuplizer.pfMETLabel=cms.InputTag("slimmedMETsModifiedMET")
#process.ggNtuplizer.ak4PFJetsCHSSrc=cms.InputTag("selectedPatJetsAK4PFCHSupdated")
#process.ggNtuplizer.ak8JetsPUPPISrc=cms.InputTag("selectedPatJetsAK8PFPUPPI")
#process.ggNtuplizer.patTriggerResults = cms.InputTag("TriggerResults", "", "RECO")
#process.ggNtuplizer.triggerEvent=cms.InputTag("slimmedPatTrigger")
##########################################################################################################################

####L1 bits
process.ggNtuplizer.triggerSelection = cms.string("L1_BptxXOR")
process.ggNtuplizer.triggerSelectionPB = cms.string("L1_BptxPlus")
process.ggNtuplizer.triggerSelectionMB = cms.string("L1_BptxMinus")
process.ggNtuplizer.triggerSelectionZB = cms.string("L1_ZeroBias")

process.ggNtuplizer.triggerConfiguration =  cms.PSet(
    hltResults = cms.InputTag('TriggerResults','','HLT'),
    l1tResults = cms.InputTag('gtStage2Digis'),
    daqPartitions = cms.uint32(1),
    #l1tIgnoreMask = cms.bool( False ),
    #l1techIgnorePrescales = cms.bool( True ),
    l1tIgnoreMaskAndPrescale = cms.bool( True ),
    throw  = cms.bool( True )
)




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
    #EcalRecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
    EcalRecHitSource = cms.InputTag("reducedEcalRecHitsEE"),
    ecalMinEt        = cms.double(50.),
    baddetEcal    = baddetEcallist,
    taggingMode = cms.bool(True),
    debug = cms.bool(False)
    )





#process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')
process.load('RecoMET.METFilters.globalSuperTightHalo2016Filter_cfi')
process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.load('CommonTools.RecoAlgos.HBHENoiseFilter_cfi')

process.load('RecoMET.METFilters.primaryVertexFilter_cfi')
##change the name from primaryVertexFilter  ---> some kind of clash happens
process.goodVertexFilter = cms.EDFilter(
    "GoodVertexFilter",
    vertexCollection = cms.InputTag('offlinePrimaryVertices'),
    minimumNDOF = cms.uint32(4) ,
    maxAbsZ = cms.double(24),
    maxd0 = cms.double(2)
)



process.HBHENoiseFilter = cms.EDFilter(
    'BooleanFlagFilter',
    #inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
    inputLabel = cms.InputTag('HBHENoiseFilterResultProducer'),
    reverseDecision = cms.bool(False)
)

###https://cmssdt.cern.ch/lxr/source/RecoMET/METFilters/test/filters_cfg.py?v=CMSSW_11_1_X_2020-03-03-2300
##https://hypernews.cern.ch/HyperNews/CMS/get/met/582.html

#process.HBHENoiseFilter.taggingMode = cms.bool(False)
process.ggNtuplizer.ecalBadCalibFilter = cms.InputTag("ecalBadCalibReducedMINIAODFilter")
process.ggNtuplizer.passedVertexFilter = cms.InputTag("goodVertexFilter")
process.ggNtuplizer.passedEcalDeadCell = cms.InputTag("EcalDeadCellTriggerPrimitiveFilter")
process.ggNtuplizer.passedGlobalHalo = cms.InputTag("globalSuperTightHalo2016Filter")
process.ggNtuplizer.passedeeBadScFilter = cms.InputTag("eeBadScFilter")
#process.ggNtuplizer.passedHBHENoiseIsoFilter = cms.InputTag("HBHENoiseFilter")
process.ggNtuplizer.passedHBHENoiseFilter = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult')

##########################################################################################################################

process.load("PhysicsTools.PatAlgos.producersLayer1.ootPhotonProducer_cff")
process.patOOTPhotons.ecalPFClusterIsoMap = cms.InputTag("ootPhotonEcalPFClusterIsolationProducer")
process.patOOTPhotons.hcalPFClusterIsoMap = cms.InputTag("ootPhotonHcalPFClusterIsolationProducer")



##########################################################################################################################
process.p = cms.Path(
    #process.goodVertexFilter * 
    process.EcalDeadCellTriggerPrimitiveFilter * 
    process.globalSuperTightHalo2016Filter * 
    process.eeBadScFilter *
    process.HBHENoiseFilterResultProducer *
    #process.HBHENoiseFilter *
    process.ecalBadCalibReducedMINIAODFilter *

    process.patDefaultSequence *
    #process.fullPatMetSequenceModifiedMET *
    process.egammaPostRecoSeq 
    * process.ggNtuplizer
    )

#process.endjob_step = cms.EndPath(process.endOfProcess)
'''
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('myOutputFile.root'),
                               SelectEvents = cms.untracked.PSet( 
                                   SelectEvents = cms.vstring("p")
                            )
                           )



process.e = cms.EndPath(process.out)
'''
#print process.dumpPython()
##########################################################################################################################
