import FWCore.ParameterSet.Config as cms
#from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector
#from CMGTools.External.pujetidproducer_cfi import stdalgos, chsalgos

ggNtuplizer = cms.EDAnalyzer("ggNtuplizer",
                 doGenParticles       = cms.bool(True),
                             #doGenJets           = cms.bool(True),
                 runOnParticleGun     = cms.bool(False),
                 runOnSherpa          = cms.bool(False),
                 dumpPFPhotons        = cms.bool(True),
                 dumpJets             = cms.bool(False),
                 dumpAK8Jets          = cms.bool(False),
                 dumpTaus             = cms.bool(False),
                 dumpPDFSystWeight    = cms.bool(False),
                 dumpHFElectrons      = cms.bool(True),
                 development          = cms.bool(False),
                 addFilterInfoAOD     = cms.bool(True),
                             addFilterInfoMINIAOD = cms.bool(False),
                 doNoHFMET            = cms.bool(False),
                 getECALprefiringWeights   = cms.bool(False),

                 doOOTphotons   = cms.bool(True),

                 year                 = cms.int32(2017),

                 trgFilterDeltaPtCut  = cms.double(0.5),
                 trgFilterDeltaRCut   = cms.double(0.3),

                 beamHaloSummary        = cms.InputTag("BeamHaloSummary"),
                             
                 triggerEvent         = cms.InputTag("hltTriggerSummaryAOD"),
                 triggerResults       = cms.InputTag("TriggerResults", "", "HLT"),
                 patTriggerResults    = cms.InputTag("TriggerResults", "", "HLT"),
                 #patTriggerResults    = cms.InputTag("TriggerResults", "", "RECO"),
                 genParticleSrc       = cms.InputTag("genParticles"),
                 generatorLabel       = cms.InputTag("generator"),
                 LHEEventLabel        = cms.InputTag("externalLHEProducer"),
                 newParticles         = cms.vint32(1000006, 1000021, 1000022, 1000024, 1000025, 1000039, 3000001, 3000002, 35),
                 pileupCollection     = cms.InputTag("addPileupInfo"),
                             VtxLabel             = cms.InputTag("offlinePrimaryVertices"), ###check if offlinePrimaryVerticesWithBS is needed
                 VtxBSLabel           = cms.InputTag("offlinePrimaryVerticesWithBS"),
                 rhoLabel             = cms.InputTag("fixedGridRhoFastjetAll"),
                 rhoCentralLabel      = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
                 pfMETLabel           = cms.InputTag("pfMet"),
                 electronSrc          = cms.InputTag("gedGsfElectrons"),
                 #calibelectronSrc     = cms.InputTag("calibratedPatElectrons"),
                 calibelectronSrc     = cms.InputTag("gedGsfElectrons"),
                 photonSrc            = cms.InputTag("gedPhotons"),
                 photonOOTSrc            = cms.InputTag("ootPhotons"),
                 #calibphotonSrc       = cms.InputTag("calibratedPatPhotons"),
                 calibphotonSrc       = cms.InputTag("gedPhotons"),
                 muonSrc              = cms.InputTag("globalMuons"),
                 gsfTrackSrc          = cms.InputTag("electronGsfTracks"),
                             ebReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
                 eeReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
                 esReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsES"),
                             ecalSCcollection       = cms.InputTag("particleFlowSuperClusterECAL", "particleFlowSuperClusterECALBarrel"),
                             #ecalSCcollectionEB       = cms.InputTag("particleFlowSuperClusterECAL", "particleFlowSuperClusterECALBarrel"),
                             #ecalSCcollectionEE       = cms.InputTag("particleFlowSuperClusterECAL", "particleFlowSuperClusterECALEndcapWithPreshower"),
                             ecalSCOOTcollection       = cms.InputTag("", "reducedOOTSuperClusters"),
                             #ecalSCOOTcollectionEB       = cms.InputTag("particleFlowSuperClusterOOTECAL", "particleFlowSuperClusterOOTECALBarrel"),
                             #ecalSCOOTcollectionEE       = cms.InputTag("particleFlowSuperClusterOOTECAL", "particleFlowSuperClusterOOTECALEndcapWithPreshower"),
                             recoPhotonSrc             = cms.InputTag("gedPhotonCore"),
                 TrackLabel                = cms.InputTag("generalTracks"),
                 gsfElectronLabel          = cms.InputTag("gedGsfElectrons"),
                 PFAllCandidates           = cms.InputTag("particleFlow"),
                 #ak4PFJetsCHS                 = cms.InputTag("updatedJets"),
                 ak4PFJetsCHSSrc           = cms.InputTag("ak4PFJetsCHS"),
                 ak4PFJetsPUPPISrc           = cms.InputTag(""), #####----> not available in AOD
                 ak4PFJetsCHSGenJetLabel      = cms.InputTag("ak4GenJets"),
                 ak8GenJetLabel                 = cms.InputTag("ak8GenJets"),
                             ak8JetsPUPPISrc                 = cms.InputTag("slimmedJetsAK8"),  #####----> not available in AOD 
                 tauSrc                    = cms.InputTag("hpsPFTauProducer"),

                 packedPFCands             = cms.InputTag("particleFlow"),
                             #BadChargedCandidateFilter = cms.InputTag("BadChargedCandidateFilter"), ####---> not available in AOD
                             #BadPFMuonFilter           = cms.InputTag("BadPFMuonFilter"),     #### ----> not available in AOD
 )
