import FWCore.ParameterSet.Config as cms
#from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector
#from CMGTools.External.pujetidproducer_cfi import stdalgos, chsalgos

ggNtuplizer = cms.EDAnalyzer("ggNtuplizer",
                 doGenParticles       = cms.bool(True),
                 doGenJets           = cms.bool(True),
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
                 addFilterInfoMINIAOD = cms.bool(True),
                 doNoHFMET            = cms.bool(False),

                 year                 = cms.int32(2017),

                 trgFilterDeltaPtCut  = cms.double(0.5),
                 trgFilterDeltaRCut   = cms.double(0.3),

                 triggerEvent         = cms.InputTag("slimmedPatTrigger", "", ""),
                 triggerResults       = cms.InputTag("TriggerResults", "", "HLT"),
                 patTriggerResults    = cms.InputTag("TriggerResults", "", "PAT"),
                 #patTriggerResults    = cms.InputTag("TriggerResults", "", "RECO"),
                 genParticleSrc       = cms.InputTag("prunedGenParticles"),
                 generatorLabel       = cms.InputTag("generator"),
                 LHEEventLabel        = cms.InputTag("externalLHEProducer"),
                 newParticles         = cms.vint32(1000006, 1000021, 1000022, 1000024, 1000025, 1000039, 3000001, 3000002, 35),
                 pileupCollection     = cms.InputTag("slimmedAddPileupInfo"),
                 VtxLabel             = cms.InputTag("offlineSlimmedPrimaryVertices"),
                 VtxBSLabel           = cms.InputTag("offlinePrimaryVerticesWithBS"),
                 rhoLabel             = cms.InputTag("fixedGridRhoFastjetAll"),
                 rhoCentralLabel      = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
                 pfMETLabel           = cms.InputTag("slimmedMETs"),
                 electronSrc          = cms.InputTag("slimmedElectrons"),
                 #calibelectronSrc     = cms.InputTag("calibratedPatElectrons"),
                 calibelectronSrc     = cms.InputTag("slimmedElectrons"),
                 photonSrc            = cms.InputTag("slimmedPhotons"),
                 #calibphotonSrc       = cms.InputTag("calibratedPatPhotons"),
                 calibphotonSrc       = cms.InputTag("slimmedPhotons"),
                 muonSrc              = cms.InputTag("slimmedMuons"),
                 gsfTrackSrc          = cms.InputTag("reducedEgamma", "reducedGsfTracks"),
                 ebReducedRecHitCollection = cms.InputTag("reducedEgamma", "reducedEBRecHits"),
                 eeReducedRecHitCollection = cms.InputTag("reducedEgamma", "reducedEERecHits"),
                 esReducedRecHitCollection = cms.InputTag("reducedEgamma", "reducedESRecHits"),
                 ecalSCcollection       = cms.InputTag("reducedEgamma", "reducedSuperClusters"),
                 recoPhotonSrc             = cms.InputTag("reducedEgamma", "reducedGedPhotonCores"),
                 TrackLabel                = cms.InputTag("generalTracks"),
                 gsfElectronLabel          = cms.InputTag("gsfElectrons"),
                 PFAllCandidates           = cms.InputTag("particleFlow"),
                 #ak4PFJetsCHS                 = cms.InputTag("updatedJets"),
                 ak4PFJetsCHSSrc           = cms.InputTag("slimmedJets"),
                 ak4PFJetsPUPPISrc           = cms.InputTag("slimmedJetsPuppi"),
                 ak4PFJetsCHSGenJetLabel      = cms.InputTag("slimmedGenJets"),
                 ak8GenJetLabel                 = cms.InputTag("slimmedGenJetsAK8"),
                 ak8JetsPUPPISrc                 = cms.InputTag("slimmedJetsAK8"),
                 tauSrc                    = cms.InputTag("slimmedTaus"),

                 packedPFCands             = cms.InputTag("packedPFCandidates"),
                 elePFClusEcalIsoProducer  = cms.InputTag("electronEcalPFClusterIsolationProducer"),
                 elePFClusHcalIsoProducer  = cms.InputTag("electronHcalPFClusterIsolationProducer"),
                 BadChargedCandidateFilter = cms.InputTag("BadChargedCandidateFilter"),
                 BadPFMuonFilter           = cms.InputTag("BadPFMuonFilter")
 )
