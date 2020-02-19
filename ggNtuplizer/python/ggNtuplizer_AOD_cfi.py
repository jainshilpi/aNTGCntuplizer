import FWCore.ParameterSet.Config as cms
#from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector
#from CMGTools.External.pujetidproducer_cfi import stdalgos, chsalgos

ggNtuplizer = cms.EDAnalyzer("ggNtuplizer",
                             doGenParticles       = cms.bool(True),
                            #doGenJets           = cms.bool(True),
                             runOnParticleGun     = cms.bool(False),
                             runOnSherpa          = cms.bool(True),
                             dumpPFPhotons        = cms.bool(True),
                             dumpJets             = cms.bool(False),
                             dumpAK8Jets          = cms.bool(False),
                             dumpTaus             = cms.bool(False),
                             dumpPDFSystWeight    = cms.bool(False),
                             dumpHFElectrons      = cms.bool(True),
                             doRecHits            = cms.bool(True),
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
                             pfMETLabel           = cms.InputTag("patMETs"),
                             electronSrc          = cms.InputTag("selectedPatElectrons"),
                             #calibelectronSrc     = cms.InputTag("calibratedPatElectrons"),
                             calibelectronSrc     = cms.InputTag("selectedPatElectrons"),
                             photonSrc            = cms.InputTag("selectedPatPhotons"),
                             photonOOTSrc            = cms.InputTag("selectedPatOOTPhotons"),
                             #calibphotonSrc       = cms.InputTag("calibratedPatPhotons"),
                             calibphotonSrc       = cms.InputTag("selectedPatPhotons"),
                             muonSrc              = cms.InputTag("selectedPatMuons"),
                             gsfTrackSrc          = cms.InputTag("electronGsfTracks"),
                             ebReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
                             eeReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
                             esReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsES"),
                             ecalSCcollection       = cms.InputTag("particleFlowSuperClusterECAL", "particleFlowSuperClusterECALBarrel"),
                             #ecalSCcollectionEB       = cms.InputTag("particleFlowSuperClusterECAL", "particleFlowSuperClusterECALBarrel"),
                             #ecalSCcollectionEE       = cms.InputTag("particleFlowSuperClusterECAL", "particleFlowSuperClusterECALEndcapWithPreshower"),
                             ecalSCOOTcollection       = cms.InputTag("", "reducedOOTSuperClusters"),
                             
                             hbheRecHitCollection      = cms.InputTag("reducedHcalRecHits","hbhereco"),
                             hoRecHitCollection      = cms.InputTag("reducedHcalRecHits","horeco"),
                             hfRecHitCollection      = cms.InputTag("reducedHcalRecHits","hfreco"),
                             cscSegmentsCollection   = cms.InputTag("cscSegments"),
                             
                             #ecalSCOOTcollectionEB       = cms.InputTag("particleFlowSuperClusterOOTECAL", "particleFlowSuperClusterOOTECALBarrel"),
                             #ecalSCOOTcollectionEE       = cms.InputTag("particleFlowSuperClusterOOTECAL", "particleFlowSuperClusterOOTECALEndcapWithPreshower"),
                             recoPhotonSrc             = cms.InputTag("gedPhotons"),
                             TrackLabel                = cms.InputTag("generalTracks"),
                             gsfElectronLabel          = cms.InputTag("gsfElectrons"),
                             PFAllCandidates           = cms.InputTag("particleFlow"),
                             #ak4PFJetsCHS                 = cms.InputTag("updatedJets"),
                             #ak4PFJetsCHSSrc           = cms.InputTag("selectedPatJetsAK4PFCHS"),
                             ak4PFJetsCHSSrc           = cms.InputTag("selectedPatJets"),
                             ak4PFJetsPUPPISrc           = cms.InputTag(""), #####----> not available in AOD
                             ak4PFJetsCHSGenJetLabel      = cms.InputTag("ak4GenJets"),
                             ak8GenJetLabel                 = cms.InputTag("ak8GenJets"),
                             ak8JetsPUPPISrc                 = cms.InputTag("slimmedJetsAK8"),  #####----> not available in AOD 
                             tauSrc                    = cms.InputTag("selectedPatTaus"),
                             
                             packedPFCands             = cms.InputTag("packedPFCandidates"),

                             offlineBeamSpot = cms.InputTag("offlineBeamSpot"),
                             #BadChargedCandidateFilter = cms.InputTag("BadChargedCandidateFilter"), ####---> not available in AOD
                             #BadPFMuonFilter           = cms.InputTag("BadPFMuonFilter"),     #### ----> not available in AOD

                             #new reso: from here: https://github.com/cms-analysis/flashgg/blob/dev_legacy_runII/MicroAOD/python/flashggDiPhotons_cfi.py
                             sigma1Pix               = cms.double(0.0125255),
                             sigma1Tib               = cms.double(0.716301),
                             sigma1Tob               = cms.double(3.17615),
                             sigma1PixFwd            = cms.double(0.0581667),
                             sigma1Tid               = cms.double(0.38521),
                             sigma1Tec               = cms.double(1.67937),
                             sigma2Pix               = cms.double(0.0298574),
                             sigma2Tib               = cms.double(0.414393),
                             sigma2Tob               = cms.double(1.06805),
                             sigma2PixFwd            = cms.double(0.180419),
                             sigma2Tid               = cms.double(0.494722),
                             sigma2Tec               = cms.double(1.21941),
                             singlelegsigma1Pix      = cms.double(0.0178107),
                             singlelegsigma1Tib      = cms.double(1.3188),
                             singlelegsigma1Tob      = cms.double(2.23662),
                             singlelegsigma1PixFwd   = cms.double(0.152157),
                             singlelegsigma1Tid      = cms.double(0.702755),
                             singlelegsigma1Tec      = cms.double(2.46599),
                             singlelegsigma2Pix      = cms.double(0.0935307),
                             singlelegsigma2Tib      = cms.double(0.756568),
                             singlelegsigma2Tob      = cms.double(0.62143),
                             singlelegsigma2PixFwd   = cms.double(0.577081),
                             singlelegsigma2Tid      = cms.double(0.892751),
                             singlelegsigma2Tec      = cms.double(1.56638),
                            
                             phoLooseIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Fall17-94X-V2-loose"),
                             phoMediumIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Fall17-94X-V2-medium"),
                             phoTightIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Fall17-94X-V2-tight"),
                             #phoMVAValuesMap = cms.InputTag("photonMVAValueMapProducer:mvaPhoID-RunIIFall17-v2-wp90"),
                             phoMVAValuesMap = cms.InputTag("egmPhotonIDs:mvaPhoID-RunIIFall17-v2-wp90"),
                             phoChargedIsolation       = cms.InputTag("photonIDValueMapProducer:phoChargedIsolation"),
                             phoNeutralHadronIsolation = cms.InputTag("photonIDValueMapProducer:phoNeutralHadronIsolation"),
                             phoPhotonIsolation        = cms.InputTag("photonIDValueMapProducer:phoPhotonIsolation"),
                             phoWorstChargedIsolation  = cms.InputTag("photonIDValueMapProducer:phoWorstChargedIsolation"),
                             eleVetoIdMap    = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-veto"),
                             eleLooseIdMap   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-loose"),
                             eleMediumIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-medium"),
                             eleTightIdMap   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-tight"),
                             eleHEEPIdMap    = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV70"),
                             #eleMVAValuesMap = cms.InputTag("electronMVAValueMapProducer:mvaEleID-Fall17-noIso-V2-wp90")
                             eleMVAValuesMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V2-wp90")

)

