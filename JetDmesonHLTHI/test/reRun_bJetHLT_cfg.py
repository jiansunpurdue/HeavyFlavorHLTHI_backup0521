import FWCore.ParameterSet.Config as cms

doRegionalClusters = True

process = cms.Process('HLT')

#process.SimpleMemoryCheck=cms.Service("SimpleMemoryCheck",
#oncePerEventMode=cms.untracked.bool(False))

process.Timing=cms.Service("Timing",
                           useJobReport = cms.untracked.bool(True)
                           )

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('HLTrigger.Configuration.HLT_User2_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
                            #secondaryFileNames = cms.untracked.vstring('file:/tmp/mnguyen/step2_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_PU.root'),
                            fileNames = cms.untracked.vstring(
        #'/store/user/mnguyen/PyquenUnquenched_Dijet_pthat80_740pre6_GEN-SIM/PyquenUnquenched_Dijet_pthat80_740pre8_MCHI2_74_V0_DIGI-RAW/ee815b27030c232e2e0a7be48a50a463/step2_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_PU_9_1_jsB.root'
#        '/store/user/mnguyen/PyquenUnquenched_BjetLO_pthat80_740pre6_GEN-SIM/PyquenUnquenched_BjetLO_pthat80_740pre8_MCHI2_74_V0_DIGI-RAW/ee815b27030c232e2e0a7be48a50a463/step2_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_PU_9_1_GT1.root'
#		'file:./step3_RAW2DIGI_L1Reco_RECO.root'
		"file:/home/sun229/store/PyquenUnquenched_D0_pthat80_740pre8_MCHI1_74_V4_GEN-SIM_MSEL1_0321/PyquenUnquenched_D0_NcollFilt_pthat80_740pre8_MCHI2_74_pre8_rerunHLTreco_MSEL1_0321/894cc73786c941885631deb58e5b7664/step3_RAW2DIGI_L1Reco_RECO_100_1_Mz7.root"
        ),
                            #eventsToProcess = cms.untracked.VEventRange("1:23-1:23")
                            #skipEvents = cms.untracked.uint32(31),
                            )

process.options = cms.untracked.PSet( 
    wantSummary = cms.untracked.bool(True) 
)


# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.381.2.28 $'),
    annotation = cms.untracked.string('bjetHLT nevts:10'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition
if doRegionalClusters:
    myFileName = 'rerunHLT_RAW2DIGI_L1Reco.root'
else:
    myFileName = 'rerunHLT_RAW2DIGI_L1Reco_noReg.root'


#myFileName = 'rerunHLT_RAW2DIGI_L1Reco_regVtx_glTrk.root'

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
                                        splitLevel = cms.untracked.int32(0),
                                        eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
                                        outputCommands = process.RAWSIMEventContent.outputCommands,
                                        #outputCommands = cms.untracked.vstring(['keep *']),       
                                        fileName = cms.untracked.string(myFileName),
                                        dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RAW')
        )
)

# Additional output definition
process.RAWSIMoutput.outputCommands.extend(['drop *_g4SimHits_*_*', 'drop *_*_*_SIM', 'drop *_*_*_L1Reco', 'drop *_*_*_RECO'])
#process.RAWSIMoutput.outputCommands.extend([
#process.RAWSIMoutput.outputCommands = ([
#        #'keep *_*_*_BJETHLT'
#        'keep *_hltIter1PFJetPixelSeeds_*_BJETHLT',
#        'keep *_hltIter1PFJetCkfTrackCandidates_*_BJETHLT',
#        'keep *_hltIter1PFJetCtfWithMaterialTracks_*_BJETHLT',
#        'keep *_hltIter1Selector_*_BJETHLT',
#        'keep *_hltIter1Merged_*_BJETHLT',
#        'keep *_hltIter2PFJetCtfWithMaterialTracks_*_BJETHLT',
#        'keep *_hltIter2Merged_*_BJETHLT',
#        'keep *_hltIter3Merged_*_BJETHLT',
#        'keep *_hltIter4Merged_*_BJETHLT',
#        'keep *_hltHISelectedVertex_*_BJETHLT',
#        'keep *_hltOnlineBeamSpot_*_BJETHLT',
#        'keep *_hltSiPixelClusters_*_BJETHLT',
#        'keep *_hltSiStripClusters_*_BJETHLT',
#        'keep *_hltSiStripRawToClustersFacility_*_BJETHLT',
#        #'keep *_hltHISiPixelClusters_*_BJETHLT',
#        #'keep *_hltHISiStripClusters_*_BJETHLT',
#        #'keep *_hltHISiStripRawToClustersFacility_*_BJETHLT',
#        'keep *_TriggerResults_*_BJETHLT'
#    ])
#
# Other statements
from HLTrigger.Configuration.CustomConfigs import ProcessName
process = ProcessName(process)

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc_hi', '')


from CondCore.DBCommon.CondDBSetup_cfi import *
process.beamspot = cms.ESSource("PoolDBESSource",CondDBSetup,
                                toGet = cms.VPSet(cms.PSet( record = cms.string('BeamSpotObjectsRcd'),
                                                            tag= cms.string('RealisticHICollisions2011_STARTHI50_mc')
                                                            )),
                                connect =cms.string('frontier://FrontierProd/CMS_COND_31X_BEAMSPOT')
                                )
process.es_prefer_beamspot = cms.ESPrefer("PoolDBESSource","beamspot")

# Path and EndPath definitions
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)


#import of pp trigger paths
process.load("UserCode.OpenHF.customHLT")
#process.hltSiPixelClusters.maxNumberOfClusters = 2000000

#jet selection
process.hltSinglePuAK4CaloJet80.MaxEta = 2.3
process.hltSelector4JetsDiJet20L1FastJet.maxNumber = 3
process.hltSelectorJets20L1FastJet.etMin = 40.
process.hltFastPixelBLifetimeL3Associator.coneSize = 0.4

process.load("CommonTools.RecoAlgos.caloJetSelector_cfi")
process.eta2CaloJets = process.caloJetSelector.clone(
    src = cms.InputTag("hltPuAK4CaloJetsCorrectedIDPassed"),
    cut = cms.string( "abs(eta)<2.3" )
    )

process.hltSelector4JetsDiJet20L1FastJet.src = "eta2CaloJets"
process.hltAK4CaloJetsCorrected.src = cms.InputTag("hltPuAK4CaloJets")
process.hltSelectorJets20L1FastJet.src = 'hltSelector4JetsDiJet20L1FastJet'
#process.hltSelectorJets20L1FastJet.src = 'eta2CaloJets'
#process.hltSelectorJets20L1FastJet.src = 'hltPuAK4CaloJetsCorrectedIDPassed'


process.reduceJetMult = process.hltSelector4JetsDiJet20L1FastJet.clone(maxNumber = 3)

process.jets4bTagger = process.hltSelectorJets20L1FastJet.clone(
    etMin = 80., 
    #src = 'hltSelectorJets20L1FastJet'
    src = 'reduceJetMult'
)

#CSV tagger

# local reco
process.load("RecoTracker.TransientTrackingRecHit.TransientTrackingRecHitBuilder_cfi")
process.load("RecoLocalTracker.SiPixelRecHits.PixelCPEGeneric_cfi")
process.load("RecoTracker.TkSeedingLayers.TTRHBuilderWithoutAngle4PixelTriplets_cfi")

#test HI ZS
#process.load("EventFilter.SiStripRawToDigi.SiStripDigis_cfi")
#process.load("RecoLocalTracker.SiStripZeroSuppression.SiStripZeroSuppression_cfi")
#process.hltSiStripRawToClustersFacility.Algorithms = process.siStripZeroSuppression.Algorithms

# TEST:  Make it look like global reco
process.hltSiPixelDigisRegForBTag.Regions = cms.PSet()
#process.hltESPPixelCPEGeneric.useLAWidthFromDB = True
#process.hltSiPixelDigisRegForBTag.ErrorList = cms.vint32(29)
#process.hltSiPixelDigisRegForBTag.IncludeErrors = cms.bool(True)
#process.hltSiPixelDigisRegForBTag.UserErrorList = cms.vint32(40)

process.HLTDoLocalPixelSequenceRegForBTag = cms.Sequence( process.hltSiPixelDigisRegForBTag + process.hltSiPixelClustersRegForBTag + process.hltSiPixelClustersRegForBTagCache + process.hltSiPixelRecHitsRegForBTag + process.hltPixelLayerPairsRegForBTag + process.hltPixelLayerTripletsRegForBTag )
    #pixels                                                                                                                                                                    

#TEST:  switch to global pixel only  
'''
process.load("RecoLocalTracker.SiPixelRecHits.SiPixelRecHits_cfi")
process.load("RecoPixelVertexing.PixelLowPtUtilities.siPixelClusterShapeCache_cfi")
process.load("EventFilter.SiPixelRawToDigi.SiPixelRawToDigi_cfi")
process.siPixelDigis.InputLabel = "rawDataCollector"
process.load("RecoLocalTracker.SiPixelClusterizer.SiPixelClusterizer_cfi")
process.HLTDoLocalPixelSequenceRegForBTag = cms.Sequence(process.siPixelDigis
                                                         + process.siPixelClusters
                                                         + process.siPixelRecHits
                                                         + process.siPixelClusterShapeCache)
'''
process.HLTDoLocalStripSequenceRegForBTag = cms.Sequence( process.hltSiStripExcludedFEDListProducer #+ process.siStripDigis + process.siStripZeroSuppression
 + process.hltSiStripRawToClustersFacility + process.hltSiStripClustersRegForBTag )
#process.hltSiPixelDigisRegForBTag.Regions.inputs = ["hltSelectorJets20L1FastJet"]

#HI vertexing, taken from 2011 single track triggers
#process.HLTHINonRegionalPixelVertexSequence = cms.Sequence(process.hltHIPixelClusterVerticesForHITrackTrigger+process.hltHIPixel3ProtoTracks+process.hltHIPixelMedianVertex+process.hltHISelectedProtoTracks+process.hltHIPixelAdaptiveVertex+process.hltHIBestAdaptiveVertex+process.hltHISelectedVertex)
# not yet available!  Use offline version
process.load("RecoHI.HiTracking.HIPixelVertices_cff")
process.load("RecoHI.HiTracking.HIPixel3PrimTracks_cfi")
process.hiPixel3ProtoTracks.FilterPSet.beamSpot = cms.InputTag("hltOnlineBeamSpot")
process.HLTHINonRegionalPixelVertexSequence = cms.Sequence(process.hiPixelClusterVertex+process.PixelLayerTriplets+process.hiPixel3ProtoTracks+process.hiPixelMedianVertex+process.hiSelectedProtoTracks+process.hiPixelAdaptiveVertex+process.hiBestAdaptiveVertex+process.hiSelectedVertex+process.hiPixel3PrimTracks
)

if doRegionalClusters:
    process.PixelLayerTriplets.BPix.HitProducer = 'hltSiPixelRecHitsRegForBTag'
    process.PixelLayerTriplets.FPix.HitProducer = 'hltSiPixelRecHitsRegForBTag'
    process.hiPixel3ProtoTracks.FilterPSet.siPixelRecHits = 'hltSiPixelRecHitsRegForBTag'
    process.hiPixel3ProtoTracks.RegionFactoryPSet.RegionPSet.siPixelRecHits = 'hltSiPixelRecHitsRegForBTag'
    process.hiPixel3PrimTracks.FilterPSet.clusterShapeCacheSrc = 'hltSiPixelClustersRegForBTagCache'
    process.hiPixelClusterVertex.pixelRecHits = "hltSiPixelRecHitsRegForBTag"

    process.localTrackerReco = cms.Sequence(process.HLTDoLocalPixelSequenceRegForBTag + process.HLTDoLocalStripSequenceRegForBTag)
    

else:
#if True == True:
    process.PixelLayerTriplets.BPix.HitProducer = 'siPixelRecHits'
    process.PixelLayerTriplets.FPix.HitProducer = 'siPixelRecHits'
    process.hiPixel3ProtoTracks.FilterPSet.siPixelRecHits = 'siPixelRecHits'
    process.hiPixel3ProtoTracks.RegionFactoryPSet.RegionPSet.siPixelRecHits = 'siPixelRecHits'
    process.hiPixel3PrimTracks.FilterPSet.clusterShapeCacheSrc = 'siPixelClusterShapeCache'

    process.hltESPTTRHBuilderPixelOnly.PixelCPE = cms.string('PixelCPE')    
    process.HLTIterativeTrackingForBTagIteration1.remove(process.hltIter1MaskedMeasurementTrackerEventForBTag)
    process.HLTIterativeTrackingForBTagIteration2.remove(process.hltIter2MaskedMeasurementTrackerEventForBTag)
    process.hltPixelLayerPairsRegForBTag.FPix.HitProducer = 'siPixelRecHits'
    process.hltPixelLayerPairsRegForBTag.BPix.HitProducer = 'siPixelRecHits'
    process.hltPixelLayerPairsRegForBTag.BPix.TTRHBuilder = cms.string('WithTrackAngle')
    process.hltPixelLayerPairsRegForBTag.FPix.TTRHBuilder = cms.string('WithTrackAngle')
    process.hltPixelLayerTripletsRegForBTag.FPix.HitProducer = 'siPixelRecHits'
    process.hltPixelLayerTripletsRegForBTag.BPix.HitProducer = 'siPixelRecHits'
    process.hltPixelLayerTripletsRegForBTag.BPix.TTRHBuilder = cms.string('WithTrackAngle')
    process.hltPixelLayerTripletsRegForBTag.FPix.TTRHBuilder = cms.string('WithTrackAngle')
    process.hltIter0PFlowPixelSeedsFromPixelTracksForBTag.TTRHBuilder = cms.string('WithTrackAngle')
    process.hltIter1PixelLayerTripletsForBTag.BPix.TTRHBuilder = cms.string('WithTrackAngle')
    process.hltIter1PixelLayerTripletsForBTag.FPix.TTRHBuilder = cms.string('WithTrackAngle')
    process.hltIter1PixelLayerTripletsForBTag.BPix.HitProducer = 'siPixelRecHits'
    process.hltIter1PixelLayerTripletsForBTag.FPix.HitProducer = 'siPixelRecHits'
    process.hltIter1ClustersRefRemovalForBTag.pixelClusters = 'siPixelClusters'
    process.hltIter1ClustersRefRemovalForBTag.stripClusters = 'siStripClusters'
    process.hltIter0PFlowCkfTrackCandidatesForBTag.MeasurementTrackerEvent = 'MeasurementTrackerEvent'
    del process.hltIter1PixelLayerTripletsForBTag.BPix.skipClusters
    del process.hltIter1PixelLayerTripletsForBTag.FPix.skipClusters
    del process.hltIter2PixelLayerPairsForBTag.BPix.skipClusters
    del process.hltIter2PixelLayerPairsForBTag.FPix.skipClusters
    process.hltIter1PFlowPixelSeedsForBTag.RegionFactoryPSet.RegionPSet.measurementTrackerName = "MeasurementTrackerEvent"
    process.hltIter1PFlowCkfTrackCandidatesForBTag.MeasurementTrackerEvent = 'MeasurementTrackerEvent'
    process.hltIter2ClustersRefRemovalForBTag.pixelClusters = 'siPixelClusters'
    process.hltIter2ClustersRefRemovalForBTag.stripClusters = 'siStripClusters'
    process.hltIter2PixelLayerPairsForBTag.BPix.TTRHBuilder = cms.string('WithTrackAngle')
    process.hltIter2PixelLayerPairsForBTag.FPix.TTRHBuilder = cms.string('WithTrackAngle')
    process.hltIter2PixelLayerPairsForBTag.BPix.HitProducer = 'siPixelRecHits'
    process.hltIter2PixelLayerPairsForBTag.FPix.HitProducer = 'siPixelRecHits'
    process.hltIter2PFlowPixelSeedsForBTag.RegionFactoryPSet.RegionPSet.measurementTrackerName = "MeasurementTrackerEvent"
    process.hltIter2PFlowCkfTrackCandidatesForBTag.MeasurementTrackerEvent = 'MeasurementTrackerEvent'
    
    #pixels
    process.load("RecoLocalTracker.SiPixelRecHits.SiPixelRecHits_cfi")
    process.load("RecoPixelVertexing.PixelLowPtUtilities.siPixelClusterShapeCache_cfi")
    process.load("EventFilter.SiPixelRawToDigi.SiPixelRawToDigi_cfi")
    process.siPixelDigis.InputLabel = 'rawDataCollector'
    process.load("RecoLocalTracker.SiPixelClusterizer.SiPixelClusterizer_cfi")
    #strips
    process.load("EventFilter.SiStripRawToDigi.SiStripDigis_cfi")
    process.load("RecoLocalTracker.SiStripClusterizer.SiStripClusterizer_cfi")
    process.load("RecoLocalTracker.SiStripZeroSuppression.SiStripZeroSuppression_cfi")
    process.load("RecoTracker.MeasurementDet.MeasurementTrackerEventProducer_cfi")
    process.load('RecoTracker.MeasurementDet.MeasurementTrackerESProducer_cff')

    process.localTrackerReco = cms.Sequence(
        process.siPixelDigis
        + process.siPixelClusters 
        + process.siPixelRecHits 
        + process.siPixelClusterShapeCache
        + process.siStripDigis
        + process.siStripZeroSuppression
        + process.siStripClusters
        + process.MeasurementTrackerEvent
        )


process.HLTFastPrimaryVertexSequence = cms.Sequence( process.eta2CaloJets + process.hltSelector4JetsDiJet20L1FastJet + 
                                                     process.hltSelectorJets20L1FastJet
                                                     + process.reduceJetMult
                                                     + process.jets4bTagger
                                                     + process.localTrackerReco                                                    
                                                     + process.HLTHINonRegionalPixelVertexSequence
                                                     )






#L3 sequences
process.HLTBtagCSVSequenceL3 = cms.Sequence( process.HLTIterativeTrackingForBTagIter02 + process.hltVerticesL3 + process.hltFastPixelBLifetimeL3Associator + process.hltFastPixelBLifetimeL3TagInfos + process.hltL3SecondaryVertexTagInfos + process.hltL3CombinedSecondaryVertexBJetTags )
process.HLTD0TagSequenceL3 = cms.Sequence( process.HLTIterativeTrackingForBTagIter02 + process.hltVerticesL3 + process.hltFastPixelBLifetimeL3Associator )



process.hltIter0PFlowPixelSeedsFromPixelTracksForBTag.useEventsWithNoVertex = False
process.hltIter0PFlowPixelSeedsFromPixelTracksForBTag.InputCollection = 'hiPixel3PrimTracks'
process.hltIter1PFlowPixelSeedsForBTag.RegionFactoryPSet.RegionPSet.input = 'hltSelectorJets20L1FastJet'
process.hltIter2PFlowPixelSeedsForBTag.RegionFactoryPSet.RegionPSet.input = 'hltSelectorJets20L1FastJet'
process.hltFastPixelBLifetimeL3Associator.jets = 'jets4bTagger'
process.hltBLifetimeL3FilterCSV.Jets = 'jets4bTagger'
process.hltIter0PFlowTrackSelectionHighPurityForBTag.min_nhits = 10
process.hltIter1PFlowTrackSelectionHighPurityLooseForBTag.min_nhits = 10
process.hltIter1PFlowTrackSelectionHighPurityTightForBTag.min_nhits = 10
process.hltIter2PFlowTrackSelectionHighPurityForBTag.min_nhits = 10
process.hltIter0PFlowTrackSelectionHighPurityForBTag.max_relpterr = 0.1
process.hltIter1PFlowTrackSelectionHighPurityLooseForBTag.max_relpterr = 0.1
process.hltIter1PFlowTrackSelectionHighPurityTightForBTag.max_relpterr = 0.1
process.hltIter2PFlowTrackSelectionHighPurityForBTag.max_relpterr = 0.1
process.hltIter0PFlowTrackSelectionHighPurityForBTag.chi2n_par = 0.3
process.hltIter1PFlowTrackSelectionHighPurityLooseForBTag.chi2n_par = 0.3
process.hltIter1PFlowTrackSelectionHighPurityTightForBTag.chi2n_par = 0.2
process.hltIter2PFlowTrackSelectionHighPurityForBTag.chi2n_par = 0.3


process.HLT_PuAK4CaloJet80PlusBJetCSV_v1 =  cms.Path(process.HLTBeginSequence*process.hltL1sL1SingleJet16BptxAND*process.hltPrePuAK4CaloJet80 * process.hltPrePuAK4CaloJet80 * process.HLTPuAK4CaloJetsSequence * process.hltSinglePuAK4CaloJet80
                                                     * process.HLTFastPrimaryVertexSequence
                                                     * process.HLTBtagCSVSequenceL3
                                                     * process.hltBLifetimeL3FilterCSV
                                                     * process.HLTEndSequence 
                                                     )

# D meson
process.dMesonTagger = cms.Sequence(process.hltDisplacedVtxProducerD0+process.hltVertexFilterD0)
#play with cuts

process.hltDisplacedVtxProducerD0.MinPtPair = cms.double( 10.0 )
process.hltVertexFilterD0.MinLxyzSignificance = 2.5
process.hltVertexFilterD0.MinCosinePointingAngle = 0.95
#process.hltVertexFilterD0.MaxNormalisedChi2 = 3.0

process.HLT_PuAK4CaloJet80PlusD0_v1 =  cms.Path(process.HLTBeginSequence*process.hltL1sL1SingleJet16BptxAND*process.hltPrePuAK4CaloJet80 * process.hltPrePuAK4CaloJet80 * process.HLTPuAK4CaloJetsSequence * process.hltSinglePuAK4CaloJet80
                                                * process.HLTFastPrimaryVertexSequence
                                                * process.HLTD0TagSequenceL3
                                                * process.dMesonTagger
                                                * process.HLTEndSequence
                                                )


# run the whole menu again
#process.HLTSchedule.extend(cms.Schedule(process.HLT_HIJet80PlusBJetCSV_v7))
#only run the triggers of interest
process.HLTSchedule = cms.Schedule(process.HLT_PuAK4CaloJet80_v1,process.HLT_PuAK4CaloJet80PlusBJetCSV_v1,process.HLT_PuAK4CaloJet80PlusD0_v1 )
#process.HLTSchedule = cms.Schedule(process.HLT_PuAK4CaloJet80PlusBJetCSV_v1 )
#process.HLTSchedule = cms.Schedule(process.HLT_PuAK4CaloJet80_v1 )



process.schedule = cms.Schedule(process.HLTSchedule)
process.schedule.extend([process.endjob_step,process.RAWSIMoutput_step])



# customisation of the process.

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.postLS1Customs
from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1 

#call to customisation function customisePostLS1 imported from SLHCUpgradeSimulations.Configuration.postLS1Customs
process = customisePostLS1(process)

# Automatic addition of the customisation function from L1Trigger.Configuration.L1Trigger_custom
from L1Trigger.Configuration.customise_overwriteL1Menu import L1Menu_CollisionsHeavyIons2015_v0

#call to customisation function customiseL1Menu_HI imported from L1Trigger.Configuration.L1Trigger_custom                                                                       
process = L1Menu_CollisionsHeavyIons2015_v0(process)

# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforMC 

#call to customisation function customizeHLTforMC imported from HLTrigger.Configuration.customizeHLTforMC
process = customizeHLTforMC(process)

process.simCaloStage1Digis.FirmwareVersion = 1
# End of customisation functions


# fix JECs in HLT.  Will be fixed in upcoming HLT menu
process.hltESPAK4CaloCorrectionNoL1 = process.hltESPAK4CaloCorrection.clone(
     correctors = cms.vstring('hltESPAK4CaloRelativeCorrectionESProducer',
'hltESPAK4CaloAbsoluteCorrectionESProducer')
)
process.hltPuAK4CaloJetsCorrectedIDPassed.correctors = cms.vstring('hltESPAK4CaloCorrectionNoL1' )

# not needed
process.HLTPuAK4CaloJetsSequence.remove(process.hltFixedGridRhoFastjetAllCalo)
process.HLTPuAK4CaloJetsSequence.remove(process.hltPuAK4CaloJetsCorrected)



# vertex swap to hiSelectedVertex
oldFastPVTag=cms.InputTag("hltFastPVPixelVertices")
newFastPVTag=cms.InputTag("hiSelectedVertex")
oldFullPVTag=cms.InputTag("hltVerticesL3")
newFullPVTag=cms.InputTag("hiSelectedVertex")
oldFullPVBSTag=cms.InputTag("hltVerticesL3","WithBS")
oldBSTag=cms.InputTag("offlineBeamSpot")
newBSTag=cms.InputTag("hltOnlineBeamSpot")
#test
#oldPixTag=cms.InputTag("hltSiPixelClustersRegForBTag")
#newPixTag=cms.InputTag("siPixelClusters")
#process.hltSiStripClustersRegForBTag.pixelClusterProducer = 'siPixelClusters'

from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag
for s in process.paths_().keys():
    massSearchReplaceAnyInputTag(getattr(process,s),oldFastPVTag,newFastPVTag)
    massSearchReplaceAnyInputTag(getattr(process,s),oldFullPVTag,newFullPVTag)
    massSearchReplaceAnyInputTag(getattr(process,s),oldFullPVBSTag,newFullPVTag)
    massSearchReplaceAnyInputTag(getattr(process,s),oldBSTag,newBSTag)
    #massSearchReplaceAnyInputTag(getattr(process,s),oldPixTag,newPixTag) #test



process.hltIter1PFlowPixelSeedsForBTag.OrderedHitsFactoryPSet.maxElement = 20000
process.hltIter2PFlowPixelSeedsForBTag.OrderedHitsFactoryPSet.maxElement = 20000
#process.hltSiPixelDigisRegForBTag.Regions.maxZ = cms.vdouble(50.)
#process.hltSiPixelDigisRegForBTag.Regions.deltaPhi = cms.vdouble(0.6)

