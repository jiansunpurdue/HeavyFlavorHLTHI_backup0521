import FWCore.ParameterSet.Config as cms

#customize the HLT menu w/ b jet stuff

#jet selctor
from HLTrigger.Configuration.HLT_FULL_cff import hltSelector4JetsDiJet20L1FastJet
hltSelector4JetsDiJet20L1FastJet.src = cms.InputTag( "hltHICaloJetCorrected")
                                                 
from HLTrigger.Configuration.HLT_FULL_cff import hltSelectorJets20L1FastJet

#on the  wishlist
#from HLTrigger.Configuration.HLT_FULL_cff import HLTPixelSeedingForHITrackTrigger

# b-taggers
from HLTrigger.Configuration.HLT_FULL_cff import hltFastPixelBLifetimeL3Associator
from HLTrigger.Configuration.HLT_FULL_cff import hltFastPixelBLifetimeL3TagInfos
from HLTrigger.Configuration.HLT_FULL_cff import hltL3SecondaryVertexTagInfos
from HLTrigger.Configuration.HLT_FULL_cff import hltL3CombinedSecondaryVertexBJetTags
from HLTrigger.Configuration.HLT_FULL_cff import hltBLifetimeL3FilterCSV

# add real iterative tracking!



#step 0 
from HLTrigger.Configuration.HLT_FULL_cff import hltIter0PFlowPixelSeedsFromPixelTracksForBTag
from HLTrigger.Configuration.HLT_FULL_cff import hltIter0PFlowCkfTrackCandidatesForBTag
from HLTrigger.Configuration.HLT_FULL_cff import hltIter0PFlowCtfWithMaterialTracksForBTag
from HLTrigger.Configuration.HLT_FULL_cff import hltIter0PFlowTrackSelectionHighPurityForBTag


HLTIterativeTrackingForBTagIteration0 = cms.Sequence( hltIter0PFlowPixelSeedsFromPixelTracksForBTag + hltIter0PFlowCkfTrackCandidatesForBTag + hltIter0PFlowCtfWithMaterialTracksForBTag + hltIter0PFlowTrackSelectionHighPurityForBTag 
)

from HLTrigger.Configuration.HLT_FULL_cff import hltIter1ClustersRefRemovalForBTag
from HLTrigger.Configuration.HLT_FULL_cff import hltIter1MaskedMeasurementTrackerEventForBTag
from HLTrigger.Configuration.HLT_FULL_cff import hltIter1PixelLayerTripletsForBTag
from HLTrigger.Configuration.HLT_FULL_cff import hltIter1PFlowPixelSeedsForBTag
from HLTrigger.Configuration.HLT_FULL_cff import hltIter1PFlowCkfTrackCandidatesForBTag
from HLTrigger.Configuration.HLT_FULL_cff import hltIter1PFlowCtfWithMaterialTracksForBTag
from HLTrigger.Configuration.HLT_FULL_cff import hltIter1PFlowTrackSelectionHighPurityLooseForBTag
from HLTrigger.Configuration.HLT_FULL_cff import hltIter1PFlowTrackSelectionHighPurityTightForBTag
from HLTrigger.Configuration.HLT_FULL_cff import hltIter1PFlowTrackSelectionHighPurityForBTag

HLTIterativeTrackingForBTagIteration1 = cms.Sequence( hltIter1ClustersRefRemovalForBTag 
+ hltIter1MaskedMeasurementTrackerEventForBTag 
+ hltIter1PixelLayerTripletsForBTag + hltIter1PFlowPixelSeedsForBTag + hltIter1PFlowCkfTrackCandidatesForBTag + hltIter1PFlowCtfWithMaterialTracksForBTag 
+ hltIter1PFlowTrackSelectionHighPurityLooseForBTag + hltIter1PFlowTrackSelectionHighPurityTightForBTag
 + hltIter1PFlowTrackSelectionHighPurityForBTag )

from HLTrigger.Configuration.HLT_FULL_cff import hltIter2ClustersRefRemovalForBTag
from HLTrigger.Configuration.HLT_FULL_cff import hltIter2MaskedMeasurementTrackerEventForBTag
from HLTrigger.Configuration.HLT_FULL_cff import hltIter2PixelLayerPairsForBTag
from HLTrigger.Configuration.HLT_FULL_cff import hltIter2PFlowPixelSeedsForBTag
from HLTrigger.Configuration.HLT_FULL_cff import hltIter2PFlowCkfTrackCandidatesForBTag
from HLTrigger.Configuration.HLT_FULL_cff import hltIter2PFlowCtfWithMaterialTracksForBTag
from HLTrigger.Configuration.HLT_FULL_cff import hltIter2PFlowTrackSelectionHighPurityForBTag

HLTIterativeTrackingForBTagIteration2 = cms.Sequence( hltIter2ClustersRefRemovalForBTag 
+ hltIter2MaskedMeasurementTrackerEventForBTag 
+ hltIter2PixelLayerPairsForBTag + hltIter2PFlowPixelSeedsForBTag + hltIter2PFlowCkfTrackCandidatesForBTag + hltIter2PFlowCtfWithMaterialTracksForBTag + hltIter2PFlowTrackSelectionHighPurityForBTag )



from HLTrigger.Configuration.HLT_FULL_cff import hltIter1MergedForBTag
from HLTrigger.Configuration.HLT_FULL_cff import hltIter2MergedForBTag



HLTIterativeTrackingForBTagIter02 = cms.Sequence( HLTIterativeTrackingForBTagIteration0 + HLTIterativeTrackingForBTagIteration1 + hltIter1MergedForBTag + HLTIterativeTrackingForBTagIteration2 + hltIter2MergedForBTag 
)

from HLTrigger.Configuration.HLT_FULL_cff import hltSiStripClustersRegForBTag
from HLTrigger.Configuration.HLT_FULL_cff import hltSiPixelDigisRegForBTag
from HLTrigger.Configuration.HLT_FULL_cff import hltSiPixelClustersRegForBTag
hltSiPixelClustersRegForBTag.maxNumberOfClusters = 999999
from HLTrigger.Configuration.HLT_FULL_cff import hltSiPixelClustersRegForBTagCache
from HLTrigger.Configuration.HLT_FULL_cff import hltSiPixelRecHitsRegForBTag
from HLTrigger.Configuration.HLT_FULL_cff import hltPixelLayerPairsRegForBTag
from HLTrigger.Configuration.HLT_FULL_cff import hltPixelLayerTripletsRegForBTag 

from HLTrigger.Configuration.HLT_FULL_cff import hltVerticesL3 

hltDisplacedVtxProducerD0 = cms.EDProducer( "HLTDisplacedDitkVtxProducer",
    Src = cms.InputTag( "hltFastPixelBLifetimeL3Associator" ),
    MinTrackPt = cms.double( 2.0 ),
    MaxEta = cms.double( 2.3 ),
    MaxInvMass = cms.double( 2.10 ),
    MinPtPair = cms.double( 10.0 ),
    MinInvMass = cms.double( 1.55 )
)

hltVertexFilterD0 = cms.EDFilter( "HLTDisplacedDmesonVtxFilter",
    saveTags = cms.bool( True ),
    PrimaryVerticesTag   = cms.InputTag("hiSelectedVertex"),
    BeamSpotTag = cms.InputTag( "hltOnlineBeamSpot" ),
	decaylength3D = cms.bool(True),
	MinLxyzSignificance = cms.double( 2.5 ),
	MinLxySignificance = cms.double( 2.0 ),
    MaxNormalisedChi2 = cms.double( 8.0 ),
    MinVtxProbability = cms.double( 0.02 ),
    MinCosinePointingAngle = cms.double( 0.95 ),
    DisplacedVertexTag = cms.InputTag( "hltDisplacedVtxProducerD0" )
)
