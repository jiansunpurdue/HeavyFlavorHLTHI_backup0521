// -*- C++ -*-
//
// Package:    JetDmesonHLTHI
// Class:      HLTDisplacedDmesonVtxFilter
//
// D meson secondary vertex filter for HI D meson triggers 
// changed from "HLTDisplacedmumuFilter"
//
// Original Author:  Jian Sun
//         Created:  Wed Nov. 5 2014
// $Id$
//
//

#include <iostream>

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerRefsCollections.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/Math/interface/Error.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/GlobalError.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "DataFormats/Common/interface/Ref.h"

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "HLTDisplacedDmesonVtxFilter.h"
#include "TMath.h"
#include <TVector3.h>

using namespace std;
using namespace reco;
using namespace edm;

//
// constructors and destructor
//
HLTDisplacedDmesonVtxFilter::HLTDisplacedDmesonVtxFilter(const edm::ParameterSet& iConfig) :
	HLTFilter(iConfig),
	decaylength3D_(iConfig.getParameter<bool>("decaylength3D")),
	minLxySignificance_ (iConfig.getParameter<double>("MinLxySignificance")),
	minLxyzSignificance_ (iConfig.getParameter<double>("MinLxyzSignificance")),
	maxNormalisedChi2_ (iConfig.getParameter<double>("MaxNormalisedChi2")), 
	minVtxProbability_ (iConfig.getParameter<double>("MinVtxProbability")),
	minCosinePointingAngle_ (iConfig.getParameter<double>("MinCosinePointingAngle")),
	DisplacedVertexTag_(iConfig.getParameter<edm::InputTag>("DisplacedVertexTag")),
	beamSpotTag_ (iConfig.getParameter<edm::InputTag> ("BeamSpotTag")),
	pvTag_ (iConfig.getParameter<edm::InputTag> ("PrimaryVerticesTag"))

{
}


HLTDisplacedDmesonVtxFilter::~HLTDisplacedDmesonVtxFilter()
{

}


void HLTDisplacedDmesonVtxFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	edm::ParameterSetDescription desc;
	makeHLTFilterDescription(desc);
	desc.add<bool>("decaylength3D",true);
	desc.add<double>("MinLxySignificance",2.0);
	desc.add<double>("MinLxyzSignificance",2.0);
	desc.add<double>("MaxNormalisedChi2",10.0);
	desc.add<double>("MinVtxProbability",0.0);
	desc.add<double>("MinCosinePointingAngle",0.9);
	desc.add<edm::InputTag>("DisplacedVertexTag",edm::InputTag("hltDisplacedVtxProducerD0"));
	desc.add<edm::InputTag>("BeamSpotTag",edm::InputTag("hltOnlineBeamSpot"));
	desc.add<edm::InputTag>("PrimaryVerticesTag",edm::InputTag("hltHISelectedVertex"));
	descriptions.add("hltDisplacedmumuFilter", desc);
}


// ------------ method called once each job just before starting event loop  ------------
void HLTDisplacedDmesonVtxFilter::beginJob()
{

}

// ------------ method called once each job just after ending the event loop  ------------
void HLTDisplacedDmesonVtxFilter::endJob() 
{

}

// ------------ method called on each new Event  ------------
bool HLTDisplacedDmesonVtxFilter::hltFilter(edm::Event& iEvent, const edm::EventSetup& iSetup, trigger::TriggerFilterObjectWithRefs & filterproduct) const
{

	// All HLT filters must create and fill an HLT filter object,
	// recording any reconstructed physics objects satisfying (or not)
	// this HLT filter, and place it in the Event.
	//  if (saveTags()) filterproduct.addCollectionTag();
	//  edm::Ref<reco::RecoChargedCandidateCollection> candref;

	// get beam spot
	reco::BeamSpot vertexBeamSpot;
	edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
	iEvent.getByLabel(beamSpotTag_,recoBeamSpotHandle);
	vertexBeamSpot = *recoBeamSpotHandle;


	// get displaced vertices
	reco::VertexCollection displacedVertexColl;
	edm::Handle<reco::VertexCollection> displacedVertexCollHandle;
	bool foundVertexColl = iEvent.getByLabel(DisplacedVertexTag_, displacedVertexCollHandle);
	if(foundVertexColl) displacedVertexColl = *displacedVertexCollHandle;


	edm::Handle<reco::VertexCollection> vertices;
	iEvent.getByLabel(pvTag_ ,vertices);
	const reco::Vertex & vtx = (*vertices)[0];    //Use the first vertex. For PbPb collision, hiSelectedVertex should be just one

	bool triggered = false;
	VertexDistance3D a3d;
	VertexDistanceXY axy;

	// loop over vertex collection
	for(reco::VertexCollection::iterator it = displacedVertexColl.begin(); it!= displacedVertexColl.end(); it++)
	{
		reco::Vertex displacedVertex = *it; 
		if(displacedVertex.tracksSize() != 2) continue;
		float normChi2 = displacedVertex.normalizedChi2();
		if (normChi2 > maxNormalisedChi2_) continue;

		double vtxProb = 0.0;
		if( (displacedVertex.chi2()>=0.0) && (displacedVertex.ndof()>0) ) vtxProb = TMath::Prob(displacedVertex.chi2(), displacedVertex.ndof() );
		if (vtxProb < minVtxProbability_) continue;

		// get the two tracks from the vertex
		reco::Vertex::trackRef_iterator trackIt =  displacedVertex.tracks_begin();
		reco::TrackRef vertextkRef1 =  (*trackIt).castTo<reco::TrackRef>() ;

		trackIt++;
		reco::TrackRef vertextkRef2 =  (*trackIt).castTo<reco::TrackRef>();

		//point angle cut
		TVector3 pperp(vertextkRef1->px() + vertextkRef2->px(), vertextkRef1->py() + vertextkRef2->py(), 0.);
		TVector3 pvtosv(displacedVertex.x() - vtx.x(), displacedVertex.y() - vtx.y(), 0.0);
		double agl = TMath::Cos(pperp.Angle(pvtosv));
		//	  cout << " agl : " << agl << endl;
		if(agl < minCosinePointingAngle_) continue;

		if( decaylength3D_ )
		{
			double ffxyz = a3d.distance(vtx, displacedVertex).value();
			double ffxyzerror = a3d.distance(vtx, displacedVertex).error();
			if( ffxyz / ffxyzerror < minLxyzSignificance_ )     continue;
		}
		else
		{
			double ffxy = axy.distance(vtx, displacedVertex).value();
			double ffxyerror = axy.distance(vtx, displacedVertex).error();
			if( ffxy / ffxyerror < minLxySignificance_ )     continue; 
		}

		triggered = true;

	}

	LogDebug("HLTDisplacedMumumuFilter") << " >>>>> Result of HLTDisplacedMumumuFilter is "<< triggered <<std::endl; 

	return triggered;
}

DEFINE_FWK_MODULE(HLTDisplacedDmesonVtxFilter);

