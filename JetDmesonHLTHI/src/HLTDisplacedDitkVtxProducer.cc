// -*- C++ -*-
//
// Package:    JetDmesonHLTHI
// Class:      HLTDisplacedDitkVtxProducer
//
// D meson secondary vertex producer for HI D meson triggers
//reco D mesons in Jet cone
//changed from "HLTDisplacedmumuVtxProducer"
//
// Original Author:  Jian Sun
//         Created:  Wed Nov. 5 2014
// $Id$
//
//

#include <iostream>

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerRefsCollections.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/RefToBase.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/JetReco/interface/JetTracksAssociation.h"
#include "HeavyFlavorHLTHI/JetDmesonHLTHI/interface/HFMasses.hh"

#include "HLTDisplacedDitkVtxProducer.h"
#include <TLorentzVector.h>

using namespace edm;
using namespace reco;
using namespace std; 
using namespace trigger;
//
// constructors and destructor
//
HLTDisplacedDitkVtxProducer::HLTDisplacedDitkVtxProducer(const edm::ParameterSet& iConfig):	
	src_ (iConfig.getParameter<edm::InputTag>("Src")),
	maxEta_ (iConfig.getParameter<double>("MaxEta")),
	mintrackPt_ (iConfig.getParameter<double>("MinTrackPt")),
	minPtPair_ (iConfig.getParameter<double>("MinPtPair")),
	minInvMass_ (iConfig.getParameter<double>("MinInvMass")),
	maxInvMass_ (iConfig.getParameter<double>("MaxInvMass"))
{
	produces<VertexCollection>();
}


HLTDisplacedDitkVtxProducer::~HLTDisplacedDitkVtxProducer()
{

}

void HLTDisplacedDitkVtxProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	edm::ParameterSetDescription desc;
	desc.add<edm::InputTag>("Src",edm::InputTag("hltFastPixelBLifetimeL3AssociatorHbb"));
	desc.add<double>("MaxEta",2.5);
	desc.add<double>("MinTrackPt",2.0);
	desc.add<double>("MinPtPair",5.0);
	desc.add<double>("MinInvMass",1.55);
	desc.add<double>("MaxInvMass",2.10);
	descriptions.add("hltDisplacedmumumuVtxProducer", desc);
}

// ------------ method called once each job just before starting event loop  ------------
void HLTDisplacedDitkVtxProducer::beginJob()
{

}

// ------------ method called once each job just after ending the event loop  ------------
void HLTDisplacedDitkVtxProducer::endJob() 
{

}

// ------------ method called on each new Event  ------------
void HLTDisplacedDitkVtxProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

	// get hold of muon trks
	Handle<reco::JetTracksAssociationCollection> jetTracksAssociation;
	iEvent.getByLabel(src_,jetTracksAssociation);
	//	iEvent.getByToken (src_,jetTracksAssociation); ???  //should change to getByToken in future, it is faster

	//get the transient track builder:
	edm::ESHandle<TransientTrackBuilder> theB;
	iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

	std::auto_ptr<VertexCollection> vertexCollection(new VertexCollection());


	for(JetTracksAssociationCollection::const_iterator it = jetTracksAssociation->begin(); it != jetTracksAssociation->end(); it++)
	{
		TrackRefVector tracks = it->second;
		const RefToBase<reco::Jet> jet = it->first;
		cout << " jet pt: " << jet->pt() << "  eta: " << jet->eta() << "  phi: " << jet->phi() << endl;

		//     math::XYZVector jetMomentum = it->first->momentum();
		TLorentzVector ka,pi,trackpair;

		for(TrackRefVector::const_iterator itTrack_pion = tracks.begin(); itTrack_pion != tracks.end(); ++itTrack_pion)
		{

			const Track & track_pion = **itTrack_pion;
			if( track_pion.pt() < mintrackPt_ )   continue;

			for(TrackRefVector::const_iterator itTrack_kaon = tracks.begin(); itTrack_kaon != tracks.end(); ++itTrack_kaon)
			{
				if( itTrack_pion == itTrack_kaon )  continue;
				const Track & track_kaon = **itTrack_kaon;
				if( track_pion.charge() * track_kaon.charge() > 0 )   continue;
				if( track_kaon.pt() < mintrackPt_ )   continue;
				pi.SetPtEtaPhiM(track_pion.pt(), track_pion.eta(), track_pion.phi(), MPION);
				ka.SetPtEtaPhiM(track_kaon.pt(), track_kaon.eta(), track_kaon.phi(), MKAON);
				trackpair = ka + pi;

				if( trackpair.Pt() < 20.0 )  continue;

                cout << "trackpair pt: " << trackpair.Pt() << "  pion pt: " << track_pion.pt() << "   kaon pt: " << track_kaon.pt() << endl;

				if (  trackpair.M() <  minInvMass_ || trackpair.M() >   maxInvMass_) continue;
				if ( trackpair.Pt() < minPtPair_)     continue;
				if ( TMath::Abs( trackpair.Eta() ) > maxEta_ )     continue;


				// do the vertex fit
				vector<TransientTrack> t_tks;
				TransientTrack ttkp1 = (*theB).build(&track_pion);
				TransientTrack ttkp2 = (*theB).build(&track_kaon);
				t_tks.push_back(ttkp1);
				t_tks.push_back(ttkp2);

				KalmanVertexFitter kvf;
				TransientVertex tv = kvf.vertex(t_tks);

				if (!tv.isValid()) continue;

				Vertex vertex = tv;
				vertex.removeTracks();    //by default, the dau tracks are not saved but the track size is saved.
				vertex.add(reco::TrackBaseRef(*itTrack_pion));// Need to reset the track size and then add the dau tracks
				vertex.add(reco::TrackBaseRef(*itTrack_kaon));
				// put vertex in the event
				vertexCollection->push_back(vertex);

			}

		}
	}


	iEvent.put(vertexCollection);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HLTDisplacedDitkVtxProducer);
