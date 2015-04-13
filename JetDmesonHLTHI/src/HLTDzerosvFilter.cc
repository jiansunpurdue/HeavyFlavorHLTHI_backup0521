// -*- C++ -*-
//
// Package:    JetDmesonHLTHI
// Class:      HLTDzerosvFilter
//
// secondary vertex filter, cutting on invariant mass of tracks from sv
// secondary vertex are produced with sv producer in B-Jet trigger sequence
//
// Original Author:  Jian Sun
//         Created:  Wed Nov. 5 2014
// $Id$
//
//
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h" 
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"
#include "DataFormats/Common/interface/Ref.h"
#include "Math/GenVector/VectorUtil.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"

#include "HeavyFlavorHLTHI/JetDmesonHLTHI/interface/HFMasses.hh"

#include "HLTDzerosvFilter.h"
#include <TLorentzVector.h>


// constructors and destructor
//
using namespace std; 
using namespace reco;
using namespace edm;

HLTDzerosvFilter::HLTDzerosvFilter(const edm::ParameterSet& iConfig) : HLTFilter(iConfig),
								       //    pixelVerticesTag_ (iConfig.getParameter<edm::InputTag>("vertexCollection")),
								       svTagInfoProducer_(iConfig.getParameter<edm::InputTag> ( "svTagInfoProducer" )),
								       hltTracksTag_ (iConfig.getParameter<edm::InputTag>("trackCollection")),
								       caloJetTag_(iConfig.getParameter<edm::InputTag>("jetCollection")),
								       min_Pt_  (iConfig.getParameter<double>("MinPt")),
								       max_Pt_  (iConfig.getParameter<double>("MaxPt")),
								       max_Eta_  (iConfig.getParameter<double>("MaxEta")),
								       max_Vz_  (iConfig.getParameter<double>("MaxVz")),
								       min_trks_  (iConfig.getParameter<int>("MinTrks")),
								       masswindow  (iConfig.getParameter<double>("Masswindow"))
{
}

HLTDzerosvFilter::~HLTDzerosvFilter()
{
}

void HLTDzerosvFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  makeHLTFilterDescription(desc);
  desc.add<double>("MinPt",0.0);
  desc.add<double>("MaxPt",1000.0);
  desc.add<double>("MaxEta",2.2);
  desc.add<double>("MaxVz",30.0);
  desc.add<int>("MinTrks",1);
  desc.add<double>("Masswindow",0.3);
  desc.add<edm::InputTag>("svTagInfoProducer",edm::InputTag("hltL3SecondaryVertexTagInfos"));
  desc.add<edm::InputTag>("trackCollection",edm::InputTag("hltIter4Merged"));
  desc.add<edm::InputTag>("jetCollection",edm::InputTag("hltSelectorJets20L1FastJet"));
  descriptions.add("hltDisplacedmumuFilter", desc);
 }

//
// member functions
//

// ------------ method called to produce the data  ------------
bool HLTDzerosvFilter::hltFilter(edm::Event& iEvent, const edm::EventSetup& iSetup, trigger::TriggerFilterObjectWithRefs & filterproduct) const
{
	edm::Handle<reco::SecondaryVertexTagInfoCollection> svTagInfoCollection;
	iEvent.getByLabel(svTagInfoProducer_, svTagInfoCollection);

	bool accept = false;

    for (std::size_t index = 0; index < svTagInfoCollection->size(); ++index)
	{
		reco::SecondaryVertexTagInfoRef svTagInfo(svTagInfoCollection, index);
//		RefToBase<Jet> jet = svTagInfo->jet();
//		cout << " jet pt: " << jet->pt() << endl;
		// Loop over the vertexes in svTagInfo
		int num_d0candidate = 0;
		for ( std::size_t vindex = 0; vindex < svTagInfo->nVertices(); ++vindex )
		{
			const Vertex& secVertex = svTagInfo->secondaryVertex(vindex);
//			Measurement1D m1D = svTagInfo->flightDistance(vindex);
//			cout << "sig: " << m1D.significance() << endl;
			TLorentzVector ka,pi,d0;
			for( std::vector< TrackBaseRef >::const_iterator track_pion = secVertex.tracks_begin(); track_pion != secVertex.tracks_end(); track_pion++ )
			{
				for( std::vector< TrackBaseRef >::const_iterator track_kaon = secVertex.tracks_begin(); track_kaon != secVertex.tracks_end(); track_kaon++ )
				{
					if( track_pion == track_kaon )   continue;
					const TrackBaseRef& trackRef_pion = *track_pion;
					const TrackBaseRef& trackRef_kaon = *track_kaon;
					if( trackRef_pion->charge() * trackRef_kaon->charge() > 0 )  continue;
					pi.SetPtEtaPhiM(trackRef_pion->pt(), trackRef_pion->eta(), trackRef_pion->phi(), MPION);
					ka.SetPtEtaPhiM(trackRef_kaon->pt(), trackRef_kaon->eta(), trackRef_kaon->phi(), MKAON);
					d0 = ka + pi;
					if ( TMath::Abs( d0.M() - MD0 ) > masswindow ) continue;
					if ( d0.Pt() < min_Pt_ || d0.Pt() > max_Pt_)     continue;

					num_d0candidate++;

				}
			}

			
		}
		if ( num_d0candidate > 0 )  accept = true;         //write trig obj here
	}


   return accept;
}

DEFINE_FWK_MODULE(HLTDzerosvFilter);
