
#ifndef HLTDzerosvFilter_h
#define HLTDzerosvFilter_h

#include "HLTrigger/HLTcore/interface/HLTFilter.h"

namespace edm {
  class ConfigurationDescriptions;
}

//
// class declaration
//

class HLTDzerosvFilter : public HLTFilter {

   public:
      explicit HLTDzerosvFilter(const edm::ParameterSet&);
      ~HLTDzerosvFilter();
      static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);
      virtual bool hltFilter(edm::Event&, const edm::EventSetup&, trigger::TriggerFilterObjectWithRefs & filterproduct) const override;

   private:
//      edm::InputTag pixelVerticesTag_;  // input tag identifying product containing Pixel-vertices
      edm::InputTag svTagInfoProducer_;
      edm::InputTag hltTracksTag_;  // input tag identifying product containing Pixel-tracks
	  edm::InputTag caloJetTag_;      // input tag for jet collection

      double min_Pt_;          // min pt cut
      double max_Pt_;          // max pt cut
      double max_Eta_;          // max eta cut
      double max_Vz_;          // max vz cut
      int min_trks_;  // minimum number of tracks from one vertex
      float masswindow;          // minimum separation of two tracks in phi-eta
};

#endif //HLTDzerosvFilter_h
