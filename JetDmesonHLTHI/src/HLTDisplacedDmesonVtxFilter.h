#ifndef HLTDisplacedDmesonVtxFilter_h
#define HLTDisplacedDmesonVtxFilter_h

#include "HLTrigger/HLTcore/interface/HLTFilter.h"

namespace edm {
  class ConfigurationDescriptions;
}

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

	
class HLTDisplacedDmesonVtxFilter : public HLTFilter {
 public:
		explicit HLTDisplacedDmesonVtxFilter(const edm::ParameterSet&);
		~HLTDisplacedDmesonVtxFilter();
	
        static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);
		virtual void beginJob() ;
		virtual bool hltFilter(edm::Event&, const edm::EventSetup&, trigger::TriggerFilterObjectWithRefs & filterproduct) const override;
		virtual void endJob() ;
		
 private:

        bool   decaylength3D_;
		double minLxySignificance_;
		double minLxyzSignificance_;
		double maxNormalisedChi2_;
		double minVtxProbability_;
		double minCosinePointingAngle_;
		edm::InputTag DisplacedVertexTag_;
		edm::InputTag beamSpotTag_;
		edm::InputTag pvTag_;
};
#endif
