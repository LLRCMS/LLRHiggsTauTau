/*
** class  : DeepMET
** author : T. Kramer (UHH)
** date   : 18 November 2020
** brief  : takes in input the met, and produces a pat::METCollection with the pt and phi set to the 
**          DeepMET values
*/

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Utilities/interface/InputTag.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/MET.h>
// #include <DataFormats/METReco/interface/PFMET.h>
// #include <DataFormats/METReco/interface/PFMETCollection.h>
// #include <DataFormats/METReco/interface/CommonMETData.h>
// #include "TLorentzVector.h"
// #include <LLRHiggsTauTau/NtupleProducer/interface/DaughterDataHelpers.h>
#include <DataFormats/Candidate/interface/Candidate.h>

#include <vector>
#include <string>

using namespace edm;
using namespace std;
using namespace reco;

class DeepMET : public edm::EDProducer {
    public: 
        /// Constructor
        explicit DeepMET(const edm::ParameterSet&);
        /// Destructor
        ~DeepMET(){};

    private:
        virtual void beginJob(){};  
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob(){};

        edm::EDGetTokenT<View<pat::MET>> theMETTag;
        int tune;
};

DeepMET::DeepMET(const edm::ParameterSet& iConfig) :
theMETTag(consumes<View<pat::MET>>(iConfig.getParameter<edm::InputTag>("srcMET"))),
tune(iConfig.getParameter<int>("tune"))
{
    produces<pat::METCollection>();
}

void DeepMET::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    // Declare ptrs to save MET variations
    std::unique_ptr<pat::METCollection> out_MET_ptr(new pat::METCollection());

    // Get the MET
    Handle<View<pat::MET> > METHandle;
    iEvent.getByToken(theMETTag, METHandle);
    const pat::MET& patMET = (*METHandle)[0];

    reco::Candidate::PolarLorentzVector deepMETP4(patMET.corPt(pat::MET::METCorrectionLevel(tune)),
                                                  0.,
                                                  patMET.corPhi(pat::MET::METCorrectionLevel(tune)), 
                                                  0.);

    pat::MET deepMET(patMET);
    deepMET.setP4(deepMETP4);

    out_MET_ptr->push_back(deepMET);

    iEvent.put(std::move(out_MET_ptr));
}

#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(DeepMET);
