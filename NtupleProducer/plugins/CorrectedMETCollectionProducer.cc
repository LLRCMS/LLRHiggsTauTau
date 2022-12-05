#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Utilities/interface/InputTag.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/MET.h>
#include <DataFormats/Candidate/interface/Candidate.h>

#include <vector>
#include <string>

using namespace edm;
using namespace std;
using namespace reco;

class CorrectedMETCollectionProducer : public edm::EDProducer {
    public: 
        /// Constructor
        explicit CorrectedMETCollectionProducer(const edm::ParameterSet&);
        /// Destructor
        ~CorrectedMETCollectionProducer(){};

    private:
        virtual void beginJob(){};  
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob(){};

        edm::EDGetTokenT<View<pat::MET>> theMETTag;
        int correctionLevel;
};

CorrectedMETCollectionProducer::CorrectedMETCollectionProducer(const edm::ParameterSet& iConfig) :
theMETTag(consumes<View<pat::MET>>(iConfig.getParameter<edm::InputTag>("srcMET"))),
correctionLevel(iConfig.getParameter<int>("correctionLevel"))
{
    produces<pat::METCollection>();
}

void CorrectedMETCollectionProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    // Declare ptrs to save MET variations
    std::unique_ptr<pat::METCollection> out_MET_ptr(new pat::METCollection());

    // Get the MET
    Handle<View<pat::MET> > METHandle;
    iEvent.getByToken(theMETTag, METHandle);
    const pat::MET& patMET = (*METHandle)[0];

    reco::Candidate::PolarLorentzVector correctedMETP4(patMET.corPt(pat::MET::METCorrectionLevel(correctionLevel)),
                                                  0.,
                                                  patMET.corPhi(pat::MET::METCorrectionLevel(correctionLevel)), 
                                                  0.);

    pat::MET correctedMET(patMET);
    correctedMET.setP4(correctedMETP4);

    out_MET_ptr->push_back(correctedMET);

    iEvent.put(std::move(out_MET_ptr));
}

#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(CorrectedMETCollectionProducer);
