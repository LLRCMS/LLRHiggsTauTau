/* This plugins unpacks the LL pair that is at a certain position inside the pair vector
** and produces a collection containing these candidates.
** This collection can be given as input to the MVA MET producer.
**
** Inputs are the pair collection and the index of the position.
** If the position index exceeds the pair collection size, nothing is created.
**
** \date:    17 December 2014
** \author:  L. Cadamuro (LLR)
*/

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Utilities/interface/InputTag.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/CutSet.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/LeptonIsoHelper.h>
#include <DataFormats/Candidate/interface/ShallowCloneCandidate.h>
#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>

#include <vector>
#include <string>

using namespace edm;
using namespace std;
using namespace reco;

// ------------------------------------------------------------------

class PairUnpacker : public edm::EDProducer {

   public:
      // ctor
      explicit PairUnpacker(const edm::ParameterSet&);
   
      //dtor
      ~PairUnpacker() {};

   private:
      virtual void beginJob(){};  
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob(){};

      const edm::InputTag _InputPairsTag;
      const unsigned int _pairIndex;
};

// ------------------------------------------------------------------


PairUnpacker::PairUnpacker (const edm::ParameterSet& iConfig) :
   _InputPairsTag (iConfig.getParameter<InputTag>("src")),
   _pairIndex     (iConfig.getParameter<int> ("pairIndex"))
{
   consumes<View<reco::CompositeCandidate> >(_InputPairsTag);

   produces<reco::CandidateCollection>();
}


void PairUnpacker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   //cout << "Unpacking pairs at index: " << _pairIndex << endl;

   // Create handle on lepton pairs
   Handle<View<reco::CompositeCandidate> > pairHandle;
   iEvent.getByLabel(_InputPairsTag, pairHandle);
   
   // output collection
   //auto_ptr<reco::CandidateCollection> result( new reco::CandidateCollection );
   std::unique_ptr<reco::CandidateCollection> result( new reco::CandidateCollection );

   // if pairIndex exceeds Pair collection size, do not do anything
   if (_pairIndex < pairHandle->size() )
   {
      const CompositeCandidate& pair = (*pairHandle)[_pairIndex];
   
      const Candidate& l1 = *(pair.daughter(0));
      const Candidate& l2 = *(pair.daughter(1));

      result->push_back(l1);
      result->push_back(l2);

      //cout << "Finished unpacking at index: " << _pairIndex << endl;
   }
   //else cout << "Unpacked nothing, returning" << endl;

   //iEvent.put(result);
   iEvent.put(std::move(result));
}


#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(PairUnpacker);
