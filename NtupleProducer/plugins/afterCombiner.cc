#ifndef AFTERCOMBINER_H
#define AFTERCOMBINER_H

#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h>
#include <FWCore/Utilities/interface/InputTag.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>

#include <iostream>

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h>
#include <FWCore/Utilities/interface/InputTag.h>
#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/MET.h>

using namespace edm;
using namespace std;
using namespace reco;

typedef edm::View<reco::Candidate> CandidateView;

class afterCombiner : public edm::EDFilter {

    public:
        afterCombiner(const edm::ParameterSet &);
        ~afterCombiner();

    private:
        bool filter(edm::Event &, edm::EventSetup const&);
        edm::InputTag srcPairs_;
};

afterCombiner::afterCombiner(const edm::ParameterSet & iConfig) //:
{
  srcPairs_ = iConfig.getParameter<edm::InputTag>("srcPairs");
  consumes<CompositeCandidateView>(srcPairs_);
}

afterCombiner::~afterCombiner()
{}

bool afterCombiner::filter(edm::Event & iEvent, edm::EventSetup const& iSetup)
{
    edm::Handle<CompositeCandidateView> pairs;
    iEvent.getByLabel(srcPairs_, pairs);
    
    unsigned int nPairs = pairs->size();
    std::cout << " ** After** Pairs size:" << nPairs << endl;
    for (unsigned int iPair = 0; iPair < nPairs; ++iPair)
    {
        cout << "     Pair " << iPair << endl;
        const CompositeCandidate& pair = (*pairs)[iPair];
        const reco::Candidate *l1 = pair.daughter(0);
        const reco::Candidate *l2 = pair.daughter(1);
        cout << "        l1: " << l1->pdgId() << endl;
        cout << "        l2: " << l2->pdgId() << endl;
    
    }
    
    return true;
}

#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(afterCombiner);

#endif
