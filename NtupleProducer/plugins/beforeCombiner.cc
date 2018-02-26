#ifndef BEFORECOMBINER_H
#define BEFORECOMBINER_H

#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h>
#include <FWCore/Utilities/interface/InputTag.h>
#include <DataFormats/PatCandidates/interface/Muon.h>

#include <iostream>


using namespace edm;
using namespace std;
// using namespace reco;

class beforeCombiner : public edm::EDFilter {

    public:
        beforeCombiner(const edm::ParameterSet &);
        ~beforeCombiner();

    private:
        bool filter(edm::Event &, edm::EventSetup const&);
        edm::EDGetTokenT<edm::View<reco::Candidate>> _candTag; // softLeptons
};

beforeCombiner::beforeCombiner(const edm::ParameterSet & iConfig) :
_candTag  (consumes<edm::View<reco::Candidate> > (iConfig.getParameter<InputTag>("src")))
{}

beforeCombiner::~beforeCombiner()
{}

bool beforeCombiner::filter(edm::Event & iEvent, edm::EventSetup const& iSetup)
{
    std::cout << " ----- New Event -----" << endl;
    edm::Handle<edm::View<reco::Candidate>> candHandle;
    iEvent.getByToken(_candTag, candHandle);
    
    const edm::View<reco::Candidate>* daus = candHandle.product();
    std::cout << " **Before** Leptons size:" << daus->size() << endl;
    
    for (uint ilep = 0; ilep < daus->size(); ++ilep)
    {
        std::cout << "     Lepton " << ilep << " pdgId: " << daus->at(ilep).pdgId() << endl;
    }

    return true;
}

#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(beforeCombiner);

#endif
