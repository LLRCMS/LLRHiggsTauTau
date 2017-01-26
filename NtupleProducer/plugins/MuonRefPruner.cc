/*
** class  : MuonRefPruner
** author : L. Cadamuro (LLR)
** date   : 24 January 2017
** brief  : select only the references in the input collection that are not in the specified collection
**          input is a View of muons
**          output is a collection of pat::Muons
*/

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Utilities/interface/InputTag.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/MuonReco/interface/Muon.h>

#include <vector>
#include <string>

using namespace edm;
using namespace std;
using namespace reco;

class MuonRefPruner : public edm::EDProducer {
    public: 
        /// Constructor
        explicit MuonRefPruner(const edm::ParameterSet&);
        /// Destructor
        ~MuonRefPruner(){};  

    private:
        virtual void beginJob(){};  
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob(){};

        edm::EDGetTokenT<pat::MuonCollection> input_;
        edm::EDGetTokenT<edm::View<reco::Muon>> toremove_;
        edm::EDGetTokenT<edm::View<reco::Muon>> toremove2_;
};

MuonRefPruner::MuonRefPruner(const edm::ParameterSet & iConfig) :
    input_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("input"))),
    toremove_(consumes<edm::View<reco::Muon>>(iConfig.getParameter<edm::InputTag>("toremove"))),
    toremove2_(consumes<edm::View<reco::Muon>>(iConfig.getParameter<edm::InputTag>("toremove2")))
{
    // produces<edm::PtrVector<pat::Muon>>();
    produces<vector<pat::Muon>>();
    produces<int>();
}

void MuonRefPruner::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    Handle<pat::MuonCollection> inputHandle;
    Handle<View<reco::Muon> > toremoveHandle;
    Handle<View<reco::Muon> > toremove2Handle;
    iEvent.getByToken(input_, inputHandle);
    iEvent.getByToken(toremove_, toremoveHandle);
    iEvent.getByToken(toremove2_, toremove2Handle);

    // std::unique_ptr<edm::PtrVector<pat::Muon>> out(new edm::PtrVector<pat::Muon>());
    std::unique_ptr<vector<pat::Muon>> out(new pat::MuonCollection());
    std::unique_ptr<int> nbad (new int(0));

    for (unsigned int i = 0; i < inputHandle->size(); ++i)
    {
        bool good = true;
        edm::Ptr<reco::Muon> myPtr(inputHandle, i);
        for (unsigned int j = 0; j < toremoveHandle->size(); ++j)
        {
            if (myPtr == toremoveHandle->ptrAt(j))
            {
                good = false;
                break;
            }
        }

        for (unsigned int j = 0; j < toremove2Handle->size(); ++j)
        {
            if (myPtr == toremove2Handle->ptrAt(j))
            {
                good = false;
                break;
            }
        }

        if(good)
        {
            pat::Muon m(inputHandle->at(i));
            out->push_back(m);
            // edm::Ptr<pat::Muon> pp = (edm::Ptr<pat::Muon>) inputHandle->ptrAt(i);
            // out->push_back(pp);
        }

        else
        {
            (*nbad) += 1;
            // std::cout << " >>> rejected muon " << i << " of pt=" << inputHandle->at(i).pt() << " eta=" << inputHandle->at(i).eta() << " phi=" << inputHandle->at(i).phi() << endl;
        }
    }

    iEvent.put(std::move(out));
    iEvent.put(std::move(nbad));
}

#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(MuonRefPruner);
