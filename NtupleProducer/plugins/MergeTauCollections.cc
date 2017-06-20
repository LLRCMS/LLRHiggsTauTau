#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <DataFormats/PatCandidates/interface/Tau.h>
//#include <DataFormats/TauReco/interface/PFTauDiscriminator.h>

//#include "DataFormats/VertexReco/interface/Vertex.h"
#include <LLRHiggsTauTau/NtupleProducer/interface/DaughterDataHelpers.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/CutSet.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/LeptonIsoHelper.h>
//#include "BDTId.h"

#include <vector>
#include <string>

using namespace edm;
using namespace std;
using namespace reco;

//bool recomputeBDT = false;

class MergeTauCollections : public edm::EDProducer {
 public:
  /// Constructor
  explicit MergeTauCollections(const edm::ParameterSet&);
    
  /// Destructor
  ~MergeTauCollections(){
    
  };  
  //ByIsolationMVA3oldDMwoLTraw
 private:
  virtual void beginJob(){};  
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};

  //const edm::InputTag theCandidateTag;
  //const edm::InputTag theCandidateTagTauUp;
  //const edm::InputTag theCandidateTagTauDown;
  edm::EDGetTokenT<pat::TauCollection> theCandidateTag;
  edm::EDGetTokenT<pat::TauCollection> theCandidateTagTauUp;
  edm::EDGetTokenT<pat::TauCollection> theCandidateTagTauDown;

  

};

MergeTauCollections::MergeTauCollections(const edm::ParameterSet& iConfig) :
  theCandidateTag(consumes<pat::TauCollection>(iConfig.getParameter<InputTag>("src"))),
  theCandidateTagTauUp(consumes<pat::TauCollection>(iConfig.getParameter<InputTag>("srcTauUp"))),
  theCandidateTagTauDown(consumes<pat::TauCollection>(iConfig.getParameter<InputTag>("srcTauDown")))
{
  produces<pat::TauCollection>();
}


void MergeTauCollections::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{  
  edm::Handle<pat::TauCollection> tauHandle;
  iEvent.getByToken(theCandidateTag, tauHandle);

  // Output collection
  //auto_ptr<pat::TauCollection> result( new pat::TauCollection() );
  std::unique_ptr<pat::TauCollection> result( new pat::TauCollection() );

  edm::Handle<pat::TauCollection> tauHandleTauUp;
  iEvent.getByToken(theCandidateTagTauUp, tauHandleTauUp);

  edm::Handle<pat::TauCollection> tauHandleTauDown;
  iEvent.getByToken(theCandidateTagTauDown, tauHandleTauDown);

  //cout<<"trying to merge..."<<endl;

  for (unsigned int i = 0; i< tauHandle->size(); ++i){
    //---Clone the pat::Tau

    const pat::Tau& l = (*tauHandle)[i];
    pat::Tau l_Nominal(l);

    if(tauHandleTauUp->size()>i)
      {
	const pat::Tau& l1 = (*tauHandleTauUp)[i];
	pat::Tau l_TauUp(l1);
	
	l_Nominal.addUserFloat("px_TauUp",l_TauUp.px());
	l_Nominal.addUserFloat("py_TauUp",l_TauUp.py());
	l_Nominal.addUserFloat("pz_TauUp",l_TauUp.pz());
	l_Nominal.addUserFloat("e_TauUp",l_TauUp.p4().T());
	if(l1.userInt("isTESShifted")) l_Nominal.addUserInt("TauUpExists",1);
	else l_Nominal.addUserInt("TauUpExists",0);
      }
    else
      {
	l_Nominal.addUserFloat("px_TauUp",-999.);
	l_Nominal.addUserFloat("py_TauUp",-999.);
	l_Nominal.addUserFloat("pz_TauUp",-999.);
	l_Nominal.addUserFloat("e_TauUp",-999.);
	l_Nominal.addUserInt("TauUpExists",0);
      }

    if(tauHandleTauDown->size()>i)
      {
	const pat::Tau& l2 = (*tauHandleTauDown)[i];
	pat::Tau l_TauDown(l2);

	l_Nominal.addUserFloat("px_TauDown",l_TauDown.px());
	l_Nominal.addUserFloat("py_TauDown",l_TauDown.py());
	l_Nominal.addUserFloat("pz_TauDown",l_TauDown.pz());
	l_Nominal.addUserFloat("e_TauDown",l_TauDown.p4().T());
	if(l2.userInt("isTESShifted")) l_Nominal.addUserInt("TauDownExists",1);
	else l_Nominal.addUserInt("TauDownExists",0);
      }
    else
      {
	l_Nominal.addUserFloat("px_TauDown",-999.);
	l_Nominal.addUserFloat("py_TauDown",-999.);
	l_Nominal.addUserFloat("pz_TauDown",-999.);
	l_Nominal.addUserFloat("e_TauDown",-999.);
	l_Nominal.addUserInt("TauDownExists",0);
      }

    result->push_back(l_Nominal);
  }

  //iEvent.put(result);
  iEvent.put(std::move(result));
}
 
#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(MergeTauCollections);
