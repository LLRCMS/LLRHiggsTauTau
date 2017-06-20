// -*- C++ -*-
//
/*
 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Giacomo Ortona
//         Created:  Tue Jan 27 11:37:27 CEST 2015
//


#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <DataFormats/PatCandidates/interface/GenericParticle.h>
#include <DataFormats/PatCandidates/interface/PackedGenParticle.h>
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

//
// class declaration
//
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"

//
// class declaration
//

class bFiller : public edm::EDProducer {
public:
  explicit bFiller(const edm::ParameterSet&);
  ~bFiller();
  
private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  //const edm::InputTag theCandidateTag;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > theCandidateTag ;
  const StringCutObjectSelector<pat::GenericParticle, true> cut;
  const CutSet<pat::GenericParticle> flags;
  
  // ----------member data ---------------------------
  bool status2_,status3_;
};

// constructors and destructor
//
bFiller::bFiller(const edm::ParameterSet& iConfig):
  //theCandidateTag(iConfig.getParameter<InputTag>("src")),
  theCandidateTag(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("src"))),
  cut(iConfig.getParameter<std::string>("cut")),
  flags(iConfig.getParameter<ParameterSet>("flags"))//, 
{
  //now do what ever initialization is needed
  produces<pat::GenericParticleCollection>();
}


bFiller::~bFiller()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
void bFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  Bool_t isb_;
  Bool_t isLast_;
  //Bool_t isStatus2_,isStatus3_;
  //Bool_t hasHdaughter_;
    
  //edm::Handle<pat::PackedGenParticle> genHandle;
  //iEvent.getByLabel(theCandidateTag, genHandle);
  Handle<edm::View<reco::GenParticle> > genHandle;
  iEvent.getByToken(theCandidateTag, genHandle);

  // Output collection
  //auto_ptr<pat::GenericParticleCollection> result( new pat::GenericParticleCollection() );
  std::unique_ptr<pat::GenericParticleCollection> result( new pat::GenericParticleCollection() );
  //  printf("size %d\n",genHandle->size());
  //for(GenParticleCollection::const_iterator genp = genHandle->begin();genp != genHandle->end(); ++ genp ) {  // loop over GEN particles
  for (unsigned int i = 0; i< genHandle->size(); ++i){
    //Candidate cand(*((*genHandle)[i]));
    const GenParticle *packb =&(*genHandle)[i]; //(*((*genHandle)[i])));
    //pat::GenericParticle b((*genHandle)[i]);
 
    isb_=false;
    isLast_ = true; 
    //isStatus2_=false;
    //isStatus3_=false;
    
    //hasHdaughter_=false;
    int bPdgId = packb->pdgId();
    if ( fabs(bPdgId)==5 || fabs(bPdgId)==6 ) isb_=true;
    //cout << "Index: " << i << " isB: " << isb_ << "(id, pt: " << bPdgId << " , " << packb->pt() << " )" << endl;
    if (isb_)
    {
        for (unsigned int iDau = 0; iDau < packb->numberOfDaughters(); iDau++)
        {
            const Candidate * Dau = packb->daughter( iDau );
            //cout << "  Dau: " << iDau << " pdgId: " << Dau->pdgId() << endl;
            if (Dau->pdgId() == bPdgId) // reject "copies" of the same particle
            {
                isLast_ = false;
                break;
            }
        }
    }
   

    //if (packb->status()==2)isStatus2_=true;
    //if (packb->status()==3)isStatus3_=true;
    //if(isStatus2_ || isStatus3_)printf("%d\n",packb->pdgId());
    if(isb_ &&  isLast_ ){
      //cout << "Number of particle passed: " << i << " pt: " << packb->pt() << " eta: " << packb->eta() << endl;
      pat::GenericParticle b(*packb);   
      if (!cut(b)) continue;
      b.addUserFloat("status",packb->status());

      //Search if b comes from a H decay
      int motmass = 0;
      //      GenParticleRefVector *mothers = packb->mothers();
      //for (int im = 0; im<nmot; ++im){
      //	if(fabs(packb->mother(im)->pdgId()==25)){
      //	  motmass=1;
      //	  break;
      //	}
      //}
      //if(motmass==0){
      //int nmot = packb->numberOfMothers();
      //for (int im = 0; im<nmot&&motmass==0; ++im){
      for(unsigned int ipruned = 0; ipruned< genHandle->size(); ++ipruned){
          if(ipruned==i)continue;
          int pdgmot = (&(*genHandle)[ipruned])->pdgId();
          if(abs(pdgmot)==25){
    	      if(userdatahelpers::isAncestor(&(*genHandle)[ipruned],packb)){
	              //reco::Candidate* mot = &(*genHandle)[ipruned];
                  //int pdgmot = mot->pdgId();
                  //motmass=(&(*genHandle)[ipruned])->mass();
                  motmass = 1;
                  break;
              }
          }
      }
      //}
      
      b.addUserFloat("fromH",motmass);
      //--- Embed flags (ie flags specified in the "flags" pset) 
      for(CutSet<pat::GenericParticle>::const_iterator flag = flags.begin(); flag != flags.end(); ++flag) {
      	b.addUserFloat(flag->first,int((*(flag->second))(b)));
      }
    
      result->push_back(b);
    }

  }//end isb
  //iEvent.put(result);
  iEvent.put(std::move(result));
}//end loop over genParticles

// ------------ method called once each job just before starting event loop  ------------
void 
bFiller::beginJob()
{
  //edm::Service<TFileService> fileService;

    //h_genParticlesSize_ = fileService->make<TH1F>("h_genParticlesSize_","size of genParticles collection",2000,-0.5,1999.5); 
}

// ------------ method called once each job just after ending the event loop  ------------
void bFiller::endJob() {}

//define this as a plug-in
#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(bFiller);
