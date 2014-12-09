/** \class TauFiller
 *
 *  No description available.
 *
 *  $Date: 2013/05/24 15:42:42 $
 *  $Revision: 1.28 $
 *  \author G. Ortona (LLR)
 */

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <DataFormats/PatCandidates/interface/Tau.h>
//#include <DataFormats/TauReco/interface/PFTauDiscriminator.h>

//#include "DataFormats/VertexReco/interface/Vertex.h"

#include <LLRHiggsTauTau/NtupleProducer/interface/CutSet.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/LeptonIsoHelper.h>
//#include "BDTId.h"

#include <vector>
#include <string>

using namespace edm;
using namespace std;
using namespace reco;

//bool recomputeBDT = false;

class TauFiller : public edm::EDProducer {
 public:
  /// Constructor
  explicit TauFiller(const edm::ParameterSet&);
    
  /// Destructor
  ~TauFiller(){
    
  };  
  //ByIsolationMVA3oldDMwoLTraw
 private:
  virtual void beginJob(){};  
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};

  const edm::InputTag theCandidateTag;
  //const edm::InputTag theDiscriminatorTag;
  const std::string theDiscriminatorTag;
  const StringCutObjectSelector<pat::Tau, true> cut;
  const CutSet<pat::Tau> flags;
  //BDTId* bdt;
};


TauFiller::TauFiller(const edm::ParameterSet& iConfig) :
  theCandidateTag(iConfig.getParameter<InputTag>("src")),
  theDiscriminatorTag(iConfig.getParameter<std::string>("discriminator")),
  cut(iConfig.getParameter<std::string>("cut")),
  flags(iConfig.getParameter<ParameterSet>("flags"))//, 
  //bdt(0)
{
  //if (recomputeBDT) bdt = new BDTId;
  produces<pat::TauCollection>();
}


void
TauFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{  

//read one PFTauDiscriminator (set discriminatorSrc_ in to an edm::InputTag before)

// loop over taus
//    for ( unsigned iTau = 0; iTau < taus->size(); ++iTau ) {
//        reco::PFTauRef tauCandidate(taus, iTau);
// check if tau candidate has passed discriminator
//        if( (*discriminator)[tauCandidate] > 0.5 ){
//        // do something with your candidate
//        }
//    }


  // Get leptons and discriminators
  //edm::Handle<pat::PFTauCollection> tauHandle;
  edm::Handle<pat::TauRefVector> tauHandle;
  iEvent.getByLabel(theCandidateTag, tauHandle);
  
  //edm::Handle<reco::PFTauDiscriminator> discriminator;
  //iEvent.getByLabel(theDiscriminatorTag, discriminator);

  //hpsPFTauDiscrimination
  
  edm::Handle<vector<Vertex> >  vertexs;
  iEvent.getByLabel("goodPrimaryVertices",vertexs);

  // Output collection
  auto_ptr<pat::TauCollection> result( new pat::TauCollection() );

  for (unsigned int i = 0; i< tauHandle->size(); ++i){
    //---Clone the pat::Tau
    pat::Tau l(*((*tauHandle)[i].get()));
    
    
    //--- PF ISO
    float PFChargedHadIso   = l.chargedHadronIso();
    float PFNeutralHadIso   = l.neutralHadronIso();
    float PFPhotonIso       = l.photonIso();

    float combRelIsoPF = LeptonIsoHelper::combRelIsoPF(l);

    //Decay mode
    //int decayMode = -1;
    //int A = l.signalPFChargedHadrCands().size();
    //int B = l.signalPFGammaCands().size();
    //if(A==1&&B==0)decayMode=1;
    //else if(A==1&&B>0)decayMode=2;
    //else if (A==3&&B==0)decayMode=3;
    float tauid=l.tauID(theDiscriminatorTag);
    //printf("A, B, tau %d %d %f \n",A,B,tauid);

    l.addUserFloat("HPSDiscriminator",tauid);
    l.addUserFloat("decayMode",l.decayMode());
    //if(decayMode<0&&tauid==0)edm::LogWarning("TauFiller: Unrecognized decay mode");
    /*

    //--- SIP, dxy, dz
    float IP      = fabs(l.dB(pat::Electron::PV3D));
    float IPError = l.edB(pat::Electron::PV3D);
    float SIP     = IP/IPError;
    */
    
    float dxy = 999.;
    float dz  = 999.;
    const Vertex* vertex = 0;
    if (vertexs->size()>0) {
      vertex = &(vertexs->front());
      dxy = l.dxy();
      if(l.hasSecondaryVertex())dz  = l.secondaryVertex().get()->z()-vertex->z();
    } 
    
    
    //--- Embed user variables
    l.addUserFloat("dxy",dxy);
    l.addUserFloat("dz",dz);
    l.addUserFloat("PFChargedHadIso",PFChargedHadIso);
    l.addUserFloat("PFNeutralHadIso",PFNeutralHadIso);
    l.addUserFloat("PFPhotonIso",PFPhotonIso);
    l.addUserFloat("combRelIsoPF",combRelIsoPF);

    //--- MC parent code 
    //     MCHistoryTools mch(iEvent);
    //     if (mch.isMC()) {
    //       int MCParentCode = 0;
    //       //      int MCParentCode = mch.getParentCode(&l); //FIXME: does not work on cmg
    //       l.addUserFloat("MCParentCode",MCParentCode);
    //     }
    
    //--- Check selection cut. Being done here, flags are not available; but this way we 
    //    avoid wasting time on rejected leptons.
    if (!cut(l)) continue;

    //--- Embed flags (ie flags specified in the "flags" pset)
    for(CutSet<pat::Tau>::const_iterator flag = flags.begin(); flag != flags.end(); ++flag) {
      l.addUserFloat(flag->first,int((*(flag->second))(l)));
    }
    
    result->push_back(l);
  }
  iEvent.put(result);
}


#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(TauFiller);

