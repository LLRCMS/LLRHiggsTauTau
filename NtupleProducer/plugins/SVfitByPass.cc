// bypasses SVfit and saves the same userfloats
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
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/MET.h>
#include <DataFormats/METReco/interface/PFMET.h>
#include <DataFormats/METReco/interface/PFMETCollection.h>
#include <DataFormats/METReco/interface/CommonMETData.h>
#include <TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h>

#include <vector>
#include <string>

using namespace edm;
using namespace std;
using namespace reco;

// ------------------------------------------------------------------

class SVfitBypass : public edm::EDProducer {
 public:
  /// Constructor
  explicit SVfitBypass(const edm::ParameterSet&);
    
  /// Destructor
  ~SVfitBypass(){
  };  

 private:
  virtual void beginJob(){};  
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};
  
  edm::InputTag theCandidateTag;
  std::vector<edm::InputTag> vtheMETTag;
  //int sampleType;
  //bool _usePairMET;
  //bool _useMVAMET;
};

// ------------------------------------------------------------------



SVfitBypass::SVfitBypass(const edm::ParameterSet& iConfig)
{
  theCandidateTag = iConfig.getParameter<InputTag>("srcPairs");
  //_usePairMET = iConfig.getUntrackedParameter<bool>("usePairMET");
  //_useMVAMET = iConfig.getUntrackedParameter<bool>("useMVAMET");

  //vtheMETTag = iConfig.getParameter<std::vector<edm::InputTag>>("srcMET");

  produces<pat::CompositeCandidateCollection>();
}  



void SVfitBypass::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{  

  // Get lepton pairs
  Handle<View<reco::CompositeCandidate> > pairHandle;
  iEvent.getByLabel(theCandidateTag, pairHandle);
  
  unsigned int elNumber = pairHandle->size();
  
    // Output collection
  auto_ptr<pat::CompositeCandidateCollection> result( new pat::CompositeCandidateCollection );

  // loop on all the pairs
  for (unsigned int i = 0; i < elNumber; ++i)
  {
    // Get the pair and the two leptons composing it
    const CompositeCandidate& pairBuf = (*pairHandle)[i];
    pat::CompositeCandidate pair(pairBuf);

    float SVfitMass = -999.;
    // add user floats: SVfit mass, met properties, etc..  
    pair.addUserFloat("SVfitMass", SVfitMass);
    pair.addUserFloat("SVfit_pt",  -999);
    pair.addUserFloat("SVfit_eta", -999);
    pair.addUserFloat("SVfit_phi", -999);
    pair.addUserFloat("SVfit_ptUnc",  -999);
    pair.addUserFloat("SVfit_etaUnc",  -999);
    pair.addUserFloat("SVfit_phiUnc",  -999);
    pair.addUserFloat("MEt_px", -999);
    pair.addUserFloat("MEt_py", -999);
    pair.addUserFloat("MEt_cov00", -999);
    pair.addUserFloat("MEt_cov01", -999);
    pair.addUserFloat("MEt_cov10", -999);
    pair.addUserFloat("MEt_cov11", -999);
    pair.addUserFloat("MEt_significance", -999);

    result->push_back(pair);
  }
  iEvent.put(result);
}

#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(SVfitBypass);

