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

  vtheMETTag = iConfig.getParameter<std::vector<edm::InputTag>>("srcMET");

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

  // get event pat MET to be saved in output
  Handle<View<pat::MET> > METHandle_PatMET;
  iEvent.getByLabel(vtheMETTag.at(0), METHandle_PatMET);
  double METx = 0.;
  double METy = 0.; 
  TMatrixD covMET(2, 2);
  float significance = -999.;
  const pat::MET& patMET = (*METHandle_PatMET)[0];
  METx = patMET.px();
  METy = patMET.py();
  Handle<double> significanceHandle;
  Handle<math::Error<2>::type> covHandle;
  iEvent.getByLabel ("METSignificance", "METSignificance", significanceHandle);
  iEvent.getByLabel ("METSignificance", "METCovariance", covHandle);
  covMET[0][0] = (*covHandle)(0,0);
  covMET[1][0] = (*covHandle)(1,0);
  covMET[0][1] = covMET[1][0]; // (1,0) is the only one saved
  covMET[1][1] = (*covHandle)(1,1);
  significance = (float) (*significanceHandle);

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
    pair.addUserFloat("SVfit_METRho", -999);
    pair.addUserFloat("SVfit_METPhi", -999);
    pair.addUserFloat("MEt_px", (float) METx);
    pair.addUserFloat("MEt_py", (float) METy);
    pair.addUserFloat("MEt_cov00", (float) covMET[0][0]);
    pair.addUserFloat("MEt_cov01", (float) covMET[0][1]);
    pair.addUserFloat("MEt_cov10", (float) covMET[1][0]);
    pair.addUserFloat("MEt_cov11", (float) covMET[1][1]);
    pair.addUserFloat("MEt_significance", significance);

    result->push_back(pair);
  }
  iEvent.put(result);
}

#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(SVfitBypass);

