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
  
  edm::EDGetTokenT<View<reco::CompositeCandidate> > theCandidateTag;
  // std::vector <edm::EDGetTokenT<View<pat::MET> > > vtheMETTag; // polymorphism of view --> good for reco::PFMET and pat::MET! 
  edm::EDGetTokenT<View<pat::MET>> theMETTag; // polymorphism of view --> good for reco::PFMET and pat::MET! 
  edm::EDGetTokenT<double> theSigTag;
  edm::EDGetTokenT<math::Error<2>::type> theCovTag;

  //int sampleType;
  bool _usePairMET;
  //bool _useMVAMET;
};

// ------------------------------------------------------------------



SVfitBypass::SVfitBypass(const edm::ParameterSet& iConfig):
theCandidateTag(consumes<View<reco::CompositeCandidate> >(iConfig.getParameter<InputTag>("srcPairs"))),
theMETTag(consumes<View<pat::MET>>(iConfig.getParameter<edm::InputTag>("srcMET"))),
theSigTag(consumes<double>(iConfig.getParameter<edm::InputTag>("srcSig"))),
theCovTag(consumes<math::Error<2>::type>(iConfig.getParameter<edm::InputTag>("srcCov")))
{
  //theCandidateTag = iConfig.getParameter<InputTag>("srcPairs");
  //_useMVAMET = iConfig.getUntrackedParameter<bool>("useMVAMET");

  _usePairMET = iConfig.getParameter<bool>("usePairMET");

  // const std::vector<edm::InputTag>& inMET = iConfig.getParameter<std::vector<edm::InputTag> >("srcMET");
  // for (std::vector<edm::InputTag>::const_iterator it = inMET.begin(); it != inMET.end(); ++it)
  // {      
  //   // vtheMETTag.emplace_back(consumes<edm::View<reco::MET> >(*it) );
  //   vtheMETTag.emplace_back(consumes<edm::View<pat::MET> >(*it) );
  // }

  produces<pat::CompositeCandidateCollection>();
}  



void SVfitBypass::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{  

  // Get lepton pairs
  Handle<View<reco::CompositeCandidate> > pairHandle;
  iEvent.getByToken(theCandidateTag, pairHandle);
  
  unsigned int elNumber = pairHandle->size();
  
    // Output collection
  //auto_ptr<pat::CompositeCandidateCollection> result( new pat::CompositeCandidateCollection );
  std::unique_ptr<pat::CompositeCandidateCollection> result( new pat::CompositeCandidateCollection );

  // get event pat MET to be saved in output
  double METx = 0.;
  double METy = 0.; 
  double uncorrMETx = -999.;
  double uncorrMETy = -999.; 
  TMatrixD covMET(2, 2);
  float significance = -999.;
  Handle<View<pat::MET> > METHandle;
  
  if (!_usePairMET)
  {
    // iEvent.getByToken(vtheMETTag.at(0), METHandle);
    iEvent.getByToken(theMETTag, METHandle);
    const pat::MET& patMET = (*METHandle)[0];
    METx = patMET.px();
    METy = patMET.py();
    Handle<double> significanceHandle;
    Handle<math::Error<2>::type> covHandle;
    iEvent.getByToken (theSigTag, significanceHandle);
    iEvent.getByToken (theCovTag, covHandle);
    covMET[0][0] = (*covHandle)(0,0);
    covMET[1][0] = (*covHandle)(1,0);
    covMET[0][1] = covMET[1][0]; // (1,0) is the only one saved
    covMET[1][1] = (*covHandle)(1,1);
    significance = (float) (*significanceHandle);
  }

  // loop on all the pairs
  for (unsigned int i = 0; i < elNumber; ++i)
  {
    if (_usePairMET)
    {
      // iEvent.getByToken(vtheMETTag.at(i), METHandle);
      iEvent.getByToken(theMETTag, METHandle);
      //metNumber = METHandle->size();

      // const PFMET* pfMET = (PFMET*) &((*METHandle)[0]) ; // all this to transform the type of the pointer!
      // const pat::MET* patMET = &((*METHandle)[0]);
      const pat::MET* patMET = &((*METHandle)[i]);
      const reco::METCovMatrix& covMETbuf = patMET->getSignificanceMatrix();
      significance = (float) patMET->significance();

      METx = patMET->px();
      METy = patMET->py();

      uncorrMETx = ( patMET->hasUserFloat("uncorrPx") ) ? patMET->userFloat("uncorrPx") : -999;
      uncorrMETy = ( patMET->hasUserFloat("uncorrPy") ) ? patMET->userFloat("uncorrPy") : -999;

      covMET[0][0] = covMETbuf(0,0);
      covMET[1][0] = covMETbuf(1,0);
      covMET[0][1] = covMETbuf(0,1);
      covMET[1][1] = covMETbuf(1,1);
    }

    // Get the pair and the two leptons composing it
    const CompositeCandidate& pairBuf = (*pairHandle)[i];
    pat::CompositeCandidate pair(pairBuf);

    float SVfitMass = -999.;
    // add user floats: SVfit mass, met properties, etc..  
    pair.addUserFloat("SVfitMass", SVfitMass);
    pair.addUserFloat("SVfitTransverseMass", -999);
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
    pair.addUserFloat("uncorrMEt_px", (float) uncorrMETx);
    pair.addUserFloat("uncorrMEt_py", (float) uncorrMETy);
    pair.addUserFloat("MEt_cov00", (float) covMET[0][0]);
    pair.addUserFloat("MEt_cov01", (float) covMET[0][1]);
    pair.addUserFloat("MEt_cov10", (float) covMET[1][0]);
    pair.addUserFloat("MEt_cov11", (float) covMET[1][1]);
    pair.addUserFloat("MEt_significance", significance);

    result->push_back(pair);
  }
  //iEvent.put(result);
  iEvent.put(std::move(result));
}

#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(SVfitBypass);

