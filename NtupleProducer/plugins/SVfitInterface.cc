/* \class SVfitInterface
**
** This class provides an interface to the SVfit standalone algorithm
** for the computation of the SVfit mass of the lepton pair candidates.
** 
** Input pairs contain a flag describing their decay mode (e, mu, tauh)
** to be used in the algorithm.
**
** input is reco::CompositeCandidate for each lepton
** that is coverted to TLorentzVector to be passed to the algorithm 
**  
** \date:    18 November 2014
** \author:  L. Cadamuro (LLR)
*/
#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/Utilities/interface/InputTag.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/CutSet.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/LeptonIsoHelper.h>
#include <DataFormats/Candidate/interface/ShallowCloneCandidate.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>

#include <TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h>

#include <vector>
#include <string>

using namespace edm;
using namespace std;
using namespace reco;

// ------------------------------------------------------------------

class SVfitInterface : public edm::EDProducer {
 public:
  /// Constructor
  explicit SVfitInterface(const edm::ParameterSet&);
    
  /// Destructor
  ~SVfitInterface(){
  };  

 private:
  virtual void beginJob(){};  
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};
  
  svFitStandalone::kDecayType GetDecayTypeFlag (int pdgId);

  const edm::InputTag theCandidateTag;
  //int sampleType;
  //int setup;
  //const CutSet<pat::Electron> flags;
};

// ------------------------------------------------------------------

// FIX THIS PART --> need to handle MET
SVfitInterface::SVfitInterface(const edm::ParameterSet& iConfig) :
  theCandidateTag(iConfig.getParameter<InputTag>("src"))
//  sampleType(iConfig.getParameter<int>("sampleType")),
//  setup(iConfig.getParameter<int>("setup")),
//  flags(iConfig.getParameter<ParameterSet>("flags")) 
{
  produces<pat::CompositeCandidateCollection>();
}


void SVfitInterface::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{  

  // Get lepton pairs
  Handle<View<reco::CompositeCandidate> > pairHandle;
  iEvent.getByLabel(theCandidateTag, pairHandle);

/*
  InputTag theRhoTag = LeptonIsoHelper::getEleRhoTag(sampleType,setup);
  edm::Handle<double> rhoHandle;
  iEvent.getByLabel(theRhoTag, rhoHandle); 
  double rho = *rhoHandle;

  edm::Handle<vector<Vertex> >  vertexs;
  iEvent.getByLabel("goodPrimaryVertices",vertexs);
*/

  // Output collection
  auto_ptr<pat::CompositeCandidateCollection> result( new pat::CompositeCandidateCollection );

  // loop on all the pairs
  for (unsigned int i = 0; i < pairHandle->size(); ++i)
  {
    // Get the pair and the two leptons composing it
    const CompositeCandidate& pairBuf = (*pairHandle)[i];
    pat::CompositeCandidate pair(pairBuf);
    
    const Candidate *l1 = pair.daughter(0);
    const Candidate *l2 = pair.daughter(1);

    svFitStandalone::kDecayType l1Type = GetDecayTypeFlag (l1->pdgId());
    svFitStandalone::kDecayType l2Type = GetDecayTypeFlag (l2->pdgId());
    
	// convert in TLorentzVector for svFit interface
	svFitStandalone::LorentzVector vecl1( l1->px(), l1->py(), l1->pz(), l1->energy() ); 
	svFitStandalone::LorentzVector vecl2( l2->px(), l2->py(), l2->pz(), l2->energy() ); 
    
    std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
    measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(l1Type, vecl1));
    measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(l2Type, vecl2));

    // define algorithm (set the debug level to 3 for testing)
    unsigned verbosity = 2;
    
    float SVfitMass = -1.;

    SVfitStandaloneAlgorithm algo(measuredTauLeptons, MET, covMET, verbosity);
    algo.addLogM(false);
    
    algo.integrateVEGAS();
    
    if ( algo.isValidSolution() )
    {    
      SVfitMass = algo.getMass(); // return value is in units of GeV
    } // otherwise mass will be -1
      
    pair.addUserFloat("SVfitMass", SVfitMass);
    result->push_back(pair);     
  }
  
  iEvent.put(result);
}

// !!!!!!!
// FIX: error message must be changed properly, no cout

svFitStandalone::kDecayType SVfitInterface::GetDecayTypeFlag (int pdgId)
{
	if (abs(pdgId) == 11) return svFitStandalone::kTauToElecDecay;
	if (abs(pdgId) == 13) return svFitStandalone::kTauToMuDecay;
	if (abs(pdgId) == 15) return svFitStandalone::kTauToHadDecay;
	
	cout << "\n***\n\nWARNING (SVfitInterface): unable to identify decay type from pdgId\n";
	cout << "     ---> Will be treated as an hadronic decay\n\n";
	return svFitStandalone::kTauToHadDecay;
}

#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(SVfitInterface);

