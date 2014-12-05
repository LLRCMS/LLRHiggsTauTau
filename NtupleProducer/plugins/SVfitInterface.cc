/* \class SVfitInterface
**
** This class provides an interface to the SVfit standalone algorithm
** for the computation of the SVfit mass of the lepton pair candidates.
** 
** The decay mode (e, mu, tauh) of each lepton is the pair is asserted
** from the pdgId associated and is used in the algorithm.
**
** input type is reco::CompositeCandidate for each lepton
** that is coverted to TLorentzVector to be passed to the algorithm
**
** output type is pat::CompositeCandidate, i.e. the original pairs
** plus some userfloats containing the SVfit mass and MET px, px, pt, phi. 
**  
** \date:    18 November 2014
** \author:  L. Cadamuro (LLR)
*/

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
  const edm::InputTag theMETTag;
  //int sampleType;
  bool _usePairMET;
};

// ------------------------------------------------------------------


SVfitInterface::SVfitInterface(const edm::ParameterSet& iConfig) :
  theCandidateTag(iConfig.getParameter<InputTag>("srcPairs")),
  theMETTag(iConfig.getParameter<InputTag>("srcMET")),
  _usePairMET(iConfig.getUntrackedParameter<bool>("usePairMET"))
{
  produces<pat::CompositeCandidateCollection>();
}


void SVfitInterface::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{  

  // Get lepton pairs
  Handle<View<reco::CompositeCandidate> > pairHandle;
  iEvent.getByLabel(theCandidateTag, pairHandle);

  // Get MET --  if usePairMet is true, need to handle a reco::PFMET
  //             else need to handle a pat::MET object
  
  // create two handles for two types
  Handle<View<reco::PFMET> > METHandlePairs;
  Handle<View<pat::MET> >    METHandleSingleMET;
  
  if (_usePairMET)   iEvent.getByLabel(theMETTag, METHandlePairs);
  else               iEvent.getByLabel(theMETTag, METHandleSingleMET);
    
  // Output collection
  auto_ptr<pat::CompositeCandidateCollection> result( new pat::CompositeCandidateCollection );



// TO DO:: guards against
// single MET has not size 1! (two cases: 0 --> skip and > 1 --> use as 1)
// single MET has not size 1 but we set usePairMET = false
  // guard against different length collections
  unsigned int elNumber = pairHandle->size();
/*
  if (pairHandle->size() != METHandle->size() )
  {
     edm::LogWarning("InputCollectionsDifferentSize")
	   << "(SVfitInterface): MET and pair collections have different sizes"
       << "    ---> Skipping this event (no output pairs produced)";
     elNumber = 0;
  }
*/
  
  // intialize MET
  svFitStandalone::Vector MET(0., 0., 0.); 
  TMatrixD covMET(2, 2);
  
  if (!_usePairMET)
  {
    const pat::MET& singleMET = (*METHandleSingleMET)[0];
    MET.SetXYZ (singleMET.px(), singleMET.py(), 0.);
    
    const reco::METCovMatrix covSingleMETbuf = singleMET.getSignificanceMatrix();
    
    covMET[0][0] = covSingleMETbuf(0,0);
	covMET[1][0] = covSingleMETbuf(1,0);
	covMET[0][1] = covSingleMETbuf(0,1);
	covMET[1][1] = covSingleMETbuf(1,1);
	
	// protection against singular matrices --> TO FIX BETTER
	if (covMET[0][0] == 0 && covMET[1][0] == 0 && covMET[0][1] == 0 && covMET[1][1] == 0)
	{
		covMET[0][0] = 2.;
		covMET[1][1] = 2.;

		covMET[0][1] = 1.;
		covMET[1][0] = 0.5;

	}
	else cout << "Urrah!!" << endl;
  }
  
  
  // loop on all the pairs
  for (unsigned int i = 0; i < elNumber; ++i)
  {
    // Get the pair and the two leptons composing it
    const CompositeCandidate& pairBuf = (*pairHandle)[i];
    pat::CompositeCandidate pair(pairBuf);
    
    const Candidate *l1 = pair.daughter(0);
    const Candidate *l2 = pair.daughter(1);

    svFitStandalone::kDecayType l1Type = GetDecayTypeFlag (l1->pdgId());
    svFitStandalone::kDecayType l2Type = GetDecayTypeFlag (l2->pdgId());
   
   	if (_usePairMET)
    {
      // Get the MET and the covariance matrix
      const PFMET& pfMET = (*METHandlePairs)[i];
      const reco::METCovMatrix& covMETbuf = pfMET.getSignificanceMatrix();
    
      MET.SetXYZ(pfMET.px(), pfMET.py(), 0.); 
    
      covMET[0][0] = covMETbuf(0,0);
      covMET[1][0] = covMETbuf(1,0);
      covMET[0][1] = covMETbuf(0,1);
      covMET[1][1] = covMETbuf(1,1);
    } 
     
    // convert in TLorentzVector for svFit interface, and into matrix object
    svFitStandalone::LorentzVector vecl1( l1->px(), l1->py(), l1->pz(), l1->energy() ); 
    svFitStandalone::LorentzVector vecl2( l2->px(), l2->py(), l2->pz(), l2->energy() ); 
    
    //cout << "l1, l2 pdgId: " << l1->pdgId() << " " << l2->pdgId() << endl;
    
    // there is a problem with the double hadronic decays ==> I'm skipping them now --> TO FIX
    if (abs(l1->pdgId()) == 15 && abs(l2->pdgId()) == 15) continue;
    
    std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
    measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(l1Type, vecl1));
    measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(l2Type, vecl2));
     
    // define algorithm (set the debug level to 3 for testing)
    unsigned int verbosity = 0;
    float SVfitMass = -1.;

    SVfitStandaloneAlgorithm algo(measuredTauLeptons, MET, covMET, verbosity);
    algo.addLogM(false); // in general, keep it false when using VEGAS integration
    
    algo.integrateVEGAS();
    
    if ( algo.isValidSolution() )
    {    
      SVfitMass = algo.getMass(); // return value is in units of GeV
    } // otherwise mass will be -1
    
    // add user floats: SVfit mass, met properties, etc..  
    pair.addUserFloat("SVfitMass", SVfitMass);
    pair.addUserFloat("MEt_px", MET.X());
    pair.addUserFloat("MEt_py", MET.Y());
    
    //pair.addUserFloat("MVAMEt_px", pfMET.px());
    //pair.addUserFloat("MVAMEt_py", pfMET.py());
    //pair.addUserFloat("MVAMEt_pt", pfMET.pt());
    //pair.addUserFloat("MVAMEt_phi", pfMET.phi());
    
    result->push_back(pair);     
  }
  
  iEvent.put(result);
}


svFitStandalone::kDecayType SVfitInterface::GetDecayTypeFlag (int pdgId)
{
	if (abs(pdgId) == 11) return svFitStandalone::kTauToElecDecay;
	if (abs(pdgId) == 13) return svFitStandalone::kTauToMuDecay;
	if (abs(pdgId) == 15) return svFitStandalone::kTauToHadDecay;
	
	edm::LogWarning("WrongDecayModePdgID")
	   << "(SVfitInterface): unable to identify decay type from pdgId"
	   << "     ---> Decay will be treated as an hadronic decay";
	return svFitStandalone::kTauToHadDecay;
}

#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(SVfitInterface);