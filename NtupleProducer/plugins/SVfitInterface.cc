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

  edm::InputTag theCandidateTag;
  std::vector<edm::InputTag> vtheMETTag;
  //int sampleType;
  bool _usePairMET;
  bool _useMVAMET;
};

// ------------------------------------------------------------------


SVfitInterface::SVfitInterface(const edm::ParameterSet& iConfig)
{
  theCandidateTag = iConfig.getParameter<InputTag>("srcPairs");
  _usePairMET = iConfig.getUntrackedParameter<bool>("usePairMET");
  _useMVAMET = iConfig.getUntrackedParameter<bool>("useMVAMET");

  vtheMETTag = iConfig.getParameter<std::vector<edm::InputTag>>("srcMET");

  produces<pat::CompositeCandidateCollection>();
}  



void SVfitInterface::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{  

  // Get lepton pairs
  Handle<View<reco::CompositeCandidate> > pairHandle;
  iEvent.getByLabel(theCandidateTag, pairHandle);
  
  unsigned int elNumber = pairHandle->size();
  unsigned int metNumber = 0;

  /*
  cout << "SVFit DEBUG: input pair collection contains: " << elNumber << " pairs" << endl;
  cout << "SVFit DEBUG: input MEt vector has :" << vtheMETTag.size() << " entries; now printing number of METs in each entry" << endl;
  for (unsigned int i = 0; i < vtheMETTag.size(); i++)
  {
     Handle<View<reco::PFMET> > METHandle_PfMET_DEBUG;
     iEvent.getByLabel(vtheMETTag.at(i), METHandle_PfMET_DEBUG);
     const reco::PFMET& pfMET_DEBUG = (*METHandle_PfMET_DEBUG)[0];
     cout << "  " << i << " | MET num: " << METHandle_PfMET_DEBUG->size() << " | MEt (px, py): " << pfMET_DEBUG.px() <<  " , " << pfMET_DEBUG.py() << endl;
  }
  cout << "FINISHED DEBUG" << endl;
  */

  // MET class type changes if using MVA MEt or 'ordinary' MEt
  
  // create two handles for two types
  Handle<View<reco::PFMET> > METHandle_PfMET;
  Handle<View<pat::MET> >    METHandle_PatMET;
  
  // intialize MET
  double METx = 0.;
  double METy = 0.; 
  TMatrixD covMET(2, 2);
  
  
  // initialize MET once if not using PairMET
  if (!_usePairMET)
  {
      // create handles      
     if (_useMVAMET)
     {
        iEvent.getByLabel(vtheMETTag.at(0), METHandle_PfMET);
        metNumber = METHandle_PfMET->size();
     }
     else
     {
        iEvent.getByLabel(vtheMETTag.at(0), METHandle_PatMET);
        metNumber = METHandle_PatMET->size();
     }     
 
     // guards against wrong number of met values provided
     if (metNumber > 1)     
        edm::LogWarning("SingleMETEntryZeroHasNotSizeOne") << "(SVfitInterface) Warning! Using single MEt, but input MEt collection provided has more than one entry"
                                                           << "   --> will use the first one";
     if (metNumber == 0)
     {
        edm::LogWarning("SingleMETEntryZeroHasNoEntries")  << "(SVfitInterface) Warning! Using single MEt, but input MEt collection provided has no entries"
                                                           << "   --> skipping this event (no output produced)";
        return;
     }

     // set MEt variables
     reco::METCovMatrix covSingleMETbuf;
  
     if (_useMVAMET)
     {   
         const reco::PFMET& pfMET = (*METHandle_PfMET)[0];
         METx = pfMET.px();
         METy = pfMET.py();
         covSingleMETbuf = pfMET.getSignificanceMatrix();
     }
    
     else
     {
        const pat::MET& patMET = (*METHandle_PatMET)[0];
        METx = patMET.px();
        METy = patMET.py();
        covSingleMETbuf = patMET.getSignificanceMatrix();
     }
    
     covMET[0][0] = covSingleMETbuf(0,0);
	 covMET[1][0] = covSingleMETbuf(1,0);
	 covMET[0][1] = covSingleMETbuf(0,1);
	 covMET[1][1] = covSingleMETbuf(1,1);
	
	 // protection against singular matrices
	 if (covMET[0][0] == 0 && covMET[1][0] == 0 && covMET[0][1] == 0 && covMET[1][1] == 0)
	    edm::LogWarning("SingularCovarianceMatrix") << "(SVfitInterface) Warning! Input covariance matrix is singular" 
                                                    << "   --> SVfit algorithm will probably crash...";
  }
  
  // Output collection
  auto_ptr<pat::CompositeCandidateCollection> result( new pat::CompositeCandidateCollection );

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
      if (_useMVAMET)
      {
         iEvent.getByLabel(vtheMETTag.at(i), METHandle_PfMET);
         metNumber = METHandle_PfMET->size();
      }
      else
      {
         iEvent.getByLabel(vtheMETTag.at(i), METHandle_PatMET);   
         metNumber = METHandle_PatMET->size();          
      }
      
      // guards against wrong met elements number
      if (metNumber > 1)     
        edm::LogWarning("SingleMETEntryZeroHasNotSizeOne") << "(SVfitInterface) Warning! Using single MEt, but input MEt collection provided has more than one entry"
                                                           << "   --> will use the first one";
      if (metNumber == 0)
      {
         edm::LogWarning("SingleMETEntryZeroHasNoEntries")  << "(SVfitInterface) Warning! Using single MEt, but input MEt collection provided has no entries"
                                                            << "   --> skipping this event (no output produced)";
         return;
      }
      
      // Get the MET and the covariance matrix
      if (!_useMVAMET) edm::LogWarning("NoMVAMETAsPairRelated") << "(SVfitInterface): Warning! Pair-related MET is available for MVA MET only" 
                                                                << "   --> using MVA MET";
      
    
      const PFMET& pfMET = (*METHandle_PfMET)[0];
      const reco::METCovMatrix& covMETbuf = pfMET.getSignificanceMatrix();
    
      METx = pfMET.px();
      METy = pfMET.py();
      
      covMET[0][0] = covMETbuf(0,0);
      covMET[1][0] = covMETbuf(1,0);
      covMET[0][1] = covMETbuf(0,1);
      covMET[1][1] = covMETbuf(1,1);
      
      // protection against singular matrices
	  if (covMET[0][0] == 0 && covMET[1][0] == 0 && covMET[0][1] == 0 && covMET[1][1] == 0)
	  {
	  	edm::LogWarning("SingularCovarianceMatrix") << "(SVfitInterface) Warning! Input covariance matrix is singular" 
                                                    << "   --> SVfit algorithm will probably crash...";
	  }
    } 
     
    //cout << "l1, l2 pdgId: " << l1->pdgId() << " " << l2->pdgId() << endl;
    
    std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
    measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(l1Type, l1->pt(), l1->eta(), l1->phi(), l1->mass() ));
    measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(l2Type, l2->pt(), l2->eta(), l2->phi(), l2->mass() ));  
     
    // define algorithm (set the debug level to 3 for testing)
    unsigned int verbosity = 0;
    float SVfitMass = -1.;

    SVfitStandaloneAlgorithm algo(measuredTauLeptons, METx, METy, covMET, verbosity);
    algo.addLogM(false); // in general, keep it false when using VEGAS integration
    
    algo.integrateVEGAS();
    
    if ( algo.isValidSolution() )
    {    
      SVfitMass = algo.getMass(); // return value is in units of GeV
    } // otherwise mass will be -1
    
    // add user floats: SVfit mass, met properties, etc..  
    pair.addUserFloat("SVfitMass", SVfitMass);
    pair.addUserFloat("MEt_px", METx);
    pair.addUserFloat("MEt_py", METy);
    
    //pair.addUserFloat("MVAMEt_px", pfMET.px());
    //pair.addUserFloat("MVAMEt_py", pfMET.py());
    //pair.addUserFloat("MVAMEt_pt", pfMET.pt());
    //pair.addUserFloat("MVAMEt_phi", pfMET.phi());
    
    result->push_back(pair);     
  }
  
  //cout << "EXITING PLUGIN..." << endl << endl;
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
