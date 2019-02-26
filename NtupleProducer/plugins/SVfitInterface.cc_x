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
** \author:  L. Cadamuro, O. Davignon (LLR)
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
#include <LLRHiggsTauTau/NtupleProducer/interface/DaughterDataHelpers.h>
#include <TLorentzVector.h>

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
  ~SVfitInterface();  

 private:
  virtual void beginJob(){};  
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};

  svFitStandalone::kDecayType GetDecayTypeFlag (int pdgId);
  bool Switch (svFitStandalone::kDecayType type1, double pt1, svFitStandalone::kDecayType type2, double pt2);
  double GetMass (svFitStandalone::kDecayType type, double candMass);
  bool IsInteresting (const reco::Candidate *l1, const reco::Candidate *l2); // if true, compute SVFit

  edm::EDGetTokenT<View<reco::CompositeCandidate> > theCandidateTag;
  // std::vector <edm::EDGetTokenT<View<pat::MET> > > vtheMETTag; // polymorphism of view --> good for reco::PFMET and pat::MET! 
  edm::EDGetTokenT<View<pat::MET> > theMETTag;
  edm::EDGetTokenT<double> theSigTag;
  edm::EDGetTokenT<math::Error<2>::type> theCovTag;
  bool _usePairMET;
  bool _computeForUpDownTES;
  TFile* inputFile_visPtResolution_;
  
  // 6,7,8 are expected to be unused
  enum pairType {
    kMuHad  = 0,
    kEHad   = 1,
    kHadHad = 2,
    kMuMu   = 3,
    kEE     = 4,
    kEMu    = 5,
    kEEPrompt = 6, // prompt Z->ee/mumu decays
    kMuMuPrompt = 7,
    kOther  = 8 // for e.g. h->bb
  };

};

// ------------------------------------------------------------------


SVfitInterface::SVfitInterface(const edm::ParameterSet& iConfig):
theCandidateTag(consumes<View<reco::CompositeCandidate> >(iConfig.getParameter<InputTag>("srcPairs"))),
//vtheMETTag(consumes<View<pat::MET>>(iConfig.getParameter<edm::InputTag>("srcMET"))),
theMETTag(consumes<View<pat::MET>>(iConfig.getParameter<edm::InputTag>("srcMET"))),
theSigTag(consumes<double>(iConfig.getParameter<edm::InputTag>("srcSig"))),
theCovTag(consumes<math::Error<2>::type>(iConfig.getParameter<edm::InputTag>("srcCov")))
{
  //theCandidateTag = iConfig.getParameter<InputTag>("srcPairs");
  _usePairMET = iConfig.getParameter<bool>("usePairMET");
  _computeForUpDownTES = iConfig.getParameter<bool>("computeForUpDownTES");
  
  // const std::vector<edm::InputTag>& inMET = iConfig.getParameter<std::vector<edm::InputTag> >("srcMET");
  // for (std::vector<edm::InputTag>::const_iterator it = inMET.begin(); it != inMET.end(); ++it)
  // {      
  //   // vtheMETTag.emplace_back(consumes<edm::View<reco::MET> >(*it) );
  //   vtheMETTag.emplace_back(consumes<edm::View<pat::MET> >(*it) );
  // }

  edm::FileInPath inputFileName_visPtResolution("TauAnalysis/SVfitStandalone/data/svFitVisMassAndPtResolutionPDF.root");
  TH1::AddDirectory(false);  
  inputFile_visPtResolution_ = new TFile(inputFileName_visPtResolution.fullPath().data());

  produces<pat::CompositeCandidateCollection>();
  
}  

SVfitInterface::~SVfitInterface()
{
    delete inputFile_visPtResolution_;
}

void SVfitInterface::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{  

  // Get lepton pairs
  Handle<View<reco::CompositeCandidate> > pairHandle;
  iEvent.getByToken(theCandidateTag, pairHandle);
  
  unsigned int pairNumber = pairHandle->size();
  unsigned int metNumber = 0;


  // MET class type changes if using MVA MEt or 'ordinary' MEt
  
  // create handle -- view makes possible to use base class type reco::MET
  // but now everything is produced as pat::MET
  Handle<View<pat::MET> > METHandle;

  // Handle<View<reco::PFMET> > METHandle_PfMET;
  // Handle<View<pat::MET> >    METHandle_PatMET;
  
  // intialize MET
  double METx = 0.;
  double METy = 0.; 
  double uncorrMETx = -999.;
  double uncorrMETy = -999.; 
  TMatrixD covMET(2, 2);
  float significance = -999.;

  iEvent.getByToken(theMETTag, METHandle);    
  // initialize MET once if not using PairMET
  if (!_usePairMET)
  {   
     metNumber = METHandle->size();
     if (metNumber != 1)     
        edm::LogWarning("pfMetHasNotSizeOne") << "(SVfitInterface) Warning! Using single pf MEt, but input MEt collection size is different from 1"
                                                           << "   --> using MET entry num. 0";
     const pat::MET& patMET = (*METHandle)[0];
     METx = patMET.px();
     METy = patMET.py();

     Handle<double> significanceHandle;
     Handle<math::Error<2>::type> covHandle;
     
     iEvent.getByToken (theSigTag, significanceHandle);
     iEvent.getByToken (theCovTag, covHandle);
     
     //cout << *significanceHandle << " " << *sig00Handle << " " << *sig01Handle << " " << *sig10Handle << " " << *sig11Handle << endl;
     covMET[0][0] = (*covHandle)(0,0);
     covMET[1][0] = (*covHandle)(1,0);
     covMET[0][1] = covMET[1][0]; // (1,0) is the only one saved
     covMET[1][1] = (*covHandle)(1,1);

     significance = (float) (*significanceHandle);
     
     // protection against singular matrices
     if (covMET[0][0] == 0 && covMET[1][0] == 0 && covMET[0][1] == 0 && covMET[1][1] == 0)
        edm::LogWarning("SingularCovarianceMatrix") << "(SVfitInterface) Warning! Input covariance matrix is singular" 
                                                    << "   --> SVfit algorithm will probably crash...";
  }
  
  // Output collection
  //auto_ptr<pat::CompositeCandidateCollection> result( new pat::CompositeCandidateCollection );
  std::unique_ptr<pat::CompositeCandidateCollection> result( new pat::CompositeCandidateCollection );

  // loop on all the pairs
  for (unsigned int i = 0; i < pairNumber; ++i)
  {    
    // Get the pair and the two leptons composing it
    const CompositeCandidate& pairBuf = (*pairHandle)[i];
    pat::CompositeCandidate pair(pairBuf);
    
    const Candidate *l1 = pair.daughter(0);
    const Candidate *l2 = pair.daughter(1);

    svFitStandalone::kDecayType l1Type = GetDecayTypeFlag (l1->pdgId());
    svFitStandalone::kDecayType l2Type = GetDecayTypeFlag (l2->pdgId());
    double mass1 = GetMass (l1Type, l1->mass());
    double mass2 = GetMass (l2Type, l2->mass());
   
    int decay1 = -1;
    int decay2 = -1;
    if (l1Type == svFitStandalone::kTauToHadDecay) decay1 = (int)(userdatahelpers::getUserFloat(l1,"decayMode"));
    if (l2Type == svFitStandalone::kTauToHadDecay) decay2 = (int)(userdatahelpers::getUserFloat(l2,"decayMode"));
   
    // // compute SVfit only if tau(s) pass the OldDM discriminator (avoids crashes...)
    // bool passDM = true;
    // if (l1Type == svFitStandalone::kTauToHadDecay && userdatahelpers::getUserInt(l1,"decayModeFindingNewDMs") != 1) passDM = false;
    // if (l2Type == svFitStandalone::kTauToHadDecay && userdatahelpers::getUserInt(l2,"decayModeFindingNewDMs") != 1) passDM = false;
    // DEBUG
    //if (l1Type == svFitStandalone::kTauToHadDecay && userdatahelpers::getUserInt(l1,"decayModeFindingNewDMs") == 1 && userdatahelpers::getUserInt(l1,"decayModeFindingOldDMs") != 1) cout << "Only new DM for tau 1!!" << endl;
    //if (l2Type == svFitStandalone::kTauToHadDecay && userdatahelpers::getUserInt(l2,"decayModeFindingNewDMs") == 1 && userdatahelpers::getUserInt(l2,"decayModeFindingOldDMs") != 1) cout << "Only new DM for tau 2!!" << endl;
    
    // // do not compute SVfit on ee or mumu pairs
    // bool isGoodPairType = true;
    // if (l1Type == svFitStandalone::kTauToElecDecay && l2Type == svFitStandalone::kTauToElecDecay) isGoodPairType = false;
    // if (l1Type == svFitStandalone::kTauToMuDecay && l2Type == svFitStandalone::kTauToMuDecay) isGoodPairType = false;

    // do not compute SVfit if the two pairs are too close one to the other
    // bool isGoodDR = true;
    // if (deltaR(l1->p4(), l2->p4()) < 0.25) isGoodDR = false;
    
    bool swi = Switch (l1Type, l1->pt(), l2Type, l2->pt());
  
    if (_usePairMET)
    {
      // iEvent.getByToken(theMETTag, METHandle);
      metNumber = METHandle->size();

      // const PFMET* pfMET = (PFMET*) &((*METHandle)[0]) ; // all this to transform the type of the pointer!
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

      // protection against singular matrices
      if (covMET[0][0] == 0 && covMET[1][0] == 0 && covMET[0][1] == 0 && covMET[1][1] == 0)
      {
          edm::LogWarning("SingularCovarianceMatrix") << "(SVfitInterface) Warning! Input covariance matrix is singular" 
                                                    << "   --> SVfit algorithm will probably crash...";
      }
    } 

    // prepare tau nominal, up, down candidates            
    std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
    std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptonsTauUp;
    std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptonsTauDown;

    // init shifted version to nominal
    TLorentzVector l1_UP (l1->px(), l1->py(), l1->pz(), l1->energy());
    TLorentzVector l1_DOWN = l1_UP;
    TLorentzVector l2_UP (l2->px(), l2->py(), l2->pz(), l2->energy());
    TLorentzVector l2_DOWN = l2_UP;

    bool l1shifted = false;
    bool l2shifted = false;
    if (userdatahelpers::hasUserInt(l1,"TauUpExists"))
      l1shifted = (userdatahelpers::getUserInt(l1,"TauUpExists") == 1 ? true : false);
    if (userdatahelpers::hasUserInt(l2,"TauUpExists"))
      l2shifted = (userdatahelpers::getUserInt(l2,"TauUpExists") == 1 ? true : false);

    if (l1shifted)
    {
      float pxUp = userdatahelpers::getUserFloat(l1,"px_TauUp");
      float pyUp = userdatahelpers::getUserFloat(l1,"py_TauUp");
      float pzUp = userdatahelpers::getUserFloat(l1,"pz_TauUp");
      float eUp  = userdatahelpers::getUserFloat(l1,"e_TauUp");
      l1_UP.SetPxPyPzE (pxUp, pyUp, pzUp, eUp);

      float pxDown = userdatahelpers::getUserFloat(l1,"px_TauDown");
      float pyDown = userdatahelpers::getUserFloat(l1,"py_TauDown");
      float pzDown = userdatahelpers::getUserFloat(l1,"pz_TauDown");
      float eDown  = userdatahelpers::getUserFloat(l1,"e_TauDown");
      l1_DOWN.SetPxPyPzE (pxDown, pyDown, pzDown, eDown);
    }

    if (l2shifted)
    {
      float pxUp = userdatahelpers::getUserFloat(l2,"px_TauUp");
      float pyUp = userdatahelpers::getUserFloat(l2,"py_TauUp");
      float pzUp = userdatahelpers::getUserFloat(l2,"pz_TauUp");
      float eUp  = userdatahelpers::getUserFloat(l2,"e_TauUp");
      l2_UP.SetPxPyPzE (pxUp, pyUp, pzUp, eUp);

      float pxDown = userdatahelpers::getUserFloat(l2,"px_TauDown");
      float pyDown = userdatahelpers::getUserFloat(l2,"py_TauDown");
      float pzDown = userdatahelpers::getUserFloat(l2,"pz_TauDown");
      float eDown  = userdatahelpers::getUserFloat(l2,"e_TauDown");
      l2_DOWN.SetPxPyPzE (pxDown, pyDown, pzDown, eDown);
    }

    // set lepton vector, ordered for SVfit
    if (swi)  // 2 first, 1 second (switch)
    {
      measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(l2Type, l2->pt(), l2->eta(), l2->phi(), mass2, decay2 ));  
      measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(l1Type, l1->pt(), l1->eta(), l1->phi(), mass1, decay1 ));

      measuredTauLeptonsTauUp.push_back(svFitStandalone::MeasuredTauLepton(l2Type, l2_UP.Pt(), l2_UP.Eta(), l2_UP.Phi(), mass2, decay2 ));  
      measuredTauLeptonsTauUp.push_back(svFitStandalone::MeasuredTauLepton(l1Type, l1_UP.Pt(), l1_UP.Eta(), l1_UP.Phi(), mass1, decay1 ));

      measuredTauLeptonsTauDown.push_back(svFitStandalone::MeasuredTauLepton(l2Type, l2_DOWN.Pt(), l2_DOWN.Eta(), l2_DOWN.Phi(), mass2, decay2 ));  
      measuredTauLeptonsTauDown.push_back(svFitStandalone::MeasuredTauLepton(l1Type, l1_DOWN.Pt(), l1_DOWN.Eta(), l1_DOWN.Phi(), mass1, decay1 ));
    }

    else
    {
      measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(l1Type, l1->pt(), l1->eta(), l1->phi(), mass1, decay1 ));
      measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(l2Type, l2->pt(), l2->eta(), l2->phi(), mass2, decay2 ));  

      measuredTauLeptonsTauUp.push_back(svFitStandalone::MeasuredTauLepton(l1Type, l1_UP.Pt(), l1_UP.Eta(), l1_UP.Phi(), mass1, decay1 ));
      measuredTauLeptonsTauUp.push_back(svFitStandalone::MeasuredTauLepton(l2Type, l2_UP.Pt(), l2_UP.Eta(), l2_UP.Phi(), mass2, decay2 ));  

      measuredTauLeptonsTauDown.push_back(svFitStandalone::MeasuredTauLepton(l1Type, l1_DOWN.Pt(), l1_DOWN.Eta(), l1_DOWN.Phi(), mass1, decay1 ));
      measuredTauLeptonsTauDown.push_back(svFitStandalone::MeasuredTauLepton(l2Type, l2_DOWN.Pt(), l2_DOWN.Eta(), l2_DOWN.Phi(), mass2, decay2 ));        
    }

    // define algorithm (set the debug level to 3 for testing)
    unsigned int verbosity = 0;
    
    double SVfitMass = -999.;
    double SVfitMassTauUp = -999.;
    double SVfitMassTauDown = -999.;

    double SVfitMassUnc = -999.;
    double SVfitMassUncTauUp = -999.;
    double SVfitMassUncTauDown = -999.;

    double SVfitTransverseMass = -999.;
    double SVfitTransverseMassTauUp = -999.;
    double SVfitTransverseMassTauDown = -999.;

    double SVfitTransverseMassUnc = -999.;
    double SVfitTransverseMassUncTauUp = -999.;
    double SVfitTransverseMassUncTauDown = -999.;

    double SVpt = -999.;
    double SVptTauUp = -999.;
    double SVptTauDown = -999.;

    double SVeta = -999.;
    double SVetaTauUp = -999.;
    double SVetaTauDown = -999.;

    double SVphi = -999.;
    double SVphiTauUp = -999.;
    double SVphiTauDown = -999.;

    double SVptUnc = -999.;
    double SVptUncTauUp = -999.;
    double SVptUncTauDown = -999.;

    double SVetaUnc = -999.;
    double SVetaUncTauUp = -999.;
    double SVetaUncTauDown = -999.;

    double SVphiUnc = -999.;
    double SVphiUncTauUp = -999.;
    double SVphiUncTauDown = -999.;

    double SVMETRho = -999.; // fitted MET
    double SVMETRhoTauUp = -999.; // fitted MET
    double SVMETRhoTauDown = -999.; // fitted MET

    double SVMETPhi = -999.;
    double SVMETPhiTauUp = -999.;
    double SVMETPhiTauDown = -999.;

    // minimal acceptance on the *nominal* pair in order to process SVFit
    // Bool_t GoodPairFlag = kTRUE;

    // Bool_t isMuTau = kFALSE;
    // Bool_t isTauMu = kFALSE;
    // Bool_t isElecTau = kFALSE;
    // Bool_t isTauElec = kFALSE;
    // Bool_t isTauTau = kFALSE;
    // Bool_t isElecMu = kFALSE;
    // Bool_t isMuElec = kFALSE;

    // if(l1Type == svFitStandalone::kTauToMuDecay && l2Type == svFitStandalone::kTauToHadDecay) isMuTau = kTRUE;
    // else if(l2Type == svFitStandalone::kTauToMuDecay && l1Type == svFitStandalone::kTauToHadDecay) isTauMu = kTRUE;
    // else if(l1Type == svFitStandalone::kTauToElecDecay && l2Type == svFitStandalone::kTauToHadDecay) isElecTau = kTRUE;
    // else if(l2Type == svFitStandalone::kTauToElecDecay && l1Type == svFitStandalone::kTauToHadDecay) isTauElec = kTRUE;
    // else if(l1Type == svFitStandalone::kTauToHadDecay && l2Type == svFitStandalone::kTauToHadDecay) isTauTau = kTRUE;
    // else if(l1Type == svFitStandalone::kTauToElecDecay && l2Type == svFitStandalone::kTauToMuDecay) isElecMu = kTRUE;
    // else if(l2Type == svFitStandalone::kTauToElecDecay && l1Type == svFitStandalone::kTauToMuDecay) isMuElec = kTRUE;
    // else GoodPairFlag = kFALSE;

    // //MuTau or TauMu
    // //mu leg
    // if(isMuTau && l1->pt()< 17.) GoodPairFlag = kFALSE;
    // if(isTauMu && l2->pt()< 17.) GoodPairFlag = kFALSE;
    // //tau leg
    // if(isMuTau && l2->pt()< 19.) GoodPairFlag = kFALSE;
    // if(isTauMu && l1->pt()< 19.) GoodPairFlag = kFALSE;

    // //ElecTau or TauElec
    // //elec leg
    // if(isElecTau && l1->pt()< 17.) GoodPairFlag = kFALSE;
    // if(isTauElec && l2->pt()< 17.) GoodPairFlag = kFALSE;
    // //tau leg
    // if(isElecTau && l2->pt()< 19.) GoodPairFlag = kFALSE;
    // if(isTauElec && l1->pt()< 19.) GoodPairFlag = kFALSE;

    // //Cuts to be applied on any lepton which is a tau
    // if(l1Type==svFitStandalone::kTauToHadDecay && userdatahelpers::getUserFloat(l1,"byCombinedIsolationDeltaBetaCorrRaw3Hits")>7.)  GoodPairFlag = kFALSE;
    // if(l1Type==svFitStandalone::kTauToHadDecay && userdatahelpers::getUserInt(l1,"decayModeFindingNewDMs") != 1)  GoodPairFlag = kFALSE;
    // if(l2Type==svFitStandalone::kTauToHadDecay && userdatahelpers::getUserFloat(l2,"byCombinedIsolationDeltaBetaCorrRaw3Hits")>7.)  GoodPairFlag = kFALSE;
    // if(l2Type==svFitStandalone::kTauToHadDecay && userdatahelpers::getUserInt(l2,"decayModeFindingNewDMs") != 1)  GoodPairFlag = kFALSE;

    // //TauTau
    // if(isTauTau && (l1->pt()< 30. || l2->pt()< 30.)) GoodPairFlag = kFALSE;

    // //ElecMu or MuElec
    // //elec leg
    // if(isElecMu && l1->pt()< 13.) GoodPairFlag = kFALSE;
    // if(isElecMu && l2->pt()< 10.) GoodPairFlag = kFALSE;
    // //mu leg
    // if(isMuElec && l2->pt()< 13.) GoodPairFlag = kFALSE;
    // if(isMuElec && l1->pt()< 10.) GoodPairFlag = kFALSE;    

    // only run SVfit if taus are passing discriminator, skip mumu and ee pairs, apply very loose quality cuts on objects
    // if (isGoodDR && GoodPairFlag)
    if (IsInteresting(l1, l2))
    {
      SVfitStandaloneAlgorithm algo(measuredTauLeptons, METx, METy, covMET, verbosity);
      algo.addLogM(false); // in general, keep it false when using VEGAS integration
      algo.shiftVisPt(true, inputFile_visPtResolution_);
      // cout << "Integrating SVFit nominal for tau " << i << " --> pt=" << l1->pt() << " eta=" << l1->eta() << " phi=" << l1->phi() << " type=" << l1Type << " -- pt=" << l2->pt() << "eta= " << l2->eta() << " phi=" << l2->phi() << " type=" << l2Type << endl;
      algo.integrateMarkovChain();
      
      if ( algo.isValidSolution() )
      {
        SVfitMass = static_cast<svFitStandalone::MCPtEtaPhiMassAdapter*>(algo.getMCQuantitiesAdapter())->getMass(); // full mass of tau lepton pair in units of GeV
        SVfitTransverseMass = static_cast<svFitStandalone::MCPtEtaPhiMassAdapter*>(algo.getMCQuantitiesAdapter())->getTransverseMass();
        SVpt      = static_cast<svFitStandalone::MCPtEtaPhiMassAdapter*>(algo.getMCQuantitiesAdapter())->getPt();
        SVeta     = static_cast<svFitStandalone::MCPtEtaPhiMassAdapter*>(algo.getMCQuantitiesAdapter())->getEta();
        SVphi     = static_cast<svFitStandalone::MCPtEtaPhiMassAdapter*>(algo.getMCQuantitiesAdapter())->getPhi();
        SVfitMassUnc = static_cast<svFitStandalone::MCPtEtaPhiMassAdapter*>(algo.getMCQuantitiesAdapter())->getMassUncert();
        SVfitTransverseMassUnc = static_cast<svFitStandalone::MCPtEtaPhiMassAdapter*>(algo.getMCQuantitiesAdapter())->getTransverseMassUncert();
        SVptUnc   = static_cast<svFitStandalone::MCPtEtaPhiMassAdapter*>(algo.getMCQuantitiesAdapter())->getPtUncert();
        SVetaUnc  = static_cast<svFitStandalone::MCPtEtaPhiMassAdapter*>(algo.getMCQuantitiesAdapter())->getEtaUncert();
        SVphiUnc  = static_cast<svFitStandalone::MCPtEtaPhiMassAdapter*>(algo.getMCQuantitiesAdapter())->getPhiUncert();
        ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > fittedDiTauSystem(SVpt, SVeta, SVphi, SVfitMass);
        Vector fittedMET = (fittedDiTauSystem.Vect() - (algo.measuredDiTauSystem()).Vect());
        SVMETRho = fittedMET.Rho();
        SVMETPhi = fittedMET.Phi();
      }
      else
        SVfitMass = -111; // -111: SVfit failed (cfr: -999: SVfit not computed)

      // compute up/down SVfit variations
      if ( (l1shifted || l2shifted) && _computeForUpDownTES)
      {
        
        // UP
        // cout << "Integrating SVFit UP for tau " << i << " --> " << l1_UP.Pt() << " " << l1_UP.Eta() << " " << l1Type << " -- " << l2_UP.Pt() << " " << l2_UP.Eta() << " " << l2Type << endl;
        SVfitStandaloneAlgorithm algoTauUp(measuredTauLeptonsTauUp, METx, METy, covMET, verbosity);
        algoTauUp.addLogM(false); // in general, keep it false when using VEGAS integration
        algoTauUp.shiftVisPt(true, inputFile_visPtResolution_);
        algoTauUp.integrateMarkovChain();
        
        if ( algoTauUp.isValidSolution() )
        {    
          SVfitMassTauUp = static_cast<svFitStandalone::MCPtEtaPhiMassAdapter*>(algoTauUp.getMCQuantitiesAdapter())->getMass(); // full mass of tau lepton pair in units of GeV
          SVfitTransverseMassTauUp = static_cast<svFitStandalone::MCPtEtaPhiMassAdapter*>(algoTauUp.getMCQuantitiesAdapter())->getTransverseMass();
          SVptTauUp      = static_cast<svFitStandalone::MCPtEtaPhiMassAdapter*>(algoTauUp.getMCQuantitiesAdapter())->getPt();
          SVetaTauUp     = static_cast<svFitStandalone::MCPtEtaPhiMassAdapter*>(algoTauUp.getMCQuantitiesAdapter())->getEta();
          SVphiTauUp     = static_cast<svFitStandalone::MCPtEtaPhiMassAdapter*>(algoTauUp.getMCQuantitiesAdapter())->getPhi();
          SVfitMassUncTauUp = static_cast<svFitStandalone::MCPtEtaPhiMassAdapter*>(algoTauUp.getMCQuantitiesAdapter())->getMassUncert();
          SVfitTransverseMassUncTauUp = static_cast<svFitStandalone::MCPtEtaPhiMassAdapter*>(algoTauUp.getMCQuantitiesAdapter())->getTransverseMassUncert();
          SVptUncTauUp   = static_cast<svFitStandalone::MCPtEtaPhiMassAdapter*>(algoTauUp.getMCQuantitiesAdapter())->getPtUncert();
          SVetaUncTauUp  = static_cast<svFitStandalone::MCPtEtaPhiMassAdapter*>(algoTauUp.getMCQuantitiesAdapter())->getEtaUncert();
          SVphiUncTauUp  = static_cast<svFitStandalone::MCPtEtaPhiMassAdapter*>(algoTauUp.getMCQuantitiesAdapter())->getPhiUncert();

          ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > fittedDiTauSystemUp(SVptTauUp, SVetaTauUp, SVphiTauUp, SVfitMassTauUp);
          Vector fittedMETUp = (fittedDiTauSystemUp.Vect() - (algoTauUp.measuredDiTauSystem()).Vect());
          SVMETRhoTauUp = fittedMETUp.Rho();
          SVMETPhiTauUp = fittedMETUp.Phi();
        }
        else
          SVfitMassTauUp = -111; // -111: SVfit failed (cfr: -999: SVfit not computed)

        // DOWN
        // cout << "Integrating SVFit DOWN for tau " << i << " --> " << l1_DOWN.Pt() << " " << l1_DOWN.Eta() << " " << l1Type << " -- " << l2_DOWN.Pt() << " " << l2_DOWN.Eta() << " " << l2Type << endl;
        SVfitStandaloneAlgorithm algoTauDown(measuredTauLeptonsTauDown, METx, METy, covMET, verbosity);
        algoTauDown.addLogM(false); // in general, keep it false when using VEGAS integration
        algoTauDown.shiftVisPt(true, inputFile_visPtResolution_);
        algoTauDown.integrateMarkovChain();

        if ( algoTauDown.isValidSolution() )
        {    
          SVfitMassTauDown = static_cast<svFitStandalone::MCPtEtaPhiMassAdapter*>(algoTauDown.getMCQuantitiesAdapter())->getMass(); // full mass of tau lepton pair in units of GeV
          SVfitTransverseMassTauDown = static_cast<svFitStandalone::MCPtEtaPhiMassAdapter*>(algoTauDown.getMCQuantitiesAdapter())->getTransverseMass();
          SVptTauDown      = static_cast<svFitStandalone::MCPtEtaPhiMassAdapter*>(algoTauDown.getMCQuantitiesAdapter())->getPt();
          SVetaTauDown     = static_cast<svFitStandalone::MCPtEtaPhiMassAdapter*>(algoTauDown.getMCQuantitiesAdapter())->getEta();
          SVphiTauDown     = static_cast<svFitStandalone::MCPtEtaPhiMassAdapter*>(algoTauDown.getMCQuantitiesAdapter())->getPhi();
          SVfitMassUncTauDown = static_cast<svFitStandalone::MCPtEtaPhiMassAdapter*>(algoTauDown.getMCQuantitiesAdapter())->getMassUncert();
          SVfitTransverseMassUncTauDown = static_cast<svFitStandalone::MCPtEtaPhiMassAdapter*>(algoTauDown.getMCQuantitiesAdapter())->getTransverseMassUncert();
          SVptUncTauDown   = static_cast<svFitStandalone::MCPtEtaPhiMassAdapter*>(algoTauDown.getMCQuantitiesAdapter())->getPtUncert();
          SVetaUncTauDown  = static_cast<svFitStandalone::MCPtEtaPhiMassAdapter*>(algoTauDown.getMCQuantitiesAdapter())->getEtaUncert();
          SVphiUncTauDown  = static_cast<svFitStandalone::MCPtEtaPhiMassAdapter*>(algoTauDown.getMCQuantitiesAdapter())->getPhiUncert();

          ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > fittedDiTauSystemDown(SVptTauDown, SVetaTauDown, SVphiTauDown, SVfitMassTauDown);
          Vector fittedMETDown = (fittedDiTauSystemDown.Vect() - (algoTauDown.measuredDiTauSystem()).Vect());
          SVMETRhoTauDown = fittedMETDown.Rho();
          SVMETPhiTauDown = fittedMETDown.Phi();
        }
        else
          SVfitMassTauDown = -111; // -111: SVfit failed (cfr: -999: SVfit not computed)
      }
      else if (_computeForUpDownTES) // if I asked to have UP/DOWN variation, but this pair has not tau shifted, simply put central value. 
      {                              // instead, if I dindn't ask for up/down, I get -999 everywhere to remember my mistakes
          SVfitMassTauUp = SVfitMassTauDown = SVfitMass ;
          SVfitTransverseMassTauUp = SVfitTransverseMassTauDown = SVfitTransverseMass ;
          SVptTauUp = SVptTauDown = SVpt ;
          SVetaTauUp = SVetaTauDown = SVeta ;
          SVphiTauUp = SVphiTauDown = SVphi ;
          SVfitMassUncTauUp = SVfitMassUncTauDown = SVfitMassUnc ;
          SVfitTransverseMassUncTauUp = SVfitTransverseMassUncTauDown = SVfitTransverseMassUnc ;
          SVptUncTauUp = SVptUncTauDown = SVptUnc ;
          SVetaUncTauUp = SVetaUncTauDown = SVetaUnc ;
          SVphiUncTauUp = SVphiUncTauDown = SVphiUnc ;
          SVMETRhoTauUp = SVMETRhoTauDown = SVMETRho ;
          SVMETPhiTauUp = SVMETPhiTauDown = SVMETPhi ;
      }
    } // end of quality checks IF
    
    // add user floats: SVfit mass, met properties, etc..  
    pair.addUserFloat("SVfitMass", (float) SVfitMass);
    pair.addUserFloat("SVfitMassTauUp", (float) SVfitMassTauUp);
    pair.addUserFloat("SVfitMassTauDown", (float) SVfitMassTauDown);

    pair.addUserFloat("SVfitTransverseMass", (float) SVfitTransverseMass);
    pair.addUserFloat("SVfitTransverseMassTauUp", (float) SVfitTransverseMassTauUp);
    pair.addUserFloat("SVfitTransverseMassTauDown", (float) SVfitTransverseMassTauDown);

    pair.addUserFloat("SVfit_pt", (float) SVpt);
    pair.addUserFloat("SVfit_ptTauUp", (float) SVptTauUp);
    pair.addUserFloat("SVfit_ptTauDown", (float) SVptTauDown);

    pair.addUserFloat("SVfit_eta", (float) SVeta);
    pair.addUserFloat("SVfit_etaTauUp", (float) SVetaTauUp);
    pair.addUserFloat("SVfit_etaTauDown", (float) SVetaTauDown);

    pair.addUserFloat("SVfit_phi", (float) SVphi);
    pair.addUserFloat("SVfit_phiTauUp", (float) SVphiTauUp);
    pair.addUserFloat("SVfit_phiTauDown", (float) SVphiTauDown);

    pair.addUserFloat("SVfitMassUnc", (float) SVfitMassUnc);
    pair.addUserFloat("SVfitMassUncTauUp", (float) SVfitMassUncTauUp);
    pair.addUserFloat("SVfitMassUncTauDown", (float) SVfitMassUncTauDown);

    pair.addUserFloat("SVfitTransverseMassUnc", (float) SVfitTransverseMassUnc);
    pair.addUserFloat("SVfitTransverseMassUncTauUp", (float) SVfitTransverseMassUncTauUp);
    pair.addUserFloat("SVfitTransverseMassUncTauDown", (float) SVfitTransverseMassUncTauDown);

    pair.addUserFloat("SVfit_ptUnc", (float) SVptUnc);
    pair.addUserFloat("SVfit_ptUncTauUp", (float) SVptUncTauUp);
    pair.addUserFloat("SVfit_ptUncTauDown", (float) SVptUncTauDown);

    pair.addUserFloat("SVfit_etaUnc", (float) SVetaUnc);
    pair.addUserFloat("SVfit_etaUncTauUp", (float) SVetaUncTauUp);
    pair.addUserFloat("SVfit_etaUncTauDown", (float) SVetaUncTauDown);

    pair.addUserFloat("SVfit_phiUnc", (float) SVphiUnc);
    pair.addUserFloat("SVfit_phiUncTauUp", (float) SVphiUncTauUp);
    pair.addUserFloat("SVfit_phiUncTauDown", (float) SVphiUncTauDown);

    pair.addUserFloat("SVfit_METRho", (float) SVMETRho);
    pair.addUserFloat("SVfit_METRhoTauUp", (float) SVMETRhoTauUp);
    pair.addUserFloat("SVfit_METRhoTauDown", (float) SVMETRhoTauDown);

    pair.addUserFloat("SVfit_METPhi", (float) SVMETPhi);
    pair.addUserFloat("SVfit_METPhiTauUp", (float) SVMETPhiTauUp);
    pair.addUserFloat("SVfit_METPhiTauDown", (float) SVMETPhiTauDown);

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

// decide if leptons 1 and 2 must be switched to respect SVfit conventions
bool SVfitInterface::Switch (svFitStandalone::kDecayType type1, double pt1, svFitStandalone::kDecayType type2, double pt2)
{
    // e e, mu mu, tau tau
    if (type1 == type2) {return (pt1 < pt2);}
    
    // e tau, mu tau
    if ( (type1 == svFitStandalone::kTauToElecDecay || type1 == svFitStandalone::kTauToMuDecay) &&
         type2 == svFitStandalone::kTauToHadDecay ) {return false;}
    if ( (type2 == svFitStandalone::kTauToElecDecay || type2 == svFitStandalone::kTauToMuDecay) &&
         type1 == svFitStandalone::kTauToHadDecay ) {return true;}

    // e mu
    if (type1 == svFitStandalone::kTauToElecDecay && type2 == svFitStandalone::kTauToMuDecay) {return false;}
    if (type2 == svFitStandalone::kTauToElecDecay && type1 == svFitStandalone::kTauToMuDecay) {return true;}
    
    cout << "SVfit Standalone: ordering not done (should never happen)" << endl;
    return false;
}

// set mass (pdg ele/mu or leave cand mass for tauh
double SVfitInterface::GetMass (svFitStandalone::kDecayType type, double candMass)
{
    if (type == svFitStandalone::kTauToElecDecay) return 0.51100e-3;
    if (type == svFitStandalone::kTauToMuDecay) return 0.10566;

    return candMass; // for tauh and all exceptions return cand mass
}

bool SVfitInterface::IsInteresting (const reco::Candidate *l1, const reco::Candidate *l2)
{
  int apdg1 = abs(l1->pdgId());
  int apdg2 = abs(l2->pdgId());

  int nmu = 0;
  int nele = 0;
  int ntau = 0;

  if (apdg1 == 13) nmu++;
  if (apdg1 == 11) nele++;
  if (apdg1 == 15) ntau++;

  if (apdg2 == 13) nmu++;
  if (apdg2 == 11) nele++;
  if (apdg2 == 15) ntau++;

  pairType pType = kOther;
  if (nmu == 1 && nele == 0 && ntau == 1) pType = kMuHad;
  if (nmu == 0 && nele == 1 && ntau == 1) pType = kEHad;
  if (nmu == 0 && nele == 0 && ntau == 2) pType = kHadHad;
  if (nmu == 2 && nele == 0 && ntau == 0) pType = kMuMu;
  if (nmu == 0 && nele == 2 && ntau == 0) pType = kEE;
  if (nmu == 1 && nele == 1 && ntau == 0) pType = kEMu;

  ///////

  // switch to apply different requirements to the objects
  if (deltaR(l1->p4(), l2->p4()) < 0.1)
    return false; // for overlap removal

  ///////

  // create pointers with usual pair ordering -- easier to apply cuts later
  const reco::Candidate* dau1;
  const reco::Candidate* dau2;

  if (pType == kMuHad)
  {
    dau1 = (apdg1 == 13 ? l1 : l2);
    dau2 = (apdg1 == 13 ? l2 : l1);

    if (dau1->pt() < 17.)
      return false;

    if (dau2->pt() < 20.)
      return false;

    if (userdatahelpers::getUserInt(l2,"decayModeFinding") != 1) // decayModeFinding == decayModeFindingOldDMs
      return false;

    bool iso1 = (userdatahelpers::getUserFloat(l1,"combRelIsoPF") < 0.3);
    //bool iso2 = (userdatahelpers::getUserInt(l2,"byVLooseIsolationMVArun2v1DBoldDMwLT") == 1);
    bool iso2 = (userdatahelpers::getUserInt(l2,"byVLooseIsolationMVArun2017v2DBoldDMwLT2017") == 1); //FRA 2017

    if (!iso1 || !iso2) 
      return false;

    return true; // passed all requirements
  }

  else if (pType == kEHad)
  {
    dau1 = (apdg1 == 11 ? l1 : l2);
    dau2 = (apdg1 == 11 ? l2 : l1);

    if (dau1->pt() < 19.)
      return false;

    if (dau2->pt() < 20.)
      return false;

    if (userdatahelpers::getUserInt(l2,"decayModeFinding") != 1)  // decayModeFinding == decayModeFindingOldDMs
      return false;

    bool iso1 = (userdatahelpers::getUserFloat(l1,"combRelIsoPF") < 0.3);
    //bool iso2 = (userdatahelpers::getUserInt(l2,"byVLooseIsolationMVArun2v1DBoldDMwLT") == 1);
    bool iso2 = (userdatahelpers::getUserInt(l2,"byVLooseIsolationMVArun2017v2DBoldDMwLT2017") == 1); //FRA 2017

    if (!iso1 || !iso2) 
      return false;

    return true; // passed all requirements
  }

  else if (pType == kHadHad)
  {
    dau1 = ((l1->pt() > l2->pt()) ? l1 : l2);
    dau2 = ((l1->pt() > l2->pt()) ? l2 : l1);

    if (dau1->pt() < 30.)
      return false;
    
    if (dau2->pt() < 30.)
      return false;
    
    if (userdatahelpers::getUserInt(l1,"decayModeFinding") != 1)  // decayModeFinding == decayModeFindingOldDMs
      return false;
    
    if (userdatahelpers::getUserInt(l2,"decayModeFinding") != 1)  // decayModeFinding == decayModeFindingOldDMs
      return false;

    //bool iso1 = (userdatahelpers::getUserInt(l1,"byVLooseIsolationMVArun2v1DBoldDMwLT") == 1);
    //bool iso2 = (userdatahelpers::getUserInt(l2,"byVLooseIsolationMVArun2v1DBoldDMwLT") == 1);
    bool iso1 = (userdatahelpers::getUserInt(l1,"byVLooseIsolationMVArun2017v2DBoldDMwLT2017") == 1); // FRA 2017
    bool iso2 = (userdatahelpers::getUserInt(l2,"byVLooseIsolationMVArun2017v2DBoldDMwLT2017") == 1); // FRA 2017

    if (!iso1 || !iso2) 
      return false;

    return true; // passed all requirements
  }

  else if (pType == kMuMu)
    return false;
  
  else if (pType == kEE)
    return false;
  
  else if (pType == kEMu)
    return false;
  
  else
  {
    // should never happen
    edm::LogWarning("Unrecognised pair") << "(SVfitInterface) Warning! could not assess the pair type, won't compute SVFit";
    return false;
  }
}


#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(SVfitInterface);
