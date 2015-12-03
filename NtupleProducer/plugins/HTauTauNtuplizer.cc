/** \class HTauTauNtuplizer
 *
 *  No description available.
 *
 *  $Date: 2014/10/31 10:08:19 $
 *  $Revision: 1.00 $
 *  \author G. Ortona (LLR) and L. Cadamuro (LLR)
 */

#define EPSIL 1.e-5 

// system include files
#include <memory>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <map>
#include <utility>
#include <TNtuple.h>
//#include <XYZTLorentzVector.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/Framework/interface/LuminosityBlock.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Common/interface/TriggerNames.h>

#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/MET.h>
#include "DataFormats/PatCandidates/interface/Jet.h"
#include <DataFormats/PatCandidates/interface/GenericParticle.h>
#include "DataFormats/VertexReco/interface/Vertex.h"
#include <DataFormats/Common/interface/View.h>
#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/METReco/interface/PFMET.h>
#include <DataFormats/METReco/interface/PFMETCollection.h>
#include <DataFormats/JetReco/interface/PFJet.h>
#include <DataFormats/JetReco/interface/PFJetCollection.h>
#include <DataFormats/Math/interface/LorentzVector.h>
#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/Common/interface/MergeableCounter.h>
#include "DataFormats/Math/interface/deltaR.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenFilterInfo.h"
#include <CommonTools/UtilAlgos/interface/TFileService.h>

#include <Muon/MuonAnalysisTools/interface/MuonEffectiveArea.h>

#include <LLRHiggsTauTau/NtupleProducer/interface/CutSet.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/LeptonIsoHelper.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/DaughterDataHelpers.h>
//#include <LLRHiggsTauTau/NtupleProducer/interface/FinalStates.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/MCHistoryTools.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/PUReweight.h>
//#include <LLRHiggsTauTau/NtupleProducer/interface/VBFCandidateJetSelector.h>
//#include <LLRHiggsTauTau/NtupleProducer/interface/bitops.h>
//#include <LLRHiggsTauTau/NtupleProducer/interface/Fisher.h>
//#include <LLRHiggsTauTau/NtupleProducer/interface/HTauTauConfigHelper.h>
//#include "HZZ4lNtupleFactory.h"
#include <LLRHiggsTauTau/NtupleProducer/interface/PhotonFwd.h>
#include "LLRHiggsTauTau/NtupleProducer/interface/triggerhelper.h"
#include "LLRHiggsTauTau/NtupleProducer/interface/PUReweight.h"
//#include "LLRHiggsTauTau/NtupleProducer/Utils/OfflineProducerHelper.h"
#include "LLRHiggsTauTau/NtupleProducer/interface/ParticleType.h"
#include "LLRHiggsTauTau/NtupleProducer/interface/GenFlags.h"
#include "LLRHiggsTauTau/NtupleProducer/interface/GenHelper.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"


#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"

#include "TLorentzVector.h"

 namespace {
//   bool writePhotons = false;  // Write photons in the tree. 
   bool writeJets = true;     // Write jets in the tree. 
   bool writeSoftLep = false;
   bool DEBUG = false;
 }

using namespace std;
using namespace edm;
using namespace reco;


// class declaration

class HTauTauNtuplizer : public edm::EDAnalyzer {
 public:
  /// Constructor
  explicit HTauTauNtuplizer(const edm::ParameterSet&);
    
  /// Destructor
  virtual ~HTauTauNtuplizer();  

 private:
  //----edm control---
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  //void InitializeBranches();
  //void InitializeVariables();
  void Initialize(); 
  Int_t FindCandIndex(const reco::Candidate&, Int_t iCand);
  //----To implement here-----
  //virtual void FillCandidate(const pat::CompositeCandidate& higgs, bool evtPass, const edm::Event&, const Int_t CRflag);
  //virtual void FillPhoton(const pat::Photon& photon);
  int FillJet(const edm::View<pat::Jet>* jet, const edm::Event&, JetCorrectionUncertainty*);
  void FillSoftLeptons(const edm::View<reco::Candidate> *dauhandler, const edm::Event& event, bool theFSR);
  //void FillbQuarks(const edm::Event&);
  void FillGenInfo(const edm::Event&);
  void FillGenJetInfo(const edm::Event&);
  int GetMatchedGen (const reco::Candidate* genL, const edm::Event& event); // return the index of the associated gen particle in the filtered gen collection, in not existing return -1
  //int CreateFlagsWord (const pat::GenericParticle* part); // build int with each bit containing some boolean flags
  static bool CompareLegs(const reco::Candidate *, const reco::Candidate *);
  float ComputeMT (math::XYZTLorentzVector visP4, float METx, float METy);
  static bool ComparePairsbyPt(pat::CompositeCandidate i, pat::CompositeCandidate j);
  static bool ComparePairsbyIso(pat::CompositeCandidate i, pat::CompositeCandidate j);

  // ----------member data ---------------------------
  //std::map <int, int> genFlagPosMap_; // to convert from input to output enum format for H/Z decays
  
  //Configs
  int theChannel;
  std::string theCandLabel;
  TString theFileName;
  bool theFSR;
  Bool_t theisMC;
  Bool_t theUseNoHFPFMet; // false: PFmet ; true: NoHFPFMet
  //Trigger
  vector<int> indexOfPath;
  vector<string> foundPaths;
  //Int_t nFoundPaths;
  //edm::InputTag triggerResultsLabel;
  edm::InputTag processName;
  //edm::InputTag triggerSet;

  HLTConfigProvider hltConfig_;
  //Output Objects
  TTree *myTree;//->See from ntuplefactory in zz4l
  TH1F *hCounter;
  TH1F *hTauIDs;
  triggerhelper* myTriggerHelper;

  PUReweight reweight;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;

  //flags
  //static const int nOutVars =14;
  bool applyTrigger;    // Only events passing trigger
  bool applySkim;       //  "     "      "     skim
  bool skipEmptyEvents; // Skip events whith no candidate in the collection
  //PUReweight reweight;

  //counters
  Int_t Nevt_Gen;
  Int_t Nevt_PassTrigger;
  Int_t Npairs;

  //Event Output variables
  Int_t _indexevents;
  Int_t _runNumber;
  Int_t _lumi;
  Long64_t _triggerbit;
  Int_t _metfilterbit;
  Float_t _met;
  Float_t _metphi;
  Float_t _MC_weight;
  Float_t _aMCatNLOweight;
  Int_t _npv;
  Float_t _lheHt;
  Int_t _npu;
  Float_t _PUNumInteractions;
  Float_t _PUReweight;
  Float_t _rho;
  Int_t _nup;
  
  // pairs
  //std::vector<TLorentzVector> _mothers;
  std::vector<Float_t> _mothers_px;
  std::vector<Float_t> _mothers_py;
  std::vector<Float_t> _mothers_pz;
  std::vector<Float_t> _mothers_e;
  
  // reco leptons
  //std::vector<TLorentzVector> _daughters;
  std::vector<string> _trigger_name;
  std::vector<Int_t> _trigger_accept;
  std::vector<Float_t> _daughters_px;
  std::vector<Float_t> _daughters_py;
  std::vector<Float_t> _daughters_pz;
  std::vector<Float_t> _daughters_e;
  std::vector<Int_t> _daughters_genindex;
  std::vector<Int_t> _daughters_charge;

  std::vector<const reco::Candidate*> _softLeptons;
  
  //std::vector<TLorentzVector> _bquarks;
  //std::vector<Float_t> _bquarks_px;
  //std::vector<Float_t> _bquarks_py;
  //std::vector<Float_t> _bquarks_pz;
  //std::vector<Float_t> _bquarks_e;
  //std::vector<Int_t> _bquarks_pdg;
  
  std::vector<Float_t> _genpart_px;
  std::vector<Float_t> _genpart_py;
  std::vector<Float_t> _genpart_pz;
  std::vector<Float_t> _genpart_e;
  std::vector<Int_t> _genpart_pdg;
  std::vector<Int_t> _genpart_status;
  //std::vector<Int_t> _genpart_mothInd;
  std::vector<Int_t> _genpart_HMothInd;
  std::vector<Int_t> _genpart_MSSMHMothInd;
  std::vector<Int_t> _genpart_TopMothInd;
  std::vector<Int_t> _genpart_TauMothInd;
  std::vector<Int_t> _genpart_ZMothInd;
  std::vector<Int_t> _genpart_WMothInd;
  std::vector<Int_t> _genpart_bMothInd;
  std::vector<Int_t> _genpart_HZDecayMode;
  std::vector<Int_t> _genpart_WDecayMode;
  std::vector<Int_t> _genpart_TopDecayMode;
  std::vector<Int_t> _genpart_TauGenDecayMode;

  std::vector<Int_t> _genpart_flags; // vector of bit flags bout gen info
  
  // gen jets
  std::vector<Float_t> _genjet_px;
  std::vector<Float_t> _genjet_py;
  std::vector<Float_t> _genjet_pz;
  std::vector<Float_t> _genjet_e;
  std::vector<Int_t> _genjet_partonFlavour; // from matched pat::Jet
  std::vector<Int_t> _genjet_hadronFlavour; // (eh yes, because it is not accessible easily from gen jets)

  //std::vector<math::XYZTLorentzVector> _daughter2;

  //Mothers output variables
  std::vector<Int_t> _indexDau1;
  std::vector<Int_t> _indexDau2;
  std::vector<Float_t> _daughters_HLTpt;
  //std::vector<Int_t> _genDaughters;
  std::vector<Bool_t> _isOSCand;
  std::vector<Float_t> _SVmass;
  std::vector<Float_t> _SVpt;
  std::vector<Float_t> _SVptUnc;
  std::vector<Float_t> _SVeta;
  std::vector<Float_t> _SVetaUnc;
  std::vector<Float_t> _SVphi;
  std::vector<Float_t> _SVphiUnc;
  std::vector<Float_t> _SVMetRho;
  std::vector<Float_t> _SVMetPhi;
  std::vector<Float_t> _metx;
  std::vector<Float_t> _mety;
  std::vector<Float_t> _metCov00;
  std::vector<Float_t> _metCov01;
  std::vector<Float_t> _metCov10;
  std::vector<Float_t> _metCov11;
  std::vector<Float_t> _metSignif;
  std::vector<Float_t> _mTDau1;
  std::vector<Float_t> _mTDau2;
  //std::vector<Float_t> _bmotmass;
   
  //Leptons variables
  std::vector<Int_t> _pdgdau;
  std::vector<Int_t> _particleType;//0=muon, 1=e, 2=tau
  std::vector<Float_t> _combreliso;
  std::vector<Float_t> _discriminator;//BDT for ele, discriminator for tau, 
  std::vector<Int_t> _daughters_muonID; //bitwise (bit 0 loose, 1 soft , 2 medium, 3 tight, 4 highPT 5 tight_noVtx)
  std::vector<Int_t> _daughters_typeOfMuon; //bitwise, 0=PF, 1=Global, 2=Tracker
  std::vector<Float_t> _dxy;
  std::vector<Float_t> _dz;
  std::vector<Float_t> _SIP;
  std::vector<bool> _daughters_iseleBDT; //isBDT for ele
  std::vector<bool> _daughters_iseleWP80; //isBDT for ele
  std::vector<bool> _daughters_iseleWP90; //isBDT for ele
  std::vector<Float_t> _daughters_eleMVAnt; //isBDT for ele
  std::vector<bool> _daughters_passConversionVeto; //isBDT for ele
  std::vector<int>  _daughters_eleMissingHits;
  std::vector<int> _daughters_iseleCUT; //CUT ID for ele (0=veto,1=loose,2=medium,3=tight)
  std::vector<Int_t> _decayType;//for taus only
  std::vector<Long64_t> _daughters_tauID; //bitwise. check h_tauID for histogram list 
  static const int ntauIds = 26;
  TString tauIDStrings[ntauIds] = {
   "byLooseIsolationMVA3newDMwLT",
   "byLooseIsolationMVA3oldDMwLT",
   "byLoosePileupWeightedIsolation3Hits",
   "byMediumIsolationMVA3newDMwLT",
   "byMediumIsolationMVA3oldDMwLT",
   "byMediumPileupWeightedIsolation3Hits",
   "byPhotonPtSumOutsideSignalCone",
   "byTightIsolationMVA3newDMwLT",
   "byTightIsolationMVA3oldDMwLT",
   "byTightPileupWeightedIsolation3Hits",
   "byVLooseIsolationMVA3newDMwLT",
   "byVLooseIsolationMVA3oldDMwLT",
   "byVTightIsolationMVA3newDMwLT",
   "byVTightIsolationMVA3oldDMwLT",
   "byVVTightIsolationMVA3newDMwLT",
   "byVVTightIsolationMVA3oldDMwLT",
   "byLooseCombinedIsolationDeltaBetaCorr3Hits",
   "byMediumCombinedIsolationDeltaBetaCorr3Hits",
   "byTightCombinedIsolationDeltaBetaCorr3Hits",
   "againstMuonLoose3",
   "againstMuonTight3",
   "againstElectronVLooseMVA5",
   "againstElectronLooseMVA5",
   "againstElectronMediumMVA5",
   "againstElectronTightMVA5",
   "againstElectronVTightMVA5"
  };
  std::vector<Float_t> _daughters_IetaIeta;
  std::vector<Float_t> _daughters_deltaPhiSuperClusterTrackAtVtx;
  std::vector<Float_t> _daughters_SCeta;
  std::vector<Float_t> _daughters_depositR03_tracker;
  std::vector<Float_t> _daughters_depositR03_ecal;
  std::vector<Float_t> _daughters_depositR03_hcal;
  std::vector<Int_t> _daughters_decayModeFindingOldDMs;

  std::vector<Float_t> _daughters_againstElectronMVA5category;
  std::vector<Float_t> _daughters_againstElectronMVA5raw;
  std::vector<Float_t> _daughters_byPileupWeightedIsolationRaw3Hits;
  std::vector<Float_t> _daughters_footprintCorrection;
  std::vector<Float_t> _daughters_neutralIsoPtSumWeight;
  std::vector<Float_t> _daughters_photonPtSumOutsideSignalCone;

  std::vector<Int_t> _daughters_decayModeFindingNewDMs;
  std::vector<Float_t> _daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits;
  std::vector<Float_t> _daughters_byIsolationMVA3oldDMwoLTraw;
  std::vector<Float_t> _daughters_byIsolationMVA3oldDMwLTraw;
  std::vector<Float_t> _daughters_byIsolationMVA3newDMwoLTraw;
  std::vector<Float_t> _daughters_byIsolationMVA3newDMwLTraw;
  std::vector<Float_t> _daughters_chargedIsoPtSum;
  std::vector<Float_t> _daughters_neutralIsoPtSum;
  std::vector<Float_t> _daughters_puCorrPtSum;
  std::vector<Int_t> _daughters_numChargedParticlesSignalCone;
  std::vector<Int_t> _daughters_numNeutralHadronsSignalCone;
  std::vector<Int_t> _daughters_numPhotonsSignalCone;
  std::vector<Int_t> _daughters_numParticlesSignalCone;
  std::vector<Int_t> _daughters_numChargedParticlesIsoCone;
  std::vector<Int_t> _daughters_numNeutralHadronsIsoCone;
  std::vector<Int_t> _daughters_numPhotonsIsoCone;
  std::vector<Int_t> _daughters_numParticlesIsoCone;
  std::vector<Float_t> _daughters_leadChargedParticlePt;
  std::vector<Float_t> _daughters_trackRefPt;
  std::vector<Int_t> _daughters_LFtrigger;
  std::vector<Int_t> _daughters_L3trigger;
  std::vector<Int_t> _daughters_FilterFired;
  std::vector<Int_t> _daughters_isGoodTriggerType;
  std::vector<Int_t> _daughters_L3FilterFired;
  std::vector<Int_t> _daughters_L3FilterFiredLast;

  //Jets variables
  Int_t _numberOfJets;
  //std::vector<TLorentzVector> _jets;
  std::vector<Float_t> _jets_px;
  std::vector<Float_t> _jets_py;
  std::vector<Float_t> _jets_pz;
  std::vector<Float_t> _jets_e;
  std::vector<Float_t> _jets_rawPt;
  std::vector<Float_t> _jets_mT;
  std::vector<Float_t> _jets_PUJetID;
  std::vector<Float_t> _jets_vtxPt;
  std::vector<Float_t> _jets_vtxMass;
  std::vector<Float_t> _jets_vtx3dL;
  std::vector<Float_t> _jets_vtxNtrk;
  std::vector<Float_t> _jets_vtx3deL;
  std::vector<Float_t> _jets_leadTrackPt;
  std::vector<Float_t> _jets_leptonPtRel; 
  std::vector<Float_t> _jets_leptonPt;    
  std::vector<Float_t> _jets_leptonDeltaR;
  std::vector<Float_t> _jets_chEmEF;
  std::vector<Float_t> _jets_chHEF;
  std::vector<Float_t> _jets_nEmEF;
  std::vector<Float_t> _jets_nHEF;
  std::vector<Int_t> _jets_chMult;
  std::vector<Float_t> _jets_jecUnc;

  std::vector<Int_t> _jets_Flavour; // parton flavour
  std::vector<Int_t> _jets_HadronFlavour; // hadron flavour
  std::vector<Int_t> _jets_genjetIndex; // index of matched gen jet in genjet vector
  std::vector<Float_t> _bdiscr;
  std::vector<Float_t> _bdiscr2;
  std::vector<Int_t> _jetID; //1=loose, 2=tight, 3=tightlepveto
  std::vector<Float_t> _jetrawf;

  //genH
  //std::vector<Float_t> _genH_px;
  //std::vector<Float_t> _genH_py;
  //std::vector<Float_t> _genH_pz;
  //std::vector<Float_t> _genH_e;
};

// ----Constructor and Destructor -----
HTauTauNtuplizer::HTauTauNtuplizer(const edm::ParameterSet& pset) : reweight(),
  triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(pset.getParameter<edm::InputTag>("triggerSet"))),
  triggerBits_(consumes<edm::TriggerResults>(pset.getParameter<edm::InputTag>("triggerResultsLabel")))

 {
  theCandLabel = pset.getUntrackedParameter<string>("CandCollection");
  //theChannel = myHelper.channel();
  theFileName = pset.getUntrackedParameter<string>("fileName");
  skipEmptyEvents = pset.getParameter<bool>("skipEmptyEvents");
  theFSR = pset.getParameter<bool>("applyFSR");
  theisMC = pset.getParameter<bool>("IsMC");
  theUseNoHFPFMet = pset.getParameter<bool>("useNOHFMet");
  //writeBestCandOnly = pset.getParameter<bool>("onlyBestCandidate");
  //sampleName = pset.getParameter<string>("sampleName");
  Nevt_Gen=0;
  Nevt_PassTrigger = 0;
  Npairs=0;
  
  ///// TRIGGER
  //triggerResultsLabel = InputTag("TriggerResults","","HLT");
  processName= pset.getParameter<edm::InputTag>("triggerResultsLabel");
  std::vector<edm::ParameterSet> HLTList = pset.getParameter <std::vector<edm::ParameterSet> > ("triggerList");
  //std::vector<std::string> 

  myTriggerHelper = new triggerhelper();// (HLTList);
  for (std::vector<edm::ParameterSet>::const_iterator iPSet = HLTList.begin();iPSet != HLTList.end(); ++iPSet) {
    const std::string& hlt = iPSet->getParameter<std::string>("HLT");
    const std::vector<std::string>& path1 = iPSet->getParameter<std::vector<std::string>>("path1");
    const std::vector<std::string>& path2 = iPSet->getParameter<std::vector<std::string>>("path2");
    const int& chan = iPSet->getParameter<int>("channel");
    // Build the mape
    myTriggerHelper->addTriggerMap(hlt,path1,path2,chan);
  }

  //triggerSet= pset.getParameter<edm::InputTag>("triggerSet");




  /*
  // init map for flags
  genFlagPosMap_.insert (std::pair<int,int> ( static_cast<int> (genhelper::HZDecay::MuHad)  , static_cast<int> (GenFlags::HZToMuHad) ) );
  genFlagPosMap_.insert (std::pair<int,int> ( static_cast<int> (genhelper::HZDecay::EHad)   , static_cast<int> (GenFlags::HZToEHad) ) );
  genFlagPosMap_.insert (std::pair<int,int> ( static_cast<int> (genhelper::HZDecay::HadHad) , static_cast<int> (GenFlags::HZToHadHad) ) );
  genFlagPosMap_.insert (std::pair<int,int> ( static_cast<int> (genhelper::HZDecay::MuMu)   , static_cast<int> (GenFlags::HZToMuMu) ) );
  genFlagPosMap_.insert (std::pair<int,int> ( static_cast<int> (genhelper::HZDecay::EE)     , static_cast<int> (GenFlags::HZToEE) ) );
  genFlagPosMap_.insert (std::pair<int,int> ( static_cast<int> (genhelper::HZDecay::EMu)    , static_cast<int> (GenFlags::HZToEMu) ) );
  genFlagPosMap_.insert (std::pair<int,int> ( static_cast<int> (genhelper::HZDecay::EEPrompt)   , static_cast<int> (GenFlags::HZToEEPrompt) ) );
  genFlagPosMap_.insert (std::pair<int,int> ( static_cast<int> (genhelper::HZDecay::MuMuPrompt) , static_cast<int> (GenFlags::HZToMuMuPrompt) ) );
  genFlagPosMap_.insert (std::pair<int,int> ( static_cast<int> (genhelper::HZDecay::Other)      , static_cast<int> (GenFlags::HZToOther) ) );
  */

  Initialize();
}

HTauTauNtuplizer::~HTauTauNtuplizer(){
  delete myTriggerHelper;
}
//

void HTauTauNtuplizer::Initialize(){
  //_mothers.clear();
  _mothers_px.clear();
  _mothers_py.clear();
  _mothers_pz.clear();
  _mothers_e.clear();
  
  //_daughters.clear();
  if(DEBUG){
  _trigger_name.clear();
  _trigger_accept.clear();
  }
  _daughters_px.clear();
  _daughters_py.clear();
  _daughters_pz.clear();
  _daughters_e.clear();
  _daughters_charge.clear();
  _daughters_genindex.clear();
  _daughters_IetaIeta.clear();
  _daughters_deltaPhiSuperClusterTrackAtVtx.clear();
  _daughters_SCeta.clear();
  _daughters_depositR03_tracker.clear();  
  _daughters_depositR03_ecal.clear();  
  _daughters_depositR03_hcal.clear();  
  _daughters_decayModeFindingOldDMs.clear();
  _daughters_decayModeFindingNewDMs.clear();
  _daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits.clear();
  _daughters_byIsolationMVA3oldDMwoLTraw.clear();
  _daughters_byIsolationMVA3oldDMwLTraw.clear();
  _daughters_byIsolationMVA3newDMwoLTraw.clear();

  _daughters_againstElectronMVA5category.clear();
  _daughters_againstElectronMVA5raw.clear();
  _daughters_byPileupWeightedIsolationRaw3Hits.clear();
  _daughters_footprintCorrection.clear();
  _daughters_neutralIsoPtSumWeight.clear();
  _daughters_photonPtSumOutsideSignalCone.clear();

  _daughters_byIsolationMVA3newDMwLTraw.clear();
  _daughters_chargedIsoPtSum.clear();
  _daughters_neutralIsoPtSum.clear();
  _daughters_puCorrPtSum.clear();
  _daughters_numChargedParticlesSignalCone.clear();
  _daughters_numNeutralHadronsSignalCone.clear();
  _daughters_numPhotonsSignalCone.clear();
  _daughters_numParticlesSignalCone.clear();
  _daughters_numChargedParticlesIsoCone.clear();
  _daughters_numNeutralHadronsIsoCone.clear();
  _daughters_numPhotonsIsoCone.clear();
  _daughters_numParticlesIsoCone.clear();
  _daughters_leadChargedParticlePt.clear();
  _daughters_trackRefPt.clear();
  _daughters_LFtrigger.clear();
  _daughters_L3trigger.clear();
  _daughters_FilterFired.clear();
  _daughters_isGoodTriggerType.clear();
  _daughters_L3FilterFired.clear();
  _daughters_L3FilterFiredLast.clear();
  _daughters_iseleBDT.clear();
  _daughters_iseleWP80.clear();
  _daughters_iseleWP90.clear();
  _daughters_eleMVAnt.clear();
  _daughters_passConversionVeto.clear();
  _daughters_eleMissingHits.clear();
  _daughters_iseleCUT.clear();
  //_daughter2.clear();
  _softLeptons.clear();
  //_genDaughters.clear();
  
  //_bquarks.clear();
  //_bquarks_px.clear();
  //_bquarks_py.clear();
  //_bquarks_pz.clear();
  //_bquarks_e.clear();
  //_bquarks_pdg.clear();
  //_bmotmass.clear();
  
  _genpart_px.clear();
  _genpart_py.clear();
  _genpart_pz.clear();
  _genpart_e.clear();
  _genpart_pdg.clear();
  _genpart_status.clear();
  //_genpart_mothInd.clear();
  _genpart_HMothInd.clear();
  _genpart_MSSMHMothInd.clear();
  _genpart_TopMothInd.clear();
  _genpart_TauMothInd.clear();
  _genpart_ZMothInd.clear();
  _genpart_WMothInd.clear();
  _genpart_bMothInd.clear();
  _genpart_HZDecayMode.clear();
  _genpart_TopDecayMode.clear();
  _genpart_WDecayMode.clear();
  _genpart_TauGenDecayMode.clear();
  _genpart_flags.clear();
  
  _genjet_px.clear();
  _genjet_py.clear();
  _genjet_pz.clear();
  _genjet_e.clear();
  _genjet_partonFlavour.clear();
  _genjet_hadronFlavour.clear();
  
  _indexDau1.clear();
  _indexDau2.clear();
  _pdgdau.clear();
  _SVmass.clear();
  _SVpt.clear();
  _SVptUnc.clear();
  _SVeta.clear();
  _SVetaUnc.clear();
  _SVphi.clear();
  _SVphiUnc.clear();
  _SVMetRho.clear();
  _SVMetPhi.clear();
  _isOSCand.clear();
  _daughters_HLTpt.clear();
  _metx.clear();
  _mety.clear();
  _metCov00.clear();
  _metCov01.clear();
  _metCov10.clear();
  _metCov11.clear();
  _metSignif.clear();
  _mTDau1.clear();
  _mTDau2.clear();
  _particleType.clear();
  _discriminator.clear();
  _daughters_typeOfMuon.clear();
  _daughters_muonID.clear();
  _dxy.clear();
  _dz.clear();
  _SIP.clear();
  _decayType.clear();
  _daughters_tauID.clear();
  _combreliso.clear();
  _indexevents=0;
  _runNumber=0;
  _lumi=0;
  _triggerbit=0;
  _metfilterbit=0;
  _met=0;
  _metphi=0.;
  _MC_weight=0.;
  _npv=0;
  _lheHt=0;
  _npu=0;
  _PUNumInteractions=0;
  _PUReweight=0.;
  _rho=0;
  _nup=-999;

//  _jets.clear();
  _jets_px.clear();
  _jets_py.clear();
  _jets_pz.clear();
  _jets_e.clear();
  _jets_rawPt.clear();
  _jets_mT.clear();
  _jets_PUJetID.clear();
  _jets_vtxPt.clear();
  _jets_vtxMass.clear();
  _jets_vtx3dL.clear();
  _jets_vtxNtrk.clear();
  _jets_vtx3deL.clear();
  _jets_leadTrackPt.clear();
  _jets_leptonPtRel.clear();
  _jets_leptonPt.clear();
  _jets_leptonDeltaR.clear();
  _jets_chEmEF.clear();
  _jets_chHEF.clear();
  _jets_nEmEF.clear();
  _jets_nHEF.clear();
  _jets_chMult.clear();
  _jets_Flavour.clear();
  _jets_HadronFlavour.clear();
  _jets_genjetIndex.clear();
  _jets_jecUnc.clear();
  _numberOfJets=0;
  _bdiscr.clear();
  _bdiscr2.clear();
  _jetID.clear();
  _jetrawf.clear();

  //_genH_px.clear();
  //_genH_py.clear();
  //_genH_pz.clear();
  //_genH_e.clear();
}

void HTauTauNtuplizer::beginJob(){
  edm::Service<TFileService> fs;
  myTree = fs->make<TTree>("HTauTauTree","HTauTauTree");
  int nbins=3+(myTriggerHelper->GetNTriggers());
  hCounter = fs->make<TH1F>("Counters","Counters",nbins,0,nbins);
  hTauIDs = fs->make<TH1F>("TauIDs","TauIDs",ntauIds,0,ntauIds);

  //Branches
  myTree->Branch("EventNumber",&_indexevents,"EventNumber/I");
  myTree->Branch("RunNumber",&_runNumber,"RunNumber/I");
  myTree->Branch("lumi",&_lumi,"lumi/I");
  myTree->Branch("triggerbit",&_triggerbit,"triggerbit/L");
  myTree->Branch("metfilterbit",&_metfilterbit,"metfilterbit/I");
  myTree->Branch("met",&_met,"met/F");
  myTree->Branch("metphi",&_metphi,"metphi/F");  
  myTree->Branch("npv",&_npv,"npv/I");  
  myTree->Branch("npu",&_npu,"npu/I"); 
  myTree->Branch("PUReweight",&_PUReweight,"PUReweight/F"); 
  myTree->Branch("rho",&_rho,"rho/F");  
  
  myTree->Branch("mothers_px",&_mothers_px);
  myTree->Branch("mothers_py",&_mothers_py);
  myTree->Branch("mothers_pz",&_mothers_pz);
  myTree->Branch("mothers_e",&_mothers_e);
  if(DEBUG){
  myTree->Branch("trigger_name",&_trigger_name);
  myTree->Branch("trigger_accept",&_trigger_accept);
  }
  myTree->Branch("daughters_px",&_daughters_px);
  myTree->Branch("daughters_py",&_daughters_py);
  myTree->Branch("daughters_pz",&_daughters_pz);
  myTree->Branch("daughters_e",&_daughters_e);
  myTree->Branch("daughters_charge",&_daughters_charge);

  if(writeSoftLep)myTree->Branch("softLeptons",&_softLeptons);
  if(theisMC){
    //myTree->Branch("genDaughters",&_genDaughters);
    //myTree->Branch("quarks_px",&_bquarks_px);
    //myTree->Branch("quarks_py",&_bquarks_py);
    //myTree->Branch("quarks_pz",&_bquarks_pz);
    //myTree->Branch("quarks_e",&_bquarks_e);
    //myTree->Branch("quarks_pdg",&_bquarks_pdg);
    //myTree->Branch("motmass",&_bmotmass);
    //myTree->Branch("genH_px",&_genH_px);
    //myTree->Branch("genH_py",&_genH_py);
    //myTree->Branch("genH_pz",&_genH_pz);
    //myTree->Branch("genH_e",&_genH_e);
    myTree->Branch("PUNumInteractions",&_PUNumInteractions,"PUNumInteractions/F");  
    myTree->Branch("daughters_genindex",&_daughters_genindex);
    myTree->Branch("MC_weight",&_MC_weight,"MC_weight/F");
    myTree->Branch("lheHt",&_lheHt,"lheHt/F");  
    myTree->Branch("aMCatNLOweight",&_aMCatNLOweight,"aMCatNLOweight/F");    
    myTree->Branch("genpart_px", &_genpart_px);
    myTree->Branch("genpart_py", &_genpart_py);
    myTree->Branch("genpart_pz", &_genpart_pz);
    myTree->Branch("genpart_e", &_genpart_e);
    myTree->Branch("genpart_pdg", &_genpart_pdg);
    myTree->Branch("genpart_status", &_genpart_status);
    //myTree->Branch("genpart_mothInd", _genpart_mothInd);
    myTree->Branch("genpart_HMothInd", &_genpart_HMothInd);
    myTree->Branch("genpart_MSSMHMothInd", &_genpart_MSSMHMothInd);
    myTree->Branch("genpart_TopMothInd", &_genpart_TopMothInd);
    myTree->Branch("genpart_TauMothInd", &_genpart_TauMothInd);
    myTree->Branch("genpart_ZMothInd", &_genpart_ZMothInd);
    myTree->Branch("genpart_WMothInd", &_genpart_WMothInd);
    myTree->Branch("genpart_bMothInd", &_genpart_bMothInd);
    myTree->Branch("genpart_HZDecayMode", &_genpart_HZDecayMode);
    myTree->Branch("genpart_TopDecayMode", &_genpart_TopDecayMode);
    myTree->Branch("genpart_WDecayMode", &_genpart_WDecayMode);
    myTree->Branch("genpart_TauGenDecayMode", &_genpart_TauGenDecayMode);
    myTree->Branch("genpart_flags", &_genpart_flags);

    myTree->Branch("genjet_px", &_genjet_px);
    myTree->Branch("genjet_py", &_genjet_py);
    myTree->Branch("genjet_pz", &_genjet_pz);
    myTree->Branch("genjet_e" , &_genjet_e);
    myTree->Branch("genjet_partonFlavour" , &_genjet_partonFlavour);
    myTree->Branch("genjet_hadronFlavour" , &_genjet_hadronFlavour);

    myTree->Branch("NUP", &_nup,"NUP/I");
  }
  //myTree->Branch("daughters2",&_daughter2);
  myTree->Branch("SVfitMass",&_SVmass);
  myTree->Branch("SVfit_pt", &_SVpt);
  myTree->Branch("SVfit_ptUnc", &_SVptUnc);
  myTree->Branch("SVfit_eta", &_SVeta);
  myTree->Branch("SVfit_etaUnc", &_SVetaUnc);
  myTree->Branch("SVfit_phi", &_SVphi);
  myTree->Branch("SVfit_phiUnc", &_SVphiUnc);
  myTree->Branch("SVfit_fitMETRho", &_SVMetRho);
  myTree->Branch("SVfit_fitMETPhi", &_SVMetPhi);

  myTree->Branch("isOSCand",&_isOSCand);
  myTree->Branch("METx",&_metx);
  myTree->Branch("METy",&_mety);
  myTree->Branch("MET_cov00",&_metCov00);
  myTree->Branch("MET_cov01",&_metCov01);
  myTree->Branch("MET_cov10",&_metCov10);
  myTree->Branch("MET_cov11",&_metCov11);  
  myTree->Branch("MET_significance",&_metSignif); 
  myTree->Branch("mT_Dau1",&_mTDau1); 
  myTree->Branch("mT_Dau2",&_mTDau2); 
  myTree->Branch("PDGIdDaughters",&_pdgdau);
  myTree->Branch("indexDau1",&_indexDau1);
  myTree->Branch("indexDau2",&_indexDau2);
  myTree->Branch("particleType",&_particleType);
  myTree->Branch("discriminator",&_discriminator);
  myTree->Branch("daughters_muonID",&_daughters_muonID);
  myTree->Branch("daughters_typeOfMuon",&_daughters_typeOfMuon);
  myTree->Branch("dxy",&_dxy);
  myTree->Branch("dz",&_dz);
  myTree->Branch("SIP",&_SIP);
  myTree->Branch("daughters_iseleBDT",&_daughters_iseleBDT);
  myTree->Branch("daughters_iseleWP80",&_daughters_iseleWP80);
  myTree->Branch("daughters_iseleWP90",&_daughters_iseleWP90);
  myTree->Branch("daughters_eleMVAnt",&_daughters_eleMVAnt);
  myTree->Branch("daughters_passConversionVeto",&_daughters_passConversionVeto);
  myTree->Branch("daughters_eleMissingHits",&_daughters_eleMissingHits);
  myTree->Branch("daughters_eleCUTID",&_daughters_iseleCUT);
  myTree->Branch("decayMode",&_decayType);
  myTree->Branch("tauID",&_daughters_tauID);
  myTree->Branch("combreliso",& _combreliso);
  myTree->Branch("daughters_IetaIeta",&_daughters_IetaIeta);
  myTree->Branch("daughters_deltaPhiSuperClusterTrackAtVtx",&_daughters_deltaPhiSuperClusterTrackAtVtx);
  myTree->Branch("daughters_SCeta",&_daughters_SCeta);
  myTree->Branch("daughters_depositR03_tracker",&_daughters_depositR03_tracker);
  myTree->Branch("daughters_depositR03_ecal",&_daughters_depositR03_ecal);
  myTree->Branch("daughters_depositR03_hcal",&_daughters_depositR03_hcal);
  myTree->Branch("daughters_decayModeFindingOldDMs", &_daughters_decayModeFindingOldDMs);

  myTree->Branch("againstElectronMVA5category",&_daughters_againstElectronMVA5category);
  myTree->Branch("againstElectronMVA5raw",&_daughters_againstElectronMVA5raw);
  myTree->Branch("byPileupWeightedIsolationRaw3Hits",&_daughters_byPileupWeightedIsolationRaw3Hits);
  myTree->Branch("footprintCorrection",&_daughters_footprintCorrection);
  myTree->Branch("neutralIsoPtSumWeight",&_daughters_neutralIsoPtSumWeight);
  myTree->Branch("photonPtSumOutsideSignalCone",&_daughters_photonPtSumOutsideSignalCone);

  myTree->Branch("daughters_decayModeFindingNewDMs", &_daughters_decayModeFindingNewDMs);
  myTree->Branch("daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits", &_daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits);
  myTree->Branch("daughters_byIsolationMVA3oldDMwoLTraw",&_daughters_byIsolationMVA3oldDMwoLTraw);
  myTree->Branch("daughters_byIsolationMVA3oldDMwLTraw",&_daughters_byIsolationMVA3oldDMwLTraw);
  myTree->Branch("daughters_byIsolationMVA3newDMwoLTraw",&_daughters_byIsolationMVA3newDMwoLTraw);
  myTree->Branch("daughters_byIsolationMVA3newDMwLTraw",&_daughters_byIsolationMVA3newDMwLTraw);
  myTree->Branch("daughters_chargedIsoPtSum", &_daughters_chargedIsoPtSum);
  myTree->Branch("daughters_neutralIsoPtSum", &_daughters_neutralIsoPtSum);
  myTree->Branch("daughters_puCorrPtSum", &_daughters_puCorrPtSum);
  myTree->Branch("daughters_numChargedParticlesSignalCone", &_daughters_numChargedParticlesSignalCone);
  myTree->Branch("daughters_numNeutralHadronsSignalCone", &_daughters_numNeutralHadronsSignalCone);
  myTree->Branch("daughters_numPhotonsSignalCone", &_daughters_numPhotonsSignalCone);
  myTree->Branch("daughters_daughters_numParticlesSignalCone", &_daughters_numParticlesSignalCone);
  myTree->Branch("daughters_numChargedParticlesIsoCone", &_daughters_numChargedParticlesIsoCone);
  myTree->Branch("daughters_numNeutralHadronsIsoCone", &_daughters_numNeutralHadronsIsoCone);
  myTree->Branch("daughters_numPhotonsIsoCone", &_daughters_numPhotonsIsoCone);
  myTree->Branch("daughters_numParticlesIsoCone", &_daughters_numParticlesIsoCone);
  myTree->Branch("daughters_leadChargedParticlePt", &_daughters_leadChargedParticlePt);
  myTree->Branch("daughters_trackRefPt", &_daughters_trackRefPt);
  myTree->Branch("daughters_isLastTriggerObjectforPath", &_daughters_LFtrigger);
  myTree->Branch("daughters_isTriggerObjectforPath", &_daughters_L3trigger);
  myTree->Branch("daughters_FilterFired",&_daughters_FilterFired);
  myTree->Branch("daughters_isGoodTriggerType",&_daughters_isGoodTriggerType);
  myTree->Branch("daughters_L3FilterFired",&_daughters_L3FilterFired);
  myTree->Branch("daughters_L3FilterFiredLast",&_daughters_L3FilterFiredLast);
  myTree->Branch("daughters_HLTpt",&_daughters_HLTpt);
  
  myTree->Branch("JetsNumber",&_numberOfJets,"JetsNumber/I");
  myTree->Branch("jets_px",&_jets_px);
  myTree->Branch("jets_py",&_jets_py);
  myTree->Branch("jets_pz",&_jets_pz);
  myTree->Branch("jets_e",&_jets_e);
  myTree->Branch("jets_rawPt", &_jets_rawPt);
  myTree->Branch("jets_mT", &_jets_mT);
  myTree->Branch("jets_Flavour",&_jets_Flavour);
  myTree->Branch("jets_HadronFlavour",&_jets_HadronFlavour);
  myTree->Branch("jets_genjetIndex", &_jets_genjetIndex);
  myTree->Branch("jets_PUJetID",&_jets_PUJetID);
  myTree->Branch("jets_vtxPt", &_jets_vtxPt);
  myTree->Branch("jets_vtxMass", &_jets_vtxMass);
  myTree->Branch("jets_vtx3dL", &_jets_vtx3dL);
  myTree->Branch("jets_vtxNtrk", &_jets_vtxNtrk);
  myTree->Branch("jets_vtx3deL", &_jets_vtx3deL);
  myTree->Branch("jets_leadTrackPt", &_jets_leadTrackPt);
  myTree->Branch("jets_leptonPtRel", &_jets_leptonPtRel);
  myTree->Branch("jets_leptonPt", &_jets_leptonPt);
  myTree->Branch("jets_leptonDeltaR", &_jets_leptonDeltaR);
  myTree->Branch("jets_chEmEF" , &_jets_chEmEF);
  myTree->Branch("jets_chHEF"  , &_jets_chHEF);
  myTree->Branch("jets_nEmEF"  , &_jets_nEmEF);
  myTree->Branch("jets_nHEF"   , &_jets_nHEF);
  myTree->Branch("jets_chMult" , &_jets_chMult);
  myTree->Branch("jets_jecUnc" , &_jets_jecUnc);

  myTree->Branch("bDiscriminator",&_bdiscr);
  myTree->Branch("bCSVscore",&_bdiscr2);
  myTree->Branch("PFjetID",&_jetID);
  myTree->Branch("jetRawf",&_jetrawf);
}

Int_t HTauTauNtuplizer::FindCandIndex(const reco::Candidate& cand,Int_t iCand=0){
  const reco::Candidate *daughter = cand.daughter(iCand);
  for(UInt_t iLeptons=0;iLeptons<_softLeptons.size();iLeptons++){
    //if(daughter==daughterPoint[iLeptons]){
    //if(daughter==_softLeptons.at(iLeptons)){
    if(daughter->masterClone().get()==_softLeptons.at(iLeptons)){
      return iLeptons;
    }
  }
  return -1;
}
// ----Analyzer (main) ----
// ------------ method called for each event  ------------
void HTauTauNtuplizer::analyze(const edm::Event& event, const edm::EventSetup& eSetup)
{
  Initialize();

  Handle<vector<reco::Vertex> >  vertexs;
  event.getByLabel("offlineSlimmedPrimaryVertices",vertexs);
  
  //----------------------------------------------------------------------
  // Analyze MC history. THIS HAS TO BE DONE BEFORE ANY RETURN STATEMENT
  // (eg skim or trigger), in order to update the gen counters correctly!!!
  
  //  std::vector<const reco::Candidate *> genZs;
  // std::vector<const reco::Candidate *> genZLeps;
  
  _npv = vertexs->size();
   if (theisMC) {
    Handle<std::vector< PileupSummaryInfo > >  PupInfo;
    event.getByLabel(edm::InputTag("slimmedAddPileupInfo"), PupInfo);    
    std::vector<PileupSummaryInfo>::const_iterator PVI;
    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      if(PVI->getBunchCrossing() == 0) { 
        _PUNumInteractions  = PVI->getPU_NumInteractions();
        int nTrueInt = PVI->getTrueNumInteractions();
        _npu = nTrueInt;        
        _PUReweight = reweight.weight(2012,2012,nTrueInt);
        break;
      } 
    }
  }


  // pile up information -- rho
  edm::Handle<double> rhoHandle;
  event.getByLabel("fixedGridRhoFastjetAll", rhoHandle);
  _rho = *rhoHandle;


  edm::Handle<edm::TriggerResults> triggerBits;

  event.getByToken(triggerBits_, triggerBits);
  const edm::TriggerNames &names = event.triggerNames(*triggerBits);
  if(DEBUG){
    for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
          _trigger_name.push_back( (string) names.triggerName(i));
          _trigger_accept.push_back( (int) triggerBits->accept(i));
    }
  } 
  if (theisMC) {
     Handle<LHEEventProduct> lheEventProduct;
     
    try {event.getByLabel("externalLHEProducer", lheEventProduct);} catch (...) {;}
    if (lheEventProduct.isValid())
    {
      const lhef::HEPEUP& lheEvent = lheEventProduct->hepeup();
      std::vector<lhef::HEPEUP::FiveVector> lheParticles = lheEvent.PUP;
      double lheHt = 0.;
      size_t numParticles = lheParticles.size();
      for ( size_t idxParticle = 0; idxParticle < numParticles; ++idxParticle ) {
        int absPdgId = TMath::Abs(lheEvent.IDUP[idxParticle]);
        int status = lheEvent.ISTUP[idxParticle];
        if ( status == 1 && ((absPdgId >= 1 &&  absPdgId<= 6) ||  absPdgId== 21) ) { // quarks and gluons
            lheHt += TMath::Sqrt((lheParticles[idxParticle][0])*(lheParticles[idxParticle][0]) + (lheParticles[idxParticle][1])*(lheParticles[idxParticle][1])); // first entry is px, second py
        }
      }
       _lheHt = lheHt;
     //cout<<"lheHt = "<<lheHt<<endl;
    }
    //else cout << "LHE product not found" << endl;
  }
  
  _triggerbit = myTriggerHelper->FindTriggerBit(event,foundPaths,indexOfPath);
  _metfilterbit = myTriggerHelper->FindMETBit(event);
  Long64_t tbit = _triggerbit;
  for(int itr=0;itr<myTriggerHelper->GetNTriggers();itr++) {
    if(myTriggerHelper->IsTriggerFired(tbit,itr)) hCounter->Fill(itr+3);
  }

  //Get candidate collection
  edm::Handle<edm::View<pat::CompositeCandidate>>candHandle;
  edm::Handle<edm::View<reco::Candidate>>dauHandle;
  edm::Handle<edm::View<pat::Jet>>jetHandle;
  edm::Handle<pat::METCollection> metHandle;
  edm::Handle<GenFilterInfo> embeddingWeightHandle;
  edm::Handle<edm::TriggerResults> triggerResults;
 
   // Accessing the JEC uncertainties 
  //ak4  
  edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  eSetup.get<JetCorrectionsRecord>().get("AK4PFchs",JetCorParColl); 
  JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  JetCorrectionUncertainty *jecUnc=0;
  jecUnc = new JetCorrectionUncertainty(JetCorPar);


  // protect in case of events where trigger hasn't fired --> no collection created 
  event.getByLabel(theCandLabel,candHandle);
  if (!candHandle.isValid()) return;
  event.getByLabel(theCandLabel,candHandle);
  event.getByLabel("jets",jetHandle);
  event.getByLabel("softLeptons",dauHandle);
  
  if (theUseNoHFPFMet) event.getByLabel("slimmedMETsNoHF",metHandle);
  else event.getByLabel("slimmedMETs",metHandle);
  
  if(theisMC){
    edm::Handle<LHEEventProduct> lheeventinfo;
    event.getByLabel("LHEEventProduct",lheeventinfo);
    if (lheeventinfo.isValid()) {
      _nup=lheeventinfo->hepeup().NUP;
    }
    edm::Handle<GenEventInfoProduct> genEvt;
    event.getByLabel("generator",genEvt);
    _aMCatNLOweight=genEvt->weight();
    _MC_weight = _aMCatNLOweight; // duplicated
  }

  const edm::View<pat::CompositeCandidate>* cands = candHandle.product();
  const edm::View<reco::Candidate>* daus = dauHandle.product();
  const edm::View<pat::Jet>* jets = jetHandle.product();
  const pat::MET &met = metHandle->front();
  //myNtuple->InitializeVariables();
    
  _indexevents = event.id().event();
  _runNumber = event.id().run();
  _lumi=event.luminosityBlock();
  // _met = met.sumEt(); // scalar sum of the pf candidates
  _met = met.pt();
  _metphi = met.phi();
    
  //Do all the stuff here
  //Compute the variables needed for the output and store them in the ntuple
  if(DEBUG)printf("===New Event===\n");

  //Loop over generated b quarks
  //if(theisMC)FillbQuarks(event);
  if(theisMC)
  {
    FillGenInfo(event); // gen particles
    FillGenJetInfo(event); // gen jets
  }
  //Loop of softleptons and fill them
  FillSoftLeptons(daus,event,theFSR);

  //Loop on Jets
  _numberOfJets = 0;
  if(writeJets)_numberOfJets = FillJet(jets, event,jecUnc);

  //Loop on pairs
  std::vector<pat::CompositeCandidate> candVector;
  for(edm::View<pat::CompositeCandidate>::const_iterator candi = cands->begin(); candi!=cands->end();++candi){
    Npairs++;
    const pat::CompositeCandidate& cand = (*candi); 
    candVector.push_back(cand);
  }
  std::sort(candVector.begin(),candVector.end(),ComparePairsbyIso);
  //std::sort(candVector.begin(),candVector.end(),ComparePairsbyPt);
  for(int iPair=0;iPair<int(candVector.size());iPair++){
    const pat::CompositeCandidate& cand = candVector.at(iPair);
    math::XYZTLorentzVector candp4 = cand.p4();

    float thisMETpx = cand.userFloat("MEt_px");
    float thisMETpy = cand.userFloat("MEt_py");
    
    _SVmass.push_back(cand.userFloat("SVfitMass"));
    _SVpt.push_back(cand.userFloat("SVfit_pt"));
    _SVptUnc.push_back(cand.userFloat("SVfit_ptUnc"));
    _SVeta.push_back(cand.userFloat("SVfit_eta"));
    _SVetaUnc.push_back(cand.userFloat("SVfit_etaUnc"));
    _SVphi.push_back(cand.userFloat("SVfit_phi"));
    _SVphiUnc.push_back(cand.userFloat("SVfit_phiUnc"));
    _SVMetRho.push_back(cand.userFloat("SVfit_METRho"));
    _SVMetPhi.push_back(cand.userFloat("SVfit_METPhi"));

    _metx.push_back(thisMETpx);
    _mety.push_back(thisMETpy);    
    _metCov00.push_back(cand.userFloat("MEt_cov00"));
    _metCov01.push_back(cand.userFloat("MEt_cov01"));
    _metCov10.push_back(cand.userFloat("MEt_cov10"));
    _metCov11.push_back(cand.userFloat("MEt_cov11"));
    _metSignif.push_back(cand.userFloat("MEt_significance"));
    
    //if(DEBUG){
      //motherPoint[iMot]=dynamic_cast<const reco::Candidate*>(&*candi);
      //printf("%p %p %p\n",motherPoint[iMot],cand.daughter(0),cand.daughter(1));
    //}
    int index1=-1,index2=-1;
    float mT1 = -1., mT2 = -1.;
    for(int iCand=0;iCand<2;iCand++){
        //int index=-1;
        const reco::Candidate *daughter = cand.daughter(iCand);
        float mT = ComputeMT ( daughter->p4(), thisMETpx, thisMETpy);
        if(iCand==0)
        {
            index1=FindCandIndex(cand,iCand);
            mT1 = mT;     
        }
        else
        {
            index2=FindCandIndex(cand,iCand);
            mT2 = mT;
        }
        if(theFSR){
            const pat::PFParticle* fsr=0;
            double maxPT=-1;
            const PhotonPtrVector* gammas = userdatahelpers::getUserPhotons(daughter);
            if (gammas!=0) {
                for (PhotonPtrVector::const_iterator g = gammas->begin();g!= gammas->end(); ++g) {
                    //const pat::Photon* gamma = g->get();
                    const pat::PFParticle* gamma = g->get();
                    double pt = gamma->pt();
                    if (pt>maxPT) {
                        maxPT = pt;
                        fsr = gamma;
                    }    
                }
            }
    
            //cand.addUserFloat("dauWithFSR",lepWithFsr); // Index of the cand daughter with associated FSR photon
    
            if (fsr!=0) {
              // Add daughter and set p4.
              candp4+=fsr->p4();
              //pfour+=fsr->p4();
              //      myCand.addDaughter(reco::ShallowCloneCandidate(fsr->masterClone()),"FSR"); //FIXME: fsr does not have a masterClone
              //pat::PFParticle myFsr(fsr);
              //myFsr.setPdgId(22); // Fix: photons that are isFromMu have abs(pdgId)=13!!!
              //cand.addDaughter(myFsr,"FSR");
              /*
              // Recompute iso for leptons with FSR    
              const Candidate* d = cand.daughter(iCand);
              float fsrCorr = 0; // The correction to PFPhotonIso
              if (!fsr->isFromMuon()) { // Type 1 photons should be subtracted from muon iso cones
                double dR = ROOT::Math::VectorUtil::DeltaR(fsr->momentum(),d->momentum());
                if (dR<0.4 && ((d->isMuon() && dR > 0.01) ||
                       (d->isElectron() && (fabs((static_cast<const pat::Electron*>(d->masterClone().get()))->superCluster()->eta()) < 1.479 || dR > 0.08)))) {
                  fsrCorr = fsr->pt();
                }
              }

              float rho = ((d->isMuon())?rhoForMu:rhoForEle);
              float combRelIsoPFCorr =  LeptonIsoHelper::combRelIsoPF(sampleType, setup, rho, d, fsrCorr);
      
              string base;
              stringstream str;
              str << "d" << iCand << ".";
              str >> base;      
              cand.addUserFloat(base+"combRelIsoPFFSRCorr",combRelIsoPFCorr);
              */
            }

            //daughter->setP4(daughter->p4()+fsr->p4());
        }

        //if(iCand==0)_indexDau1.push_back(index);
        //else _indexDau2.push_back(index);	
    }
    if(CompareLegs(cand.daughter(0),cand.daughter(1))){
      _indexDau1.push_back(index1);
      _indexDau2.push_back(index2);
      _mTDau1.push_back (mT1);
      _mTDau2.push_back (mT2);
    }else {
      _indexDau1.push_back(index2);
      _indexDau2.push_back(index1);	    
      _mTDau1.push_back (mT2);
      _mTDau2.push_back (mT1);
    }

    // do not read pair charge: sometimes taus have charge +/- 3 and this spoils the pair charge
    //if(fabs(cand.charge())>0.5)_isOSCand.push_back(false);
    //else _isOSCand.push_back(true);
    int chDau1 = (cand.daughter(0))->charge();
    int chDau2 = (cand.daughter(1))->charge();
    bool isOS = (chDau1/abs(chDau1) != chDau2/abs(chDau2));
    _isOSCand.push_back(isOS);

//    if(cand.charge()!=cand.daughter(0)->charge()+cand.daughter(1)->charge())cout<<"charge DIVERSA!!!!!!!!! "<<cand.charge()<<" "<<cand.daughter(0)->charge()<<" "<<cand.daughter(1)->charge()<<endl;
//    else cout<<"charge uguale "<<endl;
    _mothers_px.push_back( (float) candp4.X());
    _mothers_py.push_back( (float) candp4.Y());
    _mothers_pz.push_back( (float) candp4.Z());
    _mothers_e.push_back( (float) candp4.T());
    
    /*   float fillArray[nOutVars]={
    (float)event.id().run(),
    (float)event.id().event(),
    (float)cand.p4().mass(),
      (float)cand.p4().pt(),
      (float)cand.p4().eta(),
      (float)cand.p4().phi(),
      (float)cand.daughter(0)->mass(),
      (float)cand.daughter(0)->pt(),
      (float)cand.daughter(0)->eta(),
      (float)cand.daughter(0)->phi(),
      (float)cand.daughter(1)->mass(),
      (float)cand.daughter(1)->pt(),
      (float)cand.daughter(1)->eta(),
      (float)cand.daughter(1)->phi()
    };
    myTree->Fill(fillArray);*/
    //iMot++;
  }
  myTree->Fill();
  //return;
}

//Fill jets quantities
int HTauTauNtuplizer::FillJet(const edm::View<pat::Jet> *jets, const edm::Event& event, JetCorrectionUncertainty* jecUnc){
  int nJets=0;
  vector <pair<float, int>> softLeptInJet; // pt, idx 
  for(edm::View<pat::Jet>::const_iterator ijet = jets->begin(); ijet!=jets->end();++ijet){
    nJets++;
    _jets_px.push_back( (float) ijet->px());
    _jets_py.push_back( (float) ijet->py());
    _jets_pz.push_back( (float) ijet->pz());
    _jets_e.push_back( (float) ijet->energy());    
    _jets_mT.push_back( (float) ijet->mt());
    _jets_Flavour.push_back(ijet->partonFlavour());
    _jets_HadronFlavour.push_back(ijet->hadronFlavour());
    _jets_PUJetID.push_back(ijet->userFloat("pileupJetId:fullDiscriminant"));
    float vtxPx = ijet->userFloat ("vtxPx");
    float vtxPy = ijet->userFloat ("vtxPy");
    _jets_vtxPt.  push_back(TMath::Sqrt(vtxPx*vtxPx + vtxPy*vtxPy));
    _jets_vtxMass.push_back(ijet->userFloat("vtxMass"));
    _jets_vtx3dL. push_back(ijet->userFloat("vtx3DVal"));
    _jets_vtxNtrk.push_back(ijet->userFloat("vtxNtracks"));
    _jets_vtx3deL.push_back(ijet->userFloat("vtx3DSig"));

    _bdiscr.push_back(ijet->bDiscriminator("jetBProbabilityBJetTags"));
    _bdiscr2.push_back(ijet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    //PF jet ID
    float NHF = ijet->neutralHadronEnergyFraction();
    float NEMF = ijet->neutralEmEnergyFraction();
    float CHF = ijet->chargedHadronEnergyFraction();
    float MUF = ijet->muonEnergyFraction();
    float CEMF = ijet->chargedEmEnergyFraction();
    int NumNeutralParticles =ijet->neutralMultiplicity();
    int chargedMult = ijet->chargedMultiplicity();
    int NumConst = ijet->chargedMultiplicity()+NumNeutralParticles;
    float CHM = ijet->chargedMultiplicity();
    float absjeta = fabs(ijet->eta());

    _jets_chEmEF .push_back(CEMF);  
    _jets_chHEF  .push_back(CHF); 
    _jets_nEmEF  .push_back(NEMF);
    _jets_nHEF   .push_back(NHF);
    _jets_chMult .push_back(chargedMult);  

    int jetid=0; 
    //PHYS14
    /*
    if((NHF<0.99 && NEMF<0.99 && NumConst>1 && MUF<0.8) && ((absjeta<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || absjeta>2.4)){
      jetid++;
      if( (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((absjeta<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || absjeta>2.4)  ) jetid++;
    }
    */
    //Spring15
    if(absjeta<=3.0){
      if((NHF<0.99 && NEMF<0.99 && NumConst>1) && ((absjeta<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || absjeta>2.4) ){
        jetid++;
        if( (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((absjeta<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || absjeta>2.4) ) {
          jetid++;
          if( (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((absjeta<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || absjeta>2.4)) jetid++;
        }
      }
    }else{
      if(NEMF<0.90 && NumNeutralParticles>10 ){
        jetid++;
        jetid++; //TIGHT and LOOSE are the same in this eta region
      }
    }    
    _jetID.push_back(jetid);
    float jecFactor = ijet->jecFactor("Uncorrected") ;
    float jetRawPt = jecFactor * ijet->pt();
    //float jetRawPt2 = ijet->pt() / jecFactor; // this is wrong
    _jets_rawPt.push_back ( jetRawPt );
    _jetrawf.push_back(jecFactor);
  
    // loop on jet contituents to retrieve info for b jet regression
    int nDau = ijet -> numberOfDaughters();
    //cout << "JET: " << (ijet - jets->begin()) << " N daught: " << nDau << endl;

    // TLorentzVector vJet (0,0,0,0);
    // vJet.SetPxPyPzE (ijet->px(), ijet->py(), ijet->pz(), ijet->energy());
    // TLorentzVector vDau (0,0,0,0); 
    // TLorentzVector vSum (0,0,0,0); 

    float leadTrackPt = 0.;
    softLeptInJet.clear();
    for (int iDau = 0; iDau < nDau; ++iDau)
    {
      // pdg id for packed pf candidates meaning is:
      // the particle charge and pdgId: 11, 13, 22 for ele/mu/gamma, 211 for charged hadrons, 130 for neutral hadrons, 1 and 2 for hadronic and em particles in HF. 
      const Candidate * dau = ijet->daughter(iDau);
      if (abs(dau->pdgId()) == 11 || abs(dau->pdgId()) == 13)
      {
          softLeptInJet.push_back( make_pair(dau->pt(), iDau) );
      }

      if (dau->charge() != 0 ) // tracks -> charged
      {
        float ptBuf = dau->pt();
        if (ptBuf > leadTrackPt) leadTrackPt = ptBuf;
      }
      // vDau.SetPxPyPzE (dau->px(), dau->py(), dau->pz(), dau->energy());
      // vSum += vDau;
      // cout << " - " << iDau << " pdg: " << dau->pdgId() << " pt: " << dau->pt() << " charge = " << dau->charge() << endl; 
    }

    //cout << " ## LEAD TRACK PT = " << leadTrackPt << endl;
    //cout << " ## jet eta: " << ijet->eta() << endl;
    _jets_leadTrackPt.push_back(leadTrackPt);
    float leptonPtRel = -1.;
    float leptonPt = -1.;
    float leptonDeltaR = -1.;
    int softLeptIdx = -1;
    if (softLeptInJet.size() > 0)
    {
      sort(softLeptInJet.begin(), softLeptInJet.end());
      softLeptIdx = softLeptInJet.back().second;
    }
    if (softLeptIdx >= 0)
    {
      const Candidate * dau = ijet->daughter(softLeptIdx);
      leptonPtRel = dau->pt() / ijet->pt() ;
      leptonPt = dau->pt() ;
      leptonDeltaR = deltaR(*dau, *ijet) ;
    }
    _jets_leptonPtRel .push_back (leptonPtRel);
    _jets_leptonPt    .push_back (leptonPt);
    _jets_leptonDeltaR.push_back (leptonDeltaR);

    //cout << "     --> jet pt, eta, phi: " << vJet.Pt() << " " << vJet.Eta() << " " << vJet.Phi() << endl;
    //cout << "     --> sum pt, eta, phi: " << vSum.Pt() << " " << vSum.Eta() << " " << vSum.Phi() << endl;
    //if (abs(ijet->hadronFlavour()) == 5 ) cout << "     ------------ THIS WAS A B JET ------------" << endl;
    //cout << "RAW pt: " << jetRawPt << " | " << jetRawPt2 << " --> " << vSum.Pt() << endl;

    jecUnc->setJetEta(ijet->eta());
    jecUnc->setJetPt(ijet->pt()); // here you must use the CORRECTED jet pt
    _jets_jecUnc.push_back(jecUnc->getUncertainty(true));
  }


  if ( theisMC )
  {
    edm::Handle<edm::View<reco::GenJet>> genJetHandle;
    event.getByLabel ("slimmedGenJets", genJetHandle);
    std::vector <const reco::GenJet*> ptrGenJets;
    unsigned int genJetSize = genJetHandle->size();
    for (unsigned int igj = 0; igj < genJetSize; igj++)
    {
      ptrGenJets.push_back ( &((*genJetHandle)[igj]) );
    }

    // now gather associated gen jets and save their address
    // relies on the fact that gen jets are filled in the same order as they appear in the collection
    std::vector <const reco::GenJet*>::iterator it;
    for(edm::View<pat::Jet>::const_iterator ijet = jets->begin(); ijet!=jets->end();++ijet)
    {
      const reco::GenJet * thisGenJet =  ijet->genJet ();
      int genindex = -1;
      it = std::find(ptrGenJets.begin(), ptrGenJets.end(), thisGenJet);
      if (it != ptrGenJets.end())
          genindex = std::distance (ptrGenJets.begin(), it);
      _jets_genjetIndex.push_back(genindex);
    }
  }

  return nJets;
}

//Fill all leptons (we keep them all for veto purposes
void HTauTauNtuplizer::FillSoftLeptons(const edm::View<reco::Candidate> *daus, const edm::Event& event,bool theFSR){
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;

  event.getByToken(triggerObjects_, triggerObjects);
  event.getByToken(triggerBits_, triggerBits);
  const edm::TriggerNames &names = event.triggerNames(*triggerBits);
  
  for(edm::View<reco::Candidate>::const_iterator daui = daus->begin(); daui!=daus->end();++daui){
    const reco::Candidate* cand = &(*daui);
    math::XYZTLorentzVector pfour = cand->p4();
    if(theFSR){
      const pat::PFParticle* fsr=0;
      double maxPT=-1;
      const PhotonPtrVector* gammas = userdatahelpers::getUserPhotons(cand);
      if (gammas!=0) {
        for (PhotonPtrVector::const_iterator g = gammas->begin();g!= gammas->end(); ++g) {
          //const pat::Photon* gamma = g->get();
          const pat::PFParticle* gamma = g->get();
          double pt = gamma->pt();
          if (pt>maxPT) {
            maxPT  = pt;
            fsr = gamma;
          }
        }
      }
      if (fsr!=0) {
        pfour+=fsr->p4();
      }
    } 
    
    _daughters_px.push_back( (float) pfour.X());
    _daughters_py.push_back( (float) pfour.Y());
    _daughters_pz.push_back( (float) pfour.Z());
    _daughters_e.push_back( (float) pfour.T());

    // gen info
    
    if (theisMC)
    {
        int iMatched = GetMatchedGen (cand, event);
        _daughters_genindex.push_back(iMatched);
    }

    //math::XYZTLorentzVector pfour(userdatahelpers::getUserFloat(cand,"genPx"),userdatahelpers::getUserFloat(cand,"genPy"),userdatahelpers::getUserFloat(cand,"genPz"),userdatahelpers::getUserFloat(cand,"genE"));
    //if(theisMC)_genDaughters.push_back(userdatahelpers::getUserFloat(cand,"fromH"));
    
    _softLeptons.push_back(cand);//This is needed also for FindCandIndex
    _pdgdau.push_back(cand->pdgId());
    _combreliso.push_back(userdatahelpers::getUserFloat(cand,"combRelIsoPF"));
    _dxy.push_back(userdatahelpers::getUserFloat(cand,"dxy"));
    _dz.push_back(userdatahelpers::getUserFloat(cand,"dz"));
    _SIP.push_back(userdatahelpers::getUserFloat(cand,"SIP"));
    int type = ParticleType::TAU;
    if(cand->isMuon()) type = ParticleType::MUON;
    else if(cand->isElectron()) type = ParticleType::ELECTRON;
    _particleType.push_back(type);
    
    // variables
    float discr=-1.;
    int muIDflag = 0;
    bool isgood = false;
    bool isele80=false;
    bool isele90=false;
    float elemva=-2;
    bool isconversionveto=false;
    int elemissinghits = 999;
    int decay=-1;
    float ieta=-1,superatvtx=-1,depositTracker=-1,depositEcal=-1,depositHcal=-1,SCeta=-999.;
    int decayModeFindingOldDMs=-1, decayModeFindingNewDMs=-1; // tau 13 TeV ID
    float byCombinedIsolationDeltaBetaCorrRaw3Hits=-1., chargedIsoPtSum=-1., neutralIsoPtSum=-1., puCorrPtSum=-1.; // tau 13 TeV RAW iso info
    int numChargedParticlesSignalCone=-1, numNeutralHadronsSignalCone=-1, numPhotonsSignalCone=-1, numParticlesSignalCone=-1, numChargedParticlesIsoCone=-1, numNeutralHadronsIsoCone=-1, numPhotonsIsoCone=-1, numParticlesIsoCone=-1;
    float leadChargedParticlePt=-1., trackRefPt=-1.;
    int typeOfMuon=0;
    float byIsolationMVA3oldDMwoLTraw=-1, byIsolationMVA3oldDMwLTraw=-1,  byIsolationMVA3newDMwoLTraw=-1,byIsolationMVA3newDMwLTraw=-1;
    Long64_t tauIDflag = 0;
    float   
    againstElectronMVA5category,
    againstElectronMVA5raw,
    byPileupWeightedIsolationRaw3Hits,
    footprintCorrection,
    neutralIsoPtSumWeight,
    photonPtSumOutsideSignalCone;

    if(type==ParticleType::MUON){
      muIDflag=userdatahelpers::getUserInt(cand,"muonID");
      discr = (float) muIDflag; // not really needed, will use the muonID branch in ntuples...
      if(userdatahelpers::getUserFloat(cand,"isPFMuon"))typeOfMuon |= 1 << 0;
      if(userdatahelpers::getUserFloat(cand,"isGlobalMuon"))typeOfMuon |= 1 << 1;
      if(userdatahelpers::getUserFloat(cand,"isTrackerMuon"))typeOfMuon |= 1 << 2;
      depositTracker=userdatahelpers::getUserFloat(cand,"DepositR03TrackerOfficial");
      depositEcal=userdatahelpers::getUserFloat(cand,"DepositR03Ecal");
      depositHcal=userdatahelpers::getUserFloat(cand,"DepositR03Hcal");
    }else if(type==ParticleType::ELECTRON){
      discr=userdatahelpers::getUserFloat(cand,"BDT");
      ieta=userdatahelpers::getUserFloat(cand,"sigmaIetaIeta");
      superatvtx=userdatahelpers::getUserFloat(cand,"deltaPhiSuperClusterTrackAtVtx");
      SCeta = userdatahelpers::getUserFloat(cand,"SCeta");
      if(userdatahelpers::getUserInt(cand,"isBDT") == 1)isgood=true;
      if(userdatahelpers::getUserInt(cand,"isEleID80") == 1) isele80=true;
      if(userdatahelpers::getUserInt(cand,"isEleID90") == 1) isele90=true;
      elemva=(userdatahelpers::getUserFloat(cand,"eleMVAvalue"));
      if(userdatahelpers::getUserInt(cand,"isConversionVeto") == 1)isconversionveto=true;
      elemissinghits = userdatahelpers::getUserInt(cand,"missingHit");
      //if(userdatahelpers::getUserInt(cand,"isCUT"))isgoodcut=true;
    }else if(type==ParticleType::TAU){
      discr=userdatahelpers::getUserFloat(cand,"HPSDiscriminator");
      decay = userdatahelpers::getUserFloat(cand,"decayMode");
      decayModeFindingOldDMs = userdatahelpers::getUserInt (cand, "decayModeFindingOldDMs");
      decayModeFindingNewDMs = userdatahelpers::getUserInt (cand, "decayModeFindingNewDMs");
      for (uint itau =0; itau<ntauIds; itau++){
        int id = userdatahelpers::getUserInt (cand,  tauIDStrings[itau]);
        if(id>0){
          tauIDflag |= (1 << itau);
          hTauIDs->Fill(id);
        }
      }
      againstElectronMVA5category = userdatahelpers::getUserFloat (cand, "againstElectronMVA5category");
      againstElectronMVA5raw = userdatahelpers::getUserFloat (cand, "againstElectronMVA5raw");
      byPileupWeightedIsolationRaw3Hits = userdatahelpers::getUserFloat (cand, "byPileupWeightedIsolationRaw3Hits");
      footprintCorrection = userdatahelpers::getUserFloat (cand, "footprintCorrection");
      neutralIsoPtSumWeight = userdatahelpers::getUserFloat (cand, "neutralIsoPtSumWeight");
      photonPtSumOutsideSignalCone = userdatahelpers::getUserFloat (cand, "photonPtSumOutsideSignalCone");

      byCombinedIsolationDeltaBetaCorrRaw3Hits = userdatahelpers::getUserFloat (cand, "byCombinedIsolationDeltaBetaCorrRaw3Hits");
      byIsolationMVA3oldDMwoLTraw=userdatahelpers::getUserFloat (cand, "byIsolationMVA3oldDMwoLTraw");
      byIsolationMVA3oldDMwLTraw=userdatahelpers::getUserFloat (cand, "byIsolationMVA3oldDMwLTraw");
      byIsolationMVA3newDMwoLTraw=userdatahelpers::getUserFloat (cand, "byIsolationMVA3newDMwoLTraw");
      byIsolationMVA3newDMwLTraw=userdatahelpers::getUserFloat (cand, "byIsolationMVA3newDMwLTraw");
      chargedIsoPtSum = userdatahelpers::getUserFloat (cand, "chargedIsoPtSum");
      neutralIsoPtSum = userdatahelpers::getUserFloat (cand, "neutralIsoPtSum");
      puCorrPtSum = userdatahelpers::getUserFloat (cand, "puCorrPtSum");
      numChargedParticlesSignalCone = userdatahelpers::getUserInt (cand, "numChargedParticlesSignalCone");
      numNeutralHadronsSignalCone = userdatahelpers::getUserInt (cand, "numNeutralHadronsSignalCone");
      numPhotonsSignalCone = userdatahelpers::getUserInt (cand, "numPhotonsSignalCone");
      numParticlesSignalCone = userdatahelpers::getUserInt (cand, "numParticlesSignalCone");
      numChargedParticlesIsoCone = userdatahelpers::getUserInt (cand, "numChargedParticlesIsoCone");
      numNeutralHadronsIsoCone = userdatahelpers::getUserInt (cand, "numNeutralHadronsIsoCone");
      numPhotonsIsoCone = userdatahelpers::getUserInt (cand, "numPhotonsIsoCone");
      numParticlesIsoCone = userdatahelpers::getUserInt (cand, "numParticlesIsoCone");
      leadChargedParticlePt = userdatahelpers::getUserFloat (cand, "leadChargedParticlePt");
      trackRefPt = userdatahelpers::getUserFloat (cand, "trackRefPt");

    }
    _discriminator.push_back(discr);
    _daughters_typeOfMuon.push_back(typeOfMuon);
    _daughters_muonID.push_back(muIDflag);
    _daughters_tauID.push_back(tauIDflag);
    _daughters_againstElectronMVA5category.push_back(againstElectronMVA5category);
    _daughters_againstElectronMVA5raw.push_back(againstElectronMVA5raw);
    _daughters_byPileupWeightedIsolationRaw3Hits.push_back(byPileupWeightedIsolationRaw3Hits);
    _daughters_footprintCorrection.push_back(footprintCorrection);
    _daughters_neutralIsoPtSumWeight.push_back(neutralIsoPtSumWeight);
    _daughters_photonPtSumOutsideSignalCone.push_back(photonPtSumOutsideSignalCone);

    _daughters_charge.push_back(cand->charge());
    _daughters_iseleBDT.push_back(isgood);
    _daughters_iseleWP80.push_back(isele80);
    _daughters_iseleWP90.push_back(isele90);
    _daughters_eleMVAnt.push_back(elemva);
    _daughters_passConversionVeto.push_back(isconversionveto);
    _daughters_eleMissingHits.push_back(elemissinghits);
    _daughters_iseleCUT.push_back(userdatahelpers::getUserInt(cand,"isCUT"));
    _decayType.push_back(decay);
    _daughters_IetaIeta.push_back(ieta);
    _daughters_deltaPhiSuperClusterTrackAtVtx.push_back(superatvtx);
    _daughters_SCeta.push_back(SCeta);
    _daughters_depositR03_tracker.push_back(depositTracker);
    _daughters_depositR03_ecal.push_back(depositEcal);
    _daughters_depositR03_hcal.push_back(depositHcal);
    _daughters_decayModeFindingOldDMs.push_back(decayModeFindingOldDMs);
    _daughters_decayModeFindingNewDMs.push_back(decayModeFindingNewDMs);
    _daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits.push_back(byCombinedIsolationDeltaBetaCorrRaw3Hits);
    _daughters_chargedIsoPtSum.push_back(chargedIsoPtSum);
    _daughters_neutralIsoPtSum.push_back(neutralIsoPtSum);
    _daughters_puCorrPtSum.push_back(puCorrPtSum);
    _daughters_byIsolationMVA3oldDMwoLTraw.push_back(byIsolationMVA3oldDMwoLTraw);
    _daughters_byIsolationMVA3oldDMwLTraw.push_back(byIsolationMVA3oldDMwLTraw);
    _daughters_byIsolationMVA3newDMwoLTraw.push_back(byIsolationMVA3newDMwoLTraw);
    _daughters_byIsolationMVA3newDMwLTraw.push_back(byIsolationMVA3newDMwLTraw);
    _daughters_numChargedParticlesSignalCone.push_back(numChargedParticlesSignalCone);
    _daughters_numNeutralHadronsSignalCone.push_back(numNeutralHadronsSignalCone);
    _daughters_numPhotonsSignalCone.push_back(numPhotonsSignalCone);
    _daughters_numParticlesSignalCone.push_back(numParticlesSignalCone);
    _daughters_numChargedParticlesIsoCone.push_back(numChargedParticlesIsoCone);
    _daughters_numNeutralHadronsIsoCone.push_back(numNeutralHadronsIsoCone);
    _daughters_numPhotonsIsoCone.push_back(numPhotonsIsoCone);
    _daughters_numParticlesIsoCone.push_back(numParticlesIsoCone);
    _daughters_leadChargedParticlePt.push_back(leadChargedParticlePt);
    _daughters_trackRefPt.push_back(trackRefPt);

    //TRIGGER MATCHING
    int LFtriggerbit=0,L3triggerbit=0,filterFired=0;
    int triggertypeIsGood = 0;
    float hltpt=0;
    for (pat::TriggerObjectStandAlone obj : *triggerObjects) { 
      //check if the trigger object matches cand
      //bool matchCand = false;
      //if(type == ParticleType::TAU && 
      bool triggerType=false;
      if(deltaR2(obj,*cand)<0.25){
        if (type==ParticleType::TAU && (obj.hasTriggerObjectType(trigger::TriggerTau)|| obj.hasTriggerObjectType(trigger::TriggerL1TauJet)))triggerType=true;
        if (type==ParticleType::ELECTRON && (obj.hasTriggerObjectType(trigger::TriggerElectron) || obj.hasTriggerObjectType(trigger::TriggerPhoton)))triggerType=true;
        if (type==ParticleType::MUON && (obj.hasTriggerObjectType(trigger::TriggerMuon)))triggerType=true;
        //triggerhelper myTriggerHelper;
        //check fired paths
        obj.unpackPathNames(names);
        std::vector<std::string> pathNamesAll  = obj.pathNames(false);
        std::vector<std::string> pathNamesLast = obj.pathNames(true);
        for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
          bool isLF   = obj.hasPathName( pathNamesAll[h], true, false ); 
          bool isL3   = obj.hasPathName( pathNamesAll[h], false, true );
	  if (type==ParticleType::MUON && obj.hasFilterLabel("hltL3crIsoL1sSingleMu16erL1f0L2f10QL3f17QL3trkIsoFiltered0p09") && (obj.hasPathName("HLT_IsoMu17_eta2p1_v1", false, true ) or obj.hasPathName("HLT_IsoMu17_eta2p1_v1", true, false )))
	     hltpt = (float) obj.pt();
          if (type==ParticleType::ELECTRON && obj.hasFilterLabel("hltSingleEle22WP75GsfTrackIsoFilter") && (obj.hasPathName("HLT_Ele22_eta2p1_WP75_Gsf_v1", false, true ) or obj.hasPathName("HLT_Ele22_eta2p1_WP75_Gsf_v1", true, false )))
             hltpt = (float) obj.pt();
          Long64_t triggerbit = myTriggerHelper->FindTriggerNumber(pathNamesAll[h],true);
                    
          /*
          for (int i = 0; i < myTriggerHelper->GetNTriggers(); i++)
          {
            TString nameTr = myTriggerHelper->printTriggerName(i);
            cout << i << ") " << nameTr << " - " << myTriggerHelper->FindTriggerNumber(nameTr) << endl;
          }
          cout << " =============================== " << endl;
          */
          //int triggerbit = _triggerbit;
          if(triggerbit>=0){ // get the number (position) of this trigger in the _triggerbit branch output (yes, names are too similar!)
            triggerMapper map = myTriggerHelper->GetTriggerMap(pathNamesAll[h]);
            bool isfilterGood = true;
            if(type==ParticleType::TAU){
              for(int ifilt=0;ifilt<map.GetNfiltersleg2();ifilt++){
                if(! obj.hasFilterLabel(map.Getfilter(false,ifilt).c_str()))isfilterGood=false;
              }
            }else if(type==ParticleType::ELECTRON){
              for(int ifilt=0;ifilt<map.GetNfiltersleg1();ifilt++){
                if(! obj.hasFilterLabel(map.Getfilter(true,ifilt).c_str()))isfilterGood=false;
              }
            }else{//muons
              if(map.GetTriggerChannel()==triggerMapper::kemu){
                for(int ifilt=0;ifilt<map.GetNfiltersleg1();ifilt++){
                  if(! obj.hasFilterLabel(map.Getfilter(false,ifilt).c_str()))isfilterGood=false;
                }
              }else{
                for(int ifilt=0;ifilt<map.GetNfiltersleg2();ifilt++){
                  if(! obj.hasFilterLabel(map.Getfilter(true,ifilt).c_str()))isfilterGood=false;
                }
              }
            }
          //_isFilterFiredLast;
            if(isfilterGood)filterFired |= 1 <<triggerbit;
            if(isLF)LFtriggerbit |= 1 <<triggerbit;
            if(isL3)L3triggerbit |= 1 <<triggerbit;
            if (triggerType) triggertypeIsGood |= 1 << triggerbit;
          }
        }
      }
    }
    _daughters_isGoodTriggerType.push_back(triggertypeIsGood);
    _daughters_FilterFired.push_back(filterFired);
    _daughters_L3FilterFired.push_back(LFtriggerbit);
    _daughters_L3FilterFiredLast.push_back(L3triggerbit);    
    _daughters_HLTpt.push_back(hltpt);
  }
}

/*
void HTauTauNtuplizer::FillbQuarks(const edm::Event& event){
  edm::Handle<edm::View<pat::GenericParticle>>candHandle;
  event.getByLabel("bQuarks",candHandle);
  const edm::View<pat::GenericParticle>* bs = candHandle.product();
  for(edm::View<pat::GenericParticle>::const_iterator ib = bs->begin(); ib!=bs->end();++ib){
    const pat::GenericParticle* cand = &(*ib);
    _bquarks_px.push_back( (float) cand->px());
    _bquarks_py.push_back( (float) cand->py());
    _bquarks_pz.push_back( (float) cand->px());
    _bquarks_e.push_back( (float) cand->energy());
    _bquarks_pdg.push_back( (int) cand->pdgId());
    _bmotmass.push_back(userdatahelpers::getUserFloat(cand,"motHmass"));
  }

  //Retrieve Generated H (there can be more than 1!  
  Handle<edm::View<reco::GenParticle> > prunedHandle;
  event.getByLabel("prunedGenParticles", prunedHandle);
  for(unsigned int ipruned = 0; ipruned< prunedHandle->size(); ++ipruned){
     const GenParticle *packed =&(*prunedHandle)[ipruned];
     int pdgh = packed->pdgId();
     if(abs(pdgh)==25){
        // avoid Higgs clones, save only the one not going into another H (end of showering process)
        bool isLast = true;
        for (unsigned int iDau = 0; iDau < packed->numberOfDaughters(); iDau++)
        {
            const Candidate* Dau = packed->daughter( iDau );
            if (Dau->pdgId() == pdgh)
            {
                isLast = false;
                break;
            }
        }
   
        if (isLast)
        {     
            _genH_px.push_back(packed->px());          
            _genH_py.push_back(packed->py());          
            _genH_pz.push_back(packed->pz());          
            _genH_e.push_back(packed->energy());          
        }
     }        
  } 
}
*/

void HTauTauNtuplizer::FillGenInfo(const edm::Event& event)
{
    edm::Handle<edm::View<pat::GenericParticle>>candHandle;
    event.getByLabel ("genInfo", candHandle);
    const edm::View<pat::GenericParticle>* gens = candHandle.product();
    for(edm::View<pat::GenericParticle>::const_iterator igen = gens->begin(); igen!=gens->end(); ++igen)
    {
        // fill gen particle branches
        _genpart_px.push_back(igen->px());
        _genpart_py.push_back(igen->py());
        _genpart_pz.push_back(igen->pz());
        _genpart_e.push_back(igen->energy());
        _genpart_pdg.push_back(igen->pdgId());
        _genpart_status.push_back(igen->status());
        
        int HMIndex = -1;
        int MSSMHMIndex = -1;
        int TopMIndex = -1;
        int TauMIndex = -1;
        int ZMIndex = -1;
        int WMIndex = -1;
        int bMIndex = -1;
        int HZDecayMode = -1;
        int TopDecayMode = -1;
        int WDecayMode = -1;
        int TauGenDecayMode = -1;
        
        if (igen->hasUserInt("HMothIndex"))   HMIndex = igen->userInt("HMothIndex");
        if (igen->hasUserInt("MSSMHMothIndex"))   MSSMHMIndex = igen->userInt("MSSMHMothIndex");
        if (igen->hasUserInt("TopMothIndex")) TopMIndex = igen->userInt("TopMothIndex");
        if (igen->hasUserInt("TauMothIndex")) TauMIndex = igen->userInt("TauMothIndex");
        if (igen->hasUserInt("ZMothIndex"))   ZMIndex = igen->userInt("ZMothIndex");
        if (igen->hasUserInt("WMothIndex")) WMIndex = igen->userInt("WMothIndex");
        if (igen->hasUserInt("bMothIndex")) bMIndex = igen->userInt("bMothIndex");
        if (igen->hasUserInt("HZDecayMode"))   HZDecayMode = igen->userInt("HZDecayMode");
        if (igen->hasUserInt("TopDecayMode"))   TopDecayMode = igen->userInt("TopDecayMode");
        if (igen->hasUserInt("WDecayMode"))   WDecayMode = igen->userInt("WDecayMode");
        if (igen->hasUserInt("tauGenDecayMode"))   TauGenDecayMode = igen->userInt("tauGenDecayMode");
        
        _genpart_HMothInd.push_back(HMIndex);
        _genpart_MSSMHMothInd.push_back(MSSMHMIndex);
        _genpart_TopMothInd.push_back(TopMIndex);
        _genpart_TauMothInd.push_back(TauMIndex);
        _genpart_ZMothInd.push_back(ZMIndex);
        _genpart_WMothInd.push_back(WMIndex);
        _genpart_bMothInd.push_back(bMIndex);
        _genpart_HZDecayMode.push_back(HZDecayMode);
        _genpart_TopDecayMode.push_back(TopDecayMode);
        _genpart_WDecayMode.push_back(WDecayMode);
        _genpart_TauGenDecayMode.push_back(TauGenDecayMode);
        
        //const pat::GenericParticle* genClone = &(*igen);
        //int flags = CreateFlagsWord (genClone);
        int flags = igen -> userInt ("generalGenFlags");
        _genpart_flags.push_back(flags);
    }
}


void HTauTauNtuplizer::FillGenJetInfo(const edm::Event& event)
{
    edm::Handle<edm::View<reco::GenJet>> genJetHandle;
    event.getByLabel ("slimmedGenJets", genJetHandle);
    unsigned int genJetSize = genJetHandle->size();

    // to retrieve gen jet flavour from matched gen jet
    edm::Handle<edm::View<pat::Jet>> patjetHandle;
    event.getByLabel("jets", patjetHandle);
    unsigned int jetSize = patjetHandle->size();

    for (unsigned int igj = 0; igj < genJetSize; igj++)
    {
      const reco::GenJet& genJet = (*genJetHandle)[igj];
      _genjet_px.push_back ( genJet.px() );
      _genjet_py.push_back ( genJet.py() );
      _genjet_pz.push_back ( genJet.pz() );
      _genjet_e .push_back ( genJet.energy() );

      // jet flavour
      int partFlav = -999;
      int hadrFlav = -999;

      for (unsigned int ijet = 0; ijet < jetSize; ijet++)
      {
        const pat::Jet& patjet = (*patjetHandle)[ijet];
        const reco::GenJet * thismatchedGenJet =  patjet.genJet();

        if (thismatchedGenJet == &genJet)
        {
          // no error, checked :-)
          //if (partFlav != -999 && patjet.partonFlavour() != partFlav) cout << igj << " MISMATCH! Part flav: " << partFlav << " " << patjet.partonFlavour() << endl;
          //if (hadrFlav != -999 && patjet.hadronFlavour() != hadrFlav) cout << igj << " MISMATCH! Hadr flav: " << hadrFlav << " " << patjet.hadronFlavour() << endl;
          partFlav = patjet.partonFlavour();
          hadrFlav = patjet.hadronFlavour();
          break;
        }      
      }

      _genjet_partonFlavour.push_back(partFlav);
      _genjet_hadronFlavour.push_back(hadrFlav);
    }

    return;

}

// return index of gen matched to reco lepton, and -1 if not existing or not found
int HTauTauNtuplizer::GetMatchedGen (const reco::Candidate* genL, const edm::Event& event)
{
    //cout.precision(15); // just to check real precision
    
    edm::Handle<edm::View<pat::GenericParticle>>candHandle;
    event.getByLabel ("genInfo", candHandle);
        
    int index = -1;
    int status = userdatahelpers::getUserInt(genL,"status");
    int id = userdatahelpers::getUserInt(genL,"id");
    if (status == 99999 && id == 99999) return -1; // default values in *Filler.cc when no matched gen found
        
    // get matched particle by looking at pdgId, px, py, pz, e --> not the most elegant, but should work
    for (unsigned int iGen = 0; iGen < candHandle->size(); iGen++)
    {
        const pat::GenericParticle& genP = (*candHandle)[iGen];
        // first check on pdgId only, then status to remove most of checks with a minimum amount of comparison
        if (id == genP.pdgId())
        {
            if (status == genP.status())
            {
                
                float px = userdatahelpers::getUserFloat(genL,"genPx");
                float py = userdatahelpers::getUserFloat(genL,"genPy");
                float pz = userdatahelpers::getUserFloat(genL,"genPz");
                float e = userdatahelpers::getUserFloat(genL,"genE");            
                //cout << "     ==> I'm in with " << fixed<< px << " " << fixed<<py << " " << fixed<<pz << " " << fixed<<e << "  ||| " << fixed<<genP.px() << " " << fixed<<genP.py() << " " << fixed<<genP.pz() << " " << fixed<<genP.energy() << endl;
                // cast to float helps in the EPSIL comparison, but for safety EPSIL = 1.e-5 is used (seems reasonably small for the value range we use in 0.001 -- 1000) with "meaningful" values O(10) or larger
                if ( fabs((float)genP.px() - px) < EPSIL && fabs((float)genP.py() - py) < EPSIL && fabs((float)genP.pz() - pz) < EPSIL && fabs((float)genP.energy() - e) < EPSIL )
                {
                    index = iGen;
                    break; // found, no other comparisons needed
                }
            }
        }
    }
    
    return index;
}

// not used anymore, all info already stored in genFiller
// int HTauTauNtuplizer::CreateFlagsWord (const pat::GenericParticle* part)
// {
//     int flag = 0;
//     
//     if (part->hasUserInt("HMothIndex"))      flag |= (1 << static_cast<int> (GenFlags::fromH));
//     if (part->hasUserInt("TopMothIndex"))    flag |= (1 << static_cast<int> (GenFlags::fromTop));
//     if (part->hasUserInt("TauMothIndex"))    flag |= (1 << static_cast<int> (GenFlags::fromTau));
//     if (part->hasUserInt("ZMothIndex"))      flag |= (1 << static_cast<int> (GenFlags::fromZ));
//     
//     /*
//     // H/Z decau --> was changed into a dedicated branch
//     if (part->hasUserInt("HZDecayMode"))
//     {
//         int decayMode = part->userInt ("HZDecayMode");
//         int bitToUse = genFlagPosMap_.at(decayMode);
//         flag |= (1 << bitToUse);
//     }
//     */
//     return flag;
// }


void HTauTauNtuplizer::endJob(){
  hCounter->SetBinContent(1,Nevt_Gen);
  hCounter->SetBinContent(2,Nevt_PassTrigger);
  hCounter->SetBinContent(3,Npairs);

  hCounter->GetXaxis()->SetBinLabel(1,"Nevt_Gen");
  hCounter->GetXaxis()->SetBinLabel(2,"Nevt_PassTrigger");
  hCounter->GetXaxis()->SetBinLabel(3,"Npairs");

  for(int i=0;i<myTriggerHelper->GetNTriggers();i++){
    hCounter->GetXaxis()->SetBinLabel(i+4,(myTriggerHelper->printTriggerName(i)).c_str());
  }
    for(int i=1;i<=ntauIds;i++){
    hTauIDs->GetXaxis()->SetBinLabel(i,tauIDStrings[i-1].Data());
  }

}


void HTauTauNtuplizer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){
  
  Bool_t changedConfig = false;
 
  //if(!hltConfig_.init(iRun, iSetup, triggerResultsLabel.process(), changedConfig)){
  if(!hltConfig_.init(iRun, iSetup, processName.process(), changedConfig)){
    edm::LogError("HLTMatchingFilter") << "Initialization of HLTConfigProvider failed!!"; 
    return;
  }  

  if(changedConfig || foundPaths.size()==0){
    //cout<<"The present menu is "<<hltConfig.tableName()<<endl;
    indexOfPath.clear();
    foundPaths.clear();
    //for(size_t i=0; i<triggerPaths.size(); i++){
    // bool foundThisPath = false;
    for(size_t j=0; j<hltConfig_.triggerNames().size(); j++){
      string pathName = hltConfig_.triggerNames()[j];
      //TString tempo= hltConfig_.triggerNames()[j];
      //printf("%s\n",tempo.Data());
      //if(pathName==triggerPaths[i]){
      //foundThisPath = true;
      indexOfPath.push_back(j);
      foundPaths.push_back(pathName);

      cout << j << " - TTT: " << pathName << endl;
	  //	  edm::LogInfo("AnalyzeRates")<<"Added path "<<pathName<<" to foundPaths";
    } 
  }
  
}


void HTauTauNtuplizer::endRun(edm::Run const&, edm::EventSetup const&){
}
void HTauTauNtuplizer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){
}
void HTauTauNtuplizer::endLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const& iSetup)
{
  // Total number of events is the sum of the events in each of these luminosity blocks
  edm::Handle<edm::MergeableCounter> nEventsTotalCounter;
  iLumi.getByLabel("nEventsTotal", nEventsTotalCounter);
  Nevt_Gen += nEventsTotalCounter->value;

  edm::Handle<edm::MergeableCounter> nEventsPassTrigCounter;
  iLumi.getByLabel("nEventsPassTrigger", nEventsPassTrigCounter);
  Nevt_PassTrigger += nEventsPassTrigCounter->value;
}


bool HTauTauNtuplizer::CompareLegs(const reco::Candidate *i, const reco::Candidate *j){
  int iType=2,jType=2;
  
  if(i->isElectron())iType=1;
  else if(i->isMuon())iType=0;
  
  if(j->isElectron())jType=1;
  else if(j->isMuon())jType=0;
  
  if(iType>jType) return false;
  else if(iType==jType && i->pt()<j->pt()) return false;
  
  return true;
}
// implements operator "<" (return i < j)
bool HTauTauNtuplizer::ComparePairsbyIso(pat::CompositeCandidate i, pat::CompositeCandidate j){
  
  //First criteria: OS<SS
  //if ( j.charge()==0 && i.charge()!=0) return false;
  //else if ( i.charge()==0 && j.charge()!=0) return true;

  //Second criteria: ISO
  float isoi=999,isoj=999;
  int cand1j=-1,cand1i=-1;

  if(CompareLegs(i.daughter(0),i.daughter(1)))cand1i=0;
  else cand1i=1;
  if(CompareLegs(j.daughter(0),j.daughter(1)))cand1j=0;
  else cand1j=1;

  //step 1, l;eg 1 ISO
  //byCombinedIsolationDeltaBetaCorrRaw3Hits
  isoi=userdatahelpers::getUserFloat(i.daughter(cand1i),"combRelIsoPF");
  isoj=userdatahelpers::getUserFloat(j.daughter(cand1j),"combRelIsoPF");
  if (!i.daughter(cand1i)->isMuon() && !i.daughter(cand1i)->isElectron()) isoi=userdatahelpers::getUserFloat(i.daughter(cand1i),"byCombinedIsolationDeltaBetaCorrRaw3Hits");
  if (!j.daughter(cand1j)->isMuon() && !j.daughter(cand1j)->isElectron()) isoj=userdatahelpers::getUserFloat(j.daughter(cand1j),"byCombinedIsolationDeltaBetaCorrRaw3Hits");

  if (isoi<isoj)return true;
  else if(isoi>isoj)return false;

  //step 2, leg 1 Pt
  if(i.daughter(cand1i)->pt()>j.daughter(cand1j)->pt()) return true;
  else if(i.daughter(cand1i)->pt()<j.daughter(cand1j)->pt()) return false;

  //step 3, leg 2 ISO
  isoi=userdatahelpers::getUserFloat(i.daughter(1-cand1i),"combRelIsoPF");
  isoj=userdatahelpers::getUserFloat(j.daughter(1-cand1j),"combRelIsoPF");
  if (!i.daughter(1-cand1i)->isMuon() && !i.daughter(1-cand1i)->isElectron()) isoi=userdatahelpers::getUserFloat(i.daughter(1-cand1i),"byCombinedIsolationDeltaBetaCorrRaw3Hits");
  if (!j.daughter(1-cand1j)->isMuon() && !j.daughter(1-cand1j)->isElectron()) isoj=userdatahelpers::getUserFloat(j.daughter(1-cand1j),"byCombinedIsolationDeltaBetaCorrRaw3Hits");

  if (isoi<isoj)return true;
  else if(isoi>isoj)return false;

  //step 4, leg 2 Pt
  if(i.daughter(1-cand1i)->pt()>j.daughter(1-cand1j)->pt()) return true;
  //else if(i.daughter(1-cand1i)->pt()<j.daughter(1-cand1j)->pt()) return false;

  return false;

}

bool HTauTauNtuplizer::ComparePairsbyPt(pat::CompositeCandidate i, pat::CompositeCandidate j){
  
  //First criteria: OS<SS
  //if ( j.charge()==0 && i.charge()!=0) return false;
  //else if ( i.charge()==0 && j.charge()!=0) return true;

  //Second criteria: legs pt
  if(i.daughter(0)->pt()+i.daughter(1)->pt()<j.daughter(0)->pt()+j.daughter(1)->pt()) return false;
  if(i.daughter(0)->pt()+i.daughter(1)->pt()>j.daughter(0)->pt()+j.daughter(1)->pt()) return true; 
  
  //Protection for duplicated taus, damn it!
  int iType=0,jType=0;
  for(int idau=0;idau<2;idau++){
    if(i.daughter(idau)->isElectron())iType+=0;
    else if(i.daughter(idau)->isMuon())iType+=1;
    else iType+=2;
    if(j.daughter(idau)->isElectron())jType+=0;
    else if(j.daughter(idau)->isMuon())jType+=1;
    else jType+=2;
  }
  
  return (iType<jType);
  
  //I need a criteria for where to put pairs wo taus and how to deal with tautau candidates
  
  //third criteria: are there taus?
  int ilegtau=-1,jlegtau=-1;
  if(!i.daughter(0)->isMuon() && !i.daughter(0)->isElectron())ilegtau=0;
  if(!j.daughter(0)->isMuon() && !j.daughter(0)->isElectron())jlegtau=0;
 
  if(!i.daughter(1)->isMuon() && !i.daughter(1)->isElectron()){
    if(ilegtau==-1 || i.daughter(1)->pt()>i.daughter(0)->pt())ilegtau=1;}
  if(!j.daughter(1)->isMuon() && !j.daughter(1)->isElectron()){
    if(ilegtau==-1 || i.daughter(1)->pt()>i.daughter(0)->pt())jlegtau=1;}

  
  if(ilegtau==-1 && jlegtau>-1) return false; //i has no tau leptons, j has
  else if(ilegtau==-1 && jlegtau == -1) return true; //no tau leptons in neither pair, leave as it is
  
  //fourth criteria: Iso x Legpt
  if(userdatahelpers::getUserFloat(i.daughter(ilegtau),"byCombinedIsolationDeltaBetaCorrRaw3Hits")*i.daughter(fabs(ilegtau-1))->pt()>userdatahelpers::getUserFloat(j.daughter(jlegtau),"byCombinedIsolationDeltaBetaCorrRaw3Hits")*j.daughter(fabs(jlegtau-1))->pt()) return false;
  
  //fifth criteria: Iso (cut based)
  if(userdatahelpers::getUserFloat(i.daughter(ilegtau),"combRelIsoPF")>userdatahelpers::getUserFloat(j.daughter(jlegtau),"combRelIsoPF")) return false;
  
  //sixth criteria: ISO (MVA)
  if(userdatahelpers::getUserFloat(i.daughter(ilegtau),"byCombinedIsolationDeltaBetaCorrRaw3Hits")*userdatahelpers::getUserFloat(i.daughter(fabs(ilegtau-1)),"combRelIsoPF")>userdatahelpers::getUserFloat(j.daughter(jlegtau),"byCombinedIsolationDeltaBetaCorrRaw3Hits")*userdatahelpers::getUserFloat(j.daughter(fabs(jlegtau-1)),"combRelIsoPF")) return false;

  return true;
}


float HTauTauNtuplizer::ComputeMT (math::XYZTLorentzVector visP4, float METx, float METy)
{
    float MET = TMath::Sqrt (METx*METx + METy*METy);
    math::XYZTLorentzVector METP4 (METx, METy, 0, MET);
    
    float scalSum = MET + visP4.pt();
    math::XYZTLorentzVector vecSum (visP4);
    vecSum += METP4;
    float vecSumPt = vecSum.pt();
    
    return TMath::Sqrt (scalSum*scalSum - vecSumPt*vecSumPt);
}

// // ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
// void HTauTauNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
//   //The following says we do not know what parameters are allowed so do no validation
//   // Please change this to state exactly what you do use, even if it is no parameters
//   edm::ParameterSetDescription desc;
//   desc.setUnknown();
//   descriptions.addDefault(desc);
// }

//define this as a plug-in
DEFINE_FWK_MODULE(HTauTauNtuplizer);
