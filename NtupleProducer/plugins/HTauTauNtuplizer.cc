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


//#include "TLorentzVector.h"

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
  int FillJet(const edm::View<pat::Jet>* jet);
  void FillSoftLeptons(const edm::View<reco::Candidate> *dauhandler, const edm::Event& event, bool theFSR);
  //void FillbQuarks(const edm::Event&);
  void FillGenInfo(const edm::Event&);
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
  triggerhelper* myTriggerHelper;

  PUReweight reweight;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;

  //flags
  static const int nOutVars =14;
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
  Int_t _triggerbit;
  Int_t _metfilterbit;
  Float_t _met;
  Float_t _metphi;
  Float_t _MC_weight;
  Float_t _aMCatNLOweight;
  Int_t _npv;
  Int_t _npu;
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
  
  //std::vector<math::XYZTLorentzVector> _daughter2;

  //Mothers output variables
  std::vector<Int_t> _indexDau1;
  std::vector<Int_t> _indexDau2;
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
  std::vector<Float_t> _daughters_IetaIeta;
  std::vector<Float_t> _daughters_deltaPhiSuperClusterTrackAtVtx;
  std::vector<Float_t> _daughters_SCeta;
  std::vector<Float_t> _daughters_depositR03_tracker;
  std::vector<Float_t> _daughters_depositR03_ecal;
  std::vector<Float_t> _daughters_depositR03_hcal;
  std::vector<Int_t> _daughters_decayModeFindingOldDMs;
  std::vector<Int_t> _daughters_decayModeFindingNewDMs;
  std::vector<Int_t> _daughters_byLooseCombinedIsolationDeltaBetaCorr3Hits;
  std::vector<Int_t> _daughters_byMediumCombinedIsolationDeltaBetaCorr3Hits;
  std::vector<Int_t> _daughters_byTightCombinedIsolationDeltaBetaCorr3Hits;
  std::vector<Float_t> _daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits;
  std::vector<Float_t> _daughters_byIsolationMVA3oldDMwoLTraw;
  std::vector<Float_t> _daughters_byIsolationMVA3oldDMwLTraw;
  std::vector<Float_t> _daughters_byIsolationMVA3newDMwoLTraw;
  std::vector<Float_t> _daughters_byIsolationMVA3newDMwLTraw;
  std::vector<Float_t> _daughters_chargedIsoPtSum;
  std::vector<Float_t> _daughters_neutralIsoPtSum;
  std::vector<Float_t> _daughters_puCorrPtSum;
  std::vector<Int_t> _daughters_againstMuonLoose3;
  std::vector<Int_t> _daughters_againstMuonTight3;
  std::vector<Int_t> _daughters_againstElectronVLooseMVA5;
  std::vector<Int_t> _daughters_againstElectronLooseMVA5;
  std::vector<Int_t> _daughters_againstElectronMediumMVA5;
  std::vector<Int_t> _daughters_againstElectronTightMVA5;
  std::vector<Int_t> _daughters_againstElectronVTightMVA5;
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
  std::vector<Float_t> _jets_PUJetID;
  std::vector<Int_t> _jets_Flavour;
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
  _daughters_byLooseCombinedIsolationDeltaBetaCorr3Hits.clear();
  _daughters_byMediumCombinedIsolationDeltaBetaCorr3Hits.clear();
  _daughters_byTightCombinedIsolationDeltaBetaCorr3Hits.clear();
  _daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits.clear();
  _daughters_byIsolationMVA3oldDMwoLTraw.clear();
  _daughters_byIsolationMVA3oldDMwLTraw.clear();
  _daughters_byIsolationMVA3newDMwoLTraw.clear();
  _daughters_byIsolationMVA3newDMwLTraw.clear();
  _daughters_chargedIsoPtSum.clear();
  _daughters_neutralIsoPtSum.clear();
  _daughters_puCorrPtSum.clear();
  _daughters_againstMuonLoose3.clear();
  _daughters_againstMuonTight3.clear();
  _daughters_againstElectronVLooseMVA5.clear();
  _daughters_againstElectronLooseMVA5.clear();
  _daughters_againstElectronMediumMVA5.clear();
  _daughters_againstElectronTightMVA5.clear();
  _daughters_againstElectronVTightMVA5.clear();
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
  _npu=0;
  _PUReweight=0.;
  _rho=0;
  _nup=-999;

//  _jets.clear();
  _jets_px.clear();
  _jets_py.clear();
  _jets_pz.clear();
  _jets_e.clear();
  _jets_PUJetID.clear();
  _jets_Flavour.clear();
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

  //Branches
  myTree->Branch("EventNumber",&_indexevents,"EventNumber/I");
  myTree->Branch("RunNumber",&_runNumber,"RunNumber/I");
  myTree->Branch("lumi",&_lumi,"lumi/I");
  myTree->Branch("triggerbit",&_triggerbit,"triggerbit/I");
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

    myTree->Branch("daughters_genindex",&_daughters_genindex);
    myTree->Branch("MC_weight",&_MC_weight,"MC_weight/F");
    myTree->Branch("aMCatNLOweight",&_aMCatNLOweight,"aMCatNLOweight/F");    
    myTree->Branch("genpart_px", &_genpart_px);
    myTree->Branch("genpart_py", &_genpart_py);
    myTree->Branch("genpart_pz", &_genpart_pz);
    myTree->Branch("genpart_e", &_genpart_e);
    myTree->Branch("genpart_pdg", &_genpart_pdg);
    myTree->Branch("genpart_status", &_genpart_status);
    //myTree->Branch("genpart_mothInd", _genpart_mothInd);
    myTree->Branch("genpart_HMothInd", &_genpart_HMothInd);
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
  myTree->Branch("combreliso",& _combreliso);
  myTree->Branch("daughters_IetaIeta",&_daughters_IetaIeta);
  myTree->Branch("daughters_deltaPhiSuperClusterTrackAtVtx",&_daughters_deltaPhiSuperClusterTrackAtVtx);
  myTree->Branch("daughters_SCeta",&_daughters_SCeta);
  myTree->Branch("daughters_depositR03_tracker",&_daughters_depositR03_tracker);
  myTree->Branch("daughters_depositR03_ecal",&_daughters_depositR03_ecal);
  myTree->Branch("daughters_depositR03_hcal",&_daughters_depositR03_hcal);
  myTree->Branch("daughters_decayModeFindingOldDMs", &_daughters_decayModeFindingOldDMs);
  myTree->Branch("daughters_decayModeFindingNewDMs", &_daughters_decayModeFindingNewDMs);
  myTree->Branch("daughters_byLooseCombinedIsolationDeltaBetaCorr3Hits", &_daughters_byLooseCombinedIsolationDeltaBetaCorr3Hits);
  myTree->Branch("daughters_byMediumCombinedIsolationDeltaBetaCorr3Hits", &_daughters_byMediumCombinedIsolationDeltaBetaCorr3Hits);
  myTree->Branch("daughters_byTightCombinedIsolationDeltaBetaCorr3Hits", &_daughters_byTightCombinedIsolationDeltaBetaCorr3Hits);
  myTree->Branch("daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits", &_daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits);
  myTree->Branch("daughters_byIsolationMVA3oldDMwoLTraw",&_daughters_byIsolationMVA3oldDMwoLTraw);
  myTree->Branch("daughters_byIsolationMVA3oldDMwLTraw",&_daughters_byIsolationMVA3oldDMwLTraw);
  myTree->Branch("daughters_byIsolationMVA3newDMwoLTraw",&_daughters_byIsolationMVA3newDMwoLTraw);
  myTree->Branch("daughters_byIsolationMVA3newDMwLTraw",&_daughters_byIsolationMVA3newDMwLTraw);
  myTree->Branch("daughters_chargedIsoPtSum", &_daughters_chargedIsoPtSum);
  myTree->Branch("daughters_neutralIsoPtSum", &_daughters_neutralIsoPtSum);
  myTree->Branch("daughters_puCorrPtSum", &_daughters_puCorrPtSum);
  myTree->Branch("daughters_againstMuonLoose3", &_daughters_againstMuonLoose3);
  myTree->Branch("daughters_againstMuonTight3", &_daughters_againstMuonTight3);
  myTree->Branch("daughters_againstElectronVLooseMVA5", &_daughters_againstElectronVLooseMVA5);
  myTree->Branch("daughters_againstElectronLooseMVA5", &_daughters_againstElectronLooseMVA5);
  myTree->Branch("daughters_againstElectronMediumMVA5", &_daughters_againstElectronMediumMVA5);
  myTree->Branch("daughters_againstElectronTightMVA5", &_daughters_againstElectronTightMVA5);
  myTree->Branch("daughters_againstElectronVTightMVA5", &_daughters_againstElectronVTightMVA5);
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
  
  myTree->Branch("JetsNumber",&_numberOfJets,"JetsNumber/I");
  myTree->Branch("jets_px",&_jets_px);
  myTree->Branch("jets_py",&_jets_py);
  myTree->Branch("jets_pz",&_jets_pz);
  myTree->Branch("jets_e",&_jets_e);
  myTree->Branch("jets_Flavour",&_jets_Flavour);
  myTree->Branch("jets_PUJetID",&_jets_PUJetID);
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
   if (theisMC) {
    Handle<std::vector< PileupSummaryInfo > >  PupInfo;
    event.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);    
    std::vector<PileupSummaryInfo>::const_iterator PVI;
    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      if(PVI->getBunchCrossing() == 0) { 
        _npv = vertexs->size();
        _rho  = PVI->getPU_NumInteractions();
        int nTrueInt = PVI->getTrueNumInteractions();
        _npu = nTrueInt;
        _PUReweight = reweight.weight(2012,2012,nTrueInt);
        break;
      } 
    }
  }
  
  _triggerbit = myTriggerHelper->FindTriggerBit(event,foundPaths,indexOfPath);
  _metfilterbit = myTriggerHelper->FindMETBit(event);
  int tbit = _triggerbit;
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
 
  // protect in case of events where trigger hasn't fired --> no collection created 
  event.getByLabel(theCandLabel,candHandle);
  if (!candHandle.isValid()) return;
  event.getByLabel(theCandLabel,candHandle);
  event.getByLabel("jets",jetHandle);
  event.getByLabel("softLeptons",dauHandle);
  event.getByLabel("slimmedMETs",metHandle);
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
  if(theisMC)FillGenInfo(event);

  //Loop of softleptons and fill them
  FillSoftLeptons(daus,event,theFSR);

  //Loop on Jets
  _numberOfJets = 0;
  if(writeJets)_numberOfJets = FillJet(jets);

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
int HTauTauNtuplizer::FillJet(const edm::View<pat::Jet> *jets){
  int nJets=0;
  for(edm::View<pat::Jet>::const_iterator ijet = jets->begin(); ijet!=jets->end();++ijet){
    nJets++;
    _jets_px.push_back( (float) ijet->px());
    _jets_py.push_back( (float) ijet->py());
    _jets_pz.push_back( (float) ijet->pz());
    _jets_e.push_back( (float) ijet->energy());
    _jets_Flavour.push_back(ijet->partonFlavour());
    _jets_PUJetID.push_back(ijet->userFloat("pileupJetId:fullDiscriminant"));
    _bdiscr.push_back(ijet->bDiscriminator("jetBProbabilityBJetTags"));
    _bdiscr2.push_back(ijet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));

    //PF jet ID
    float NHF = ijet->neutralHadronEnergyFraction();
    float NEMF = ijet->neutralEmEnergyFraction();
    float CHF = ijet->chargedHadronEnergyFraction();
    float MUF = ijet->muonEnergyFraction();
    float CEMF = ijet->chargedEmEnergyFraction();
    float NumNeutralParticles =ijet->neutralMultiplicity();
    float NumConst = ijet->chargedMultiplicity()+NumNeutralParticles;
    float CHM = ijet->chargedMultiplicity();
    float absjeta = fabs(ijet->eta());

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
    _jetrawf.push_back(ijet->jecFactor("Uncorrected"));
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
    int byLooseCombinedIsolationDeltaBetaCorr3Hits=-1, byMediumCombinedIsolationDeltaBetaCorr3Hits=-1, byTightCombinedIsolationDeltaBetaCorr3Hits=-1; // tau 13 TeV Iso
    float byCombinedIsolationDeltaBetaCorrRaw3Hits=-1., chargedIsoPtSum=-1., neutralIsoPtSum=-1., puCorrPtSum=-1.; // tau 13 TeV RAW iso info
    int againstMuonLoose3=-1, againstMuonTight3=-1; // tau 13 TeV muon rejection
    int numChargedParticlesSignalCone=-1, numNeutralHadronsSignalCone=-1, numPhotonsSignalCone=-1, numParticlesSignalCone=-1, numChargedParticlesIsoCone=-1, numNeutralHadronsIsoCone=-1, numPhotonsIsoCone=-1, numParticlesIsoCone=-1;
    float leadChargedParticlePt=-1., trackRefPt=-1.;
    int againstElectronVLooseMVA5 =-1, againstElectronLooseMVA5 = -1, againstElectronMediumMVA5 = -1, againstElectronTightMVA5 = -1, againstElectronVTightMVA5 = -1; // tau 13 TeV ele rejection
    int typeOfMuon=0;
    float byIsolationMVA3oldDMwoLTraw=-1, byIsolationMVA3oldDMwLTraw=-1,  byIsolationMVA3newDMwoLTraw=-1,byIsolationMVA3newDMwLTraw=-1;

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
      byLooseCombinedIsolationDeltaBetaCorr3Hits = userdatahelpers::getUserInt (cand, "byLooseCombinedIsolationDeltaBetaCorr3Hits");
      byMediumCombinedIsolationDeltaBetaCorr3Hits = userdatahelpers::getUserInt (cand, "byMediumCombinedIsolationDeltaBetaCorr3Hits");
      byTightCombinedIsolationDeltaBetaCorr3Hits = userdatahelpers::getUserInt (cand, "byTightCombinedIsolationDeltaBetaCorr3Hits");
      byCombinedIsolationDeltaBetaCorrRaw3Hits = userdatahelpers::getUserFloat (cand, "byCombinedIsolationDeltaBetaCorrRaw3Hits");
      byIsolationMVA3oldDMwoLTraw=userdatahelpers::getUserFloat (cand, "byIsolationMVA3oldDMwoLTraw");
      byIsolationMVA3oldDMwLTraw=userdatahelpers::getUserFloat (cand, "byIsolationMVA3oldDMwLTraw");
      byIsolationMVA3newDMwoLTraw=userdatahelpers::getUserFloat (cand, "byIsolationMVA3newDMwoLTraw");
      byIsolationMVA3newDMwLTraw=userdatahelpers::getUserFloat (cand, "byIsolationMVA3newDMwLTraw");
      chargedIsoPtSum = userdatahelpers::getUserFloat (cand, "chargedIsoPtSum");
      neutralIsoPtSum = userdatahelpers::getUserFloat (cand, "neutralIsoPtSum");
      puCorrPtSum = userdatahelpers::getUserFloat (cand, "puCorrPtSum");
      againstMuonLoose3 = userdatahelpers::getUserInt (cand, "againstMuonLoose3");
      againstMuonTight3 = userdatahelpers::getUserInt (cand, "againstMuonTight3");
      againstElectronVLooseMVA5 = userdatahelpers::getUserInt (cand, "againstElectronVLooseMVA5");
      againstElectronLooseMVA5 = userdatahelpers::getUserInt (cand, "againstElectronLooseMVA5");
      againstElectronMediumMVA5 = userdatahelpers::getUserInt (cand, "againstElectronMediumMVA5");
      againstElectronTightMVA5 = userdatahelpers::getUserInt (cand, "againstElectronTightMVA5");
      againstElectronVTightMVA5 = userdatahelpers::getUserInt (cand, "againstElectronVTightMVA5");
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
    _daughters_byLooseCombinedIsolationDeltaBetaCorr3Hits.push_back(byLooseCombinedIsolationDeltaBetaCorr3Hits);
    _daughters_byMediumCombinedIsolationDeltaBetaCorr3Hits.push_back(byMediumCombinedIsolationDeltaBetaCorr3Hits);
    _daughters_byTightCombinedIsolationDeltaBetaCorr3Hits.push_back(byTightCombinedIsolationDeltaBetaCorr3Hits);
    _daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits.push_back(byCombinedIsolationDeltaBetaCorrRaw3Hits);
    _daughters_chargedIsoPtSum.push_back(chargedIsoPtSum);
    _daughters_neutralIsoPtSum.push_back(neutralIsoPtSum);
    _daughters_puCorrPtSum.push_back(puCorrPtSum);
    _daughters_againstMuonLoose3.push_back(againstMuonLoose3);
    _daughters_againstMuonTight3.push_back(againstMuonTight3);
    _daughters_againstElectronVLooseMVA5.push_back(againstElectronVLooseMVA5);
    _daughters_againstElectronLooseMVA5.push_back(againstElectronLooseMVA5);
    _daughters_againstElectronMediumMVA5.push_back(againstElectronMediumMVA5);
    _daughters_againstElectronTightMVA5.push_back(againstElectronTightMVA5);
    _daughters_againstElectronVTightMVA5.push_back(againstElectronVTightMVA5);
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
          int triggerbit = myTriggerHelper->FindTriggerNumber(pathNamesAll[h],true);
          
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

      //cout << j << " - TTT: " << pathName << endl;
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
  isoi=userdatahelpers::getUserFloat(i.daughter(cand1i),"combRelIsoPF");
  isoj=userdatahelpers::getUserFloat(j.daughter(cand1j),"combRelIsoPF");

  if (isoi<isoj)return true;
  else if(isoi>isoj)return false;

  //step 2, leg 1 Pt
  if(i.daughter(cand1i)->pt()>j.daughter(cand1j)->pt()) return true;
  else if(i.daughter(cand1i)->pt()<j.daughter(cand1j)->pt()) return false;

  //step 3, leg 2 ISO
  isoi=userdatahelpers::getUserFloat(i.daughter(1-cand1i),"combRelIsoPF");
  isoj=userdatahelpers::getUserFloat(j.daughter(1-cand1j),"combRelIsoPF");

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
