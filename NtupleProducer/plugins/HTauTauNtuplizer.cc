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
#include <bitset>
//#include <XYZTLorentzVector.h>

// user include files
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"
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
#include <DataFormats/PatCandidates/interface/PackedCandidate.h>
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"

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
#include "Geometry/Records/interface/HcalParametersRcd.h"
#include "FWCore/Framework/interface/eventsetuprecord_registration_macro.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"
#include "CondFormats/GeometryObjects/interface/HcalParameters.h"
#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/HcalCommonData/interface/HcalParametersFromDD.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalTrajectoryExtrapolatorToLine.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"

#include "TLorentzVector.h"

#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/L1Trigger/interface/Jet.h"

 namespace {
//   bool writePhotons = false;  // Write photons in the tree. 
   bool writeJets = true;     // Write jets in the tree. 
   bool writeFatJets = true;
   bool writeSoftLep = false;
   bool writeL1 = true;
   bool DEBUG = false;
   int DETAIL=1;
 }

using namespace std;
using namespace edm;
using namespace reco;

// Map for JEC uncertainty sources
typedef std::map<std::string, std::unique_ptr<JetCorrectionUncertainty>> myJECMap;

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
  //int FillJet(const edm::View<pat::Jet>* jet, const edm::Event&, JetCorrectionUncertainty*);
  int FillJet(const edm::View<pat::Jet>* jet, const edm::Event&, edm::EventSetup const&, JetCorrectionUncertainty*, myJECMap*);
  void FillFatJet(const edm::View<pat::Jet>* fatjets, const edm::Event&);
  void FillSoftLeptons(const edm::View<reco::Candidate> *dauhandler, const edm::Event& event, const edm::EventSetup& setup, bool theFSR, const edm::View<pat::Jet>* jets, const BXVector<l1t::Tau>* l1taus);
  void VBFtrigMatch(const edm::View<pat::Jet>* jet, const edm::Event&); //FRA
  //void FillbQuarks(const edm::Event&);
  void FillGenInfo(const edm::Event&);
  void FillGenJetInfo(const edm::Event&);
  int GetMatchedGen (const reco::Candidate* genL, const edm::Event& event); // return the index of the associated gen particle in the filtered gen collection, in not existing return -1
  void FillL1Obj(const BXVector<l1t::Tau>* taus, const BXVector<l1t::Jet>* jets, const edm::Event& event); // chia
  
  //int CreateFlagsWord (const pat::GenericParticle* part); // build int with each bit containing some boolean flags
  static bool CompareLegs(const reco::Candidate *, const reco::Candidate *);
  float ComputeMT (math::XYZTLorentzVector visP4, float METx, float METy);
  static bool ComparePairsbyPt(pat::CompositeCandidate i, pat::CompositeCandidate j);
  static bool ComparePairsbyIso(pat::CompositeCandidate i, pat::CompositeCandidate j);

  bool refitPV(const edm::Event & iEvent, const edm::EventSetup & iSetup);
  bool findPrimaryVertices(const edm::Event & iEvent, const edm::EventSetup & iSetup);
  TVector3 getPCA(const edm::Event & iEvent, const edm::EventSetup & iSetup,
		  const reco::Track *aTrack, const GlobalPoint & aPoint);

  // ----------member data ---------------------------
  //std::map <int, int> genFlagPosMap_; // to convert from input to output enum format for H/Z decays
  
  //Configs
  TString theFileName;
  bool theFSR;
  Bool_t theisMC;
  Bool_t doCPVariables;
  Bool_t computeQGVar;
  string theJECName;
  // Bool_t theUseNoHFPFMet; // false: PFmet ; true: NoHFPFMet
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
  //edm::EDGetTokenT<vector<l1extra::L1JetParticle>> l1ExtraIsoTau_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<edm::TriggerResults> metFilterBits_;
  edm::EDGetTokenT<vector<Vertex>> theVtxTag;
  edm::EDGetTokenT<edm::View<reco::VertexCompositePtrCandidate>> theSecVtxTag; //FRA
  edm::EDGetTokenT<double> theRhoTag;
  edm::EDGetTokenT<double> theRhoMiniRelIsoTag;
  edm::EDGetTokenT<double> theRhoForJERTag;
  edm::EDGetTokenT<vector<PileupSummaryInfo>> thePUTag;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> thePFCandTag;
  edm::EDGetTokenT<edm::View<pat::CompositeCandidate>> theCandTag;
  edm::EDGetTokenT<edm::View<pat::Jet>> theJetTag;
  edm::EDGetTokenT<edm::View<pat::Jet>> theFatJetTag;
  edm::EDGetTokenT<edm::ValueMap<float>> theQGTaggerTag;
  edm::EDGetTokenT<edm::View<reco::Candidate>> theLepTag;
  edm::EDGetTokenT<LHEEventProduct> theLHETag;
  edm::EDGetTokenT<GenEventInfoProduct> theGenTag;
  edm::EDGetTokenT<pat::METCollection> theMetTag;
  edm::EDGetTokenT<pat::METCollection> theMetERTag;
  edm::EDGetTokenT<pat::METCollection> thePUPPIMetTag;
  edm::EDGetTokenT<math::Error<2>::type> thePFMETCovTag;
  edm::EDGetTokenT<double> thePFMETSignifTag;
  edm::EDGetTokenT<edm::View<pat::GenericParticle>> theGenericTag;
  edm::EDGetTokenT<edm::View<reco::GenJet>> theGenJetTag;
  edm::EDGetTokenT<edm::MergeableCounter> theTotTag;
  edm::EDGetTokenT<edm::MergeableCounter> thePassTag;
  edm::EDGetTokenT<LHEEventProduct> theLHEPTag;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotTag;
  edm::EDGetTokenT<BXVector<l1t::Tau> > theL1TauTag;
  edm::EDGetTokenT<BXVector<l1t::Jet> > theL1JetTag;
  edm::EDGetTokenT<int> theNBadMuTag;
  edm::EDGetTokenT<GenLumiInfoHeader> genLumiHeaderTag;
  edm::EDGetTokenT< bool >ecalBadCalibFilterUpdate_token ;

  //flags
  //static const int nOutVars =14;
  bool applyTrigger;    // Only events passing trigger
  bool applySkim;       //  "     "      "     skim
  //PUReweight reweight;

  //counters
  Int_t Nevt_Gen;
  Int_t Nevt_PassTrigger;
  Int_t Npairs;

  //Event Output variables
  ULong64_t _indexevents;
  Int_t _runNumber;
  Int_t _lumi;
  Long64_t _triggerbit;
  Int_t _metfilterbit;
  Int_t _NBadMu;
  Bool_t _passecalBadCalibFilterUpdate;
  Float_t _met;
  Float_t _metphi;
  Float_t _met_er;
  Float_t _met_er_phi;
  Float_t _PUPPImet;
  Float_t _PUPPImetphi;
  Float_t _PFMETCov00;
  Float_t _PFMETCov01;
  Float_t _PFMETCov10;
  Float_t _PFMETCov11;
  Float_t _PFMETsignif;
  Float_t _MC_weight;
  Float_t _aMCatNLOweight;
  Int_t _npv;
  Float_t _lheHt;
  Int_t   _lheNOutPartons;
  Int_t   _lheNOutB;
  Int_t   _lheNOutC;
  Float_t _npu;
  Int_t   _PUNumInteractions;
  Float_t _PUReweight;
  Float_t _rho;
  Int_t _nup;
  Float_t _MC_weight_scale_muF0p5;
  Float_t _MC_weight_scale_muF2;
  Float_t _MC_weight_scale_muR0p5;
  Float_t _MC_weight_scale_muR2;
  Float_t _pv_x=0, _pv_y=0, _pv_z=0;
  Float_t _pvGen_x=0, _pvGen_y=0, _pvGen_z=0;
  Float_t _pvRefit_x=0, _pvRefit_y=0, _pvRefit_z=0;
  bool _isRefitPV=false;
  
  // pairs
  //std::vector<TLorentzVector> _mothers;
  std::vector<Float_t> _mothers_px;
  std::vector<Float_t> _mothers_py;
  std::vector<Float_t> _mothers_pz;
  std::vector<Float_t> _mothers_e;
  std::vector<Long64_t> _mothers_trgSeparateMatch; // are the two legs matched to different HLT objs?
                                                   // stored bitwise for HLT paths as done for daughters_trgMatch

  
  // reco leptons
  //std::vector<TLorentzVector> _daughters;
  std::vector<string> _trigger_name;
  std::vector<Int_t> _trigger_accept;
  std::vector<Float_t> _daughters_px;
  std::vector<Float_t> _daughters_py;
  std::vector<Float_t> _daughters_pz;
  std::vector<Float_t> _daughters_e;
  
  std::vector<Float_t> _daughters_charged_px;
  std::vector<Float_t> _daughters_charged_py;
  std::vector<Float_t> _daughters_charged_pz;
  std::vector<Float_t> _daughters_charged_e;

  std::vector<Float_t> _daughters_neutral_px;
  std::vector<Float_t> _daughters_neutral_py;
  std::vector<Float_t> _daughters_neutral_pz;
  std::vector<Float_t> _daughters_neutral_e;

  std::vector<Int_t> _daughters_TauUpExists;
  std::vector<Float_t> _daughters_px_TauUp;
  std::vector<Float_t> _daughters_py_TauUp;
  std::vector<Float_t> _daughters_pz_TauUp;
  std::vector<Float_t> _daughters_e_TauUp;
  std::vector<Int_t> _daughters_TauDownExists;
  std::vector<Float_t> _daughters_px_TauDown;
  std::vector<Float_t> _daughters_py_TauDown;
  std::vector<Float_t> _daughters_pz_TauDown;
  std::vector<Float_t> _daughters_e_TauDown;
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

  std::vector<Float_t> _genpart_pca_x;
  std::vector<Float_t> _genpart_pca_y;
  std::vector<Float_t> _genpart_pca_z;
  
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
  std::vector<Int_t> _genpart_TauGenDetailedDecayMode;

  std::vector<Int_t> _genpart_flags; // vector of bit flags bout gen info
  
  // gen jets
  std::vector<Float_t> _genjet_px;
  std::vector<Float_t> _genjet_py;
  std::vector<Float_t> _genjet_pz;
  std::vector<Float_t> _genjet_e;
  std::vector<Int_t> _genjet_partonFlavour; // from matched pat::Jet
  std::vector<Int_t> _genjet_hadronFlavour; // (eh yes, because it is not accessible easily from gen jets)

  //L1 taus
  std::vector<Float_t> _L1_tauEt;
  std::vector<Float_t> _L1_tauEta;
  std::vector<Float_t> _L1_tauPhi;
  std::vector<short int> _L1_tauIso;

  //L1 jets
  std::vector<Float_t> _L1_jetEt;
  std::vector<Float_t> _L1_jetEta;
  std::vector<Float_t> _L1_jetPhi;
  
  //std::vector<math::XYZTLorentzVector> _daughter2;

  //Mothers output variables
  std::vector<Int_t> _indexDau1;
  std::vector<Int_t> _indexDau2;
  std::vector<Float_t> _daughters_HLTpt;
  std::vector<Bool_t>  _daughters_isL1IsoTau28Matched;
  std::vector<Float_t>  _daughters_highestEt_L1IsoTauMatched;
  //std::vector<Int_t> _genDaughters;
  std::vector<Bool_t> _isOSCand;

  std::vector<Float_t> _SVmass;
  std::vector<Float_t> _SVmassTauUp;
  std::vector<Float_t> _SVmassTauDown;
  std::vector<Float_t> _SVmassMETUp;
  std::vector<Float_t> _SVmassMETDown;

  std::vector<Float_t> _SVmassUnc;
  std::vector<Float_t> _SVmassUncTauUp;
  std::vector<Float_t> _SVmassUncTauDown;
  std::vector<Float_t> _SVmassUncMETUp;
  std::vector<Float_t> _SVmassUncMETDown;

  std::vector<Float_t> _SVmassTransverse;
  std::vector<Float_t> _SVmassTransverseTauUp;
  std::vector<Float_t> _SVmassTransverseTauDown;
  std::vector<Float_t> _SVmassTransverseMETUp;
  std::vector<Float_t> _SVmassTransverseMETDown;

  std::vector<Float_t> _SVmassTransverseUnc;
  std::vector<Float_t> _SVmassTransverseUncTauUp;
  std::vector<Float_t> _SVmassTransverseUncTauDown;
  std::vector<Float_t> _SVmassTransverseUncMETUp;
  std::vector<Float_t> _SVmassTransverseUncMETDown;

  std::vector<Float_t> _SVpt;
  std::vector<Float_t> _SVptTauUp;
  std::vector<Float_t> _SVptTauDown;
  std::vector<Float_t> _SVptMETUp;
  std::vector<Float_t> _SVptMETDown;

  std::vector<Float_t> _SVptUnc;
  std::vector<Float_t> _SVptUncTauUp;
  std::vector<Float_t> _SVptUncTauDown;
  std::vector<Float_t> _SVptUncMETUp;
  std::vector<Float_t> _SVptUncMETDown;

  std::vector<Float_t> _SVeta;
  std::vector<Float_t> _SVetaTauUp;
  std::vector<Float_t> _SVetaTauDown;
  std::vector<Float_t> _SVetaMETUp;
  std::vector<Float_t> _SVetaMETDown;

  std::vector<Float_t> _SVetaUnc;
  std::vector<Float_t> _SVetaUncTauUp;
  std::vector<Float_t> _SVetaUncTauDown;
  std::vector<Float_t> _SVetaUncMETUp;
  std::vector<Float_t> _SVetaUncMETDown;

  std::vector<Float_t> _SVphi;
  std::vector<Float_t> _SVphiTauUp;
  std::vector<Float_t> _SVphiTauDown;
  std::vector<Float_t> _SVphiMETUp;
  std::vector<Float_t> _SVphiMETDown;

  std::vector<Float_t> _SVphiUnc;
  std::vector<Float_t> _SVphiUncTauUp;
  std::vector<Float_t> _SVphiUncTauDown;
  std::vector<Float_t> _SVphiUncMETUp;
  std::vector<Float_t> _SVphiUncMETDown;

  std::vector<Float_t> _SVMetRho;
  std::vector<Float_t> _SVMetRhoTauUp;
  std::vector<Float_t> _SVMetRhoTauDown;
  std::vector<Float_t> _SVMetRhoMETUp;
  std::vector<Float_t> _SVMetRhoMETDown;

  std::vector<Float_t> _SVMetPhi;
  std::vector<Float_t> _SVMetPhiTauUp;
  std::vector<Float_t> _SVMetPhiTauDown;
  std::vector<Float_t> _SVMetPhiMETUp;
  std::vector<Float_t> _SVMetPhiMETDown;

  std::vector<Float_t> _metx;
  std::vector<Float_t> _mety;
  std::vector<Float_t> _metx_up;
  std::vector<Float_t> _mety_up;
  std::vector<Float_t> _metx_down;
  std::vector<Float_t> _mety_down;
  std::vector<Float_t> _metx_up_tes;
  std::vector<Float_t> _mety_up_tes;
  std::vector<Float_t> _metx_down_tes;
  std::vector<Float_t> _mety_down_tes;
  std::vector<Float_t> _uncorrmetx;
  std::vector<Float_t> _uncorrmety;
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
  std::vector<Float_t> _combreliso03;
  std::vector<Float_t> _discriminator;//BDT for ele, discriminator for tau, 
  std::vector<Int_t> _daughters_muonID; //bitwise (bit 0 loose, 1 soft , 2 medium, 3 tight, 4 highPT 5 tight_noVtx)
  std::vector<Int_t> _daughters_typeOfMuon; //bitwise, 0=PF, 1=Global, 2=Tracker
  std::vector<Float_t> _dxy;
  std::vector<Float_t> _dz;
  std::vector<Float_t> _dxy_innerTrack;
  std::vector<Float_t> _dz_innerTrack;
  std::vector<Float_t> _daughters_rel_error_trackpt;
  std::vector<Float_t> _SIP;
  std::vector<bool> _daughters_iseleBDT; //isBDT for ele
  std::vector<bool> _daughters_iseleWPLoose; //isBDT for ele
  std::vector<bool> _daughters_iseleWP80; //isBDT for ele
  std::vector<bool> _daughters_iseleWP90; //isBDT for ele
  std::vector<bool> _daughters_iseleNoIsoWPLoose; //isBDT for ele no Iso
  std::vector<bool> _daughters_iseleNoIsoWP80; //isBDT for ele no Iso
  std::vector<bool> _daughters_iseleNoIsoWP90; //isBDT for ele no Iso
  std::vector<Float_t> _daughters_eleMVAnt; //isBDT for ele
  std::vector<Float_t> _daughters_eleMVA_HZZ; //isBDT for ele
  std::vector<bool> _daughters_passConversionVeto; //isBDT for ele
  std::vector<int>  _daughters_eleMissingHits;
  std::vector<int>  _daughters_eleMissingLostHits;
  std::vector<bool>  _daughters_iseleChargeConsistent;
  std::vector<int> _daughters_iseleCUT; //CUT ID for ele (0=veto,1=loose,2=medium,3=tight)
  std::vector<Int_t> _decayType;//for taus only
  std::vector<Long64_t> _daughters_tauID; //bitwise. check h_tauID for histogram list 
  static const int ntauIds = 41;
  TString tauIDStrings[ntauIds] = {
   "byLooseCombinedIsolationDeltaBetaCorr3Hits",
   "byMediumCombinedIsolationDeltaBetaCorr3Hits",
   "byTightCombinedIsolationDeltaBetaCorr3Hits",
   "againstMuonLoose3",
   "againstMuonTight3",
   "againstElectronVLooseMVA6",
   "againstElectronLooseMVA6",
   "againstElectronMediumMVA6",
   "againstElectronTightMVA6",
   "againstElectronVTightMVA6",
   "byVLooseIsolationMVArun2v1DBoldDMwLT",
   "byLooseIsolationMVArun2v1DBoldDMwLT",
   "byMediumIsolationMVArun2v1DBoldDMwLT",
   "byTightIsolationMVArun2v1DBoldDMwLT",
   "byVTightIsolationMVArun2v1DBoldDMwLT",
   "byVLooseIsolationMVArun2v1DBnewDMwLT",
   "byLooseIsolationMVArun2v1DBnewDMwLT",
   "byMediumIsolationMVArun2v1DBnewDMwLT",
   "byTightIsolationMVArun2v1DBnewDMwLT",
   "byVTightIsolationMVArun2v1DBnewDMwLT",
   "byLooseIsolationMVArun2v1DBdR03oldDMwLT",
   "byMediumIsolationMVArun2v1DBdR03oldDMwLT",
   "byTightIsolationMVArun2v1DBdR03oldDMwLT",
   "byVTightIsolationMVArun2v1DBdR03oldDMwLT",
   "byVLooseIsolationMVArun2017v1DBoldDMwLT2017", //FRA syncApr2018
   "byLooseIsolationMVArun2017v1DBoldDMwLT2017",  //FRA syncApr2018
   "byMediumIsolationMVArun2017v1DBoldDMwLT2017", //FRA syncApr2018
   "byTightIsolationMVArun2017v1DBoldDMwLT2017",  //FRA syncApr2018
   "byVTightIsolationMVArun2017v1DBoldDMwLT2017", //FRA syncApr2018
   "byVVLooseIsolationMVArun2017v2DBoldDMwLT2017", //FRA syncApr2018
   "byVLooseIsolationMVArun2017v2DBoldDMwLT2017", //FRA syncApr2018
   "byLooseIsolationMVArun2017v2DBoldDMwLT2017",  //FRA syncApr2018
   "byMediumIsolationMVArun2017v2DBoldDMwLT2017", //FRA syncApr2018
   "byTightIsolationMVArun2017v2DBoldDMwLT2017",  //FRA syncApr2018
   "byVTightIsolationMVArun2017v2DBoldDMwLT2017", //FRA syncApr2018
   "byVVTightIsolationMVArun2017v2DBoldDMwLT2017", //FRA syncApr2018
   "byVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017", //FRA syncApr2018
   "byLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017",  //FRA syncApr2018
   "byMediumIsolationMVArun2017v2DBoldDMdR0p3wLT2017", //FRA syncApr2018
   "byTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017",  //FRA syncApr2018
   "byVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017", //FRA syncApr2018
  };
  std::vector<Float_t> _daughters_IetaIeta;
  std::vector<Float_t> _daughters_full5x5_IetaIeta;
  std::vector<Float_t> _daughters_hOverE;
  std::vector<Float_t> _daughters_deltaEtaSuperClusterTrackAtVtx;
  std::vector<Float_t> _daughters_deltaPhiSuperClusterTrackAtVtx;
  std::vector<Float_t> _daughters_IoEmIoP;
  std::vector<Float_t> _daughters_IoEmIoP_ttH;
  std::vector<Float_t> _daughters_SCeta;
  std::vector<Float_t> _daughters_depositR03_tracker;
  std::vector<Float_t> _daughters_depositR03_ecal;
  std::vector<Float_t> _daughters_depositR03_hcal;
  std::vector<Int_t> _daughters_decayModeFindingOldDMs;
  std::vector<Float_t> _daughters_footprintCorrection;
  std::vector<Float_t> _daughters_neutralIsoPtSumWeight;
  std::vector<Float_t> _daughters_photonPtSumOutsideSignalCone;
  std::vector<Int_t> _daughters_decayModeFindingNewDMs;
  std::vector<Float_t> _daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits;
  std::vector<Float_t> _daughters_byIsolationMVArun2v1DBoldDMwLTraw;
  std::vector<Float_t> _daughters_byIsolationMVArun2017v2DBoldDMwLTraw2017; //FRA
  std::vector<Float_t> _daughters_byIsolationMVArun2017v1DBoldDMwLTraw2017; //FRA
  std::vector<Float_t> _daughters_byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017; //FRA
  std::vector<Int_t> _daughters_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017; //FRA
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
  std::vector<Long64_t> _daughters_trgMatched;
  std::vector<Long64_t> _daughters_FilterFired;
  std::vector<Long64_t> _daughters_isGoodTriggerType;
  std::vector<Long64_t> _daughters_L3FilterFired;
  std::vector<Long64_t> _daughters_L3FilterFiredLast;

  std::vector<Int_t> _daughters_jetNDauChargedMVASel;
  std::vector<Float_t> _daughters_miniRelIsoCharged;
  std::vector<Float_t> _daughters_miniRelIsoNeutral;
  std::vector<Float_t> _daughters_jetPtRel;
  std::vector<Float_t> _daughters_jetPtRatio;
  std::vector<Float_t> _daughters_jetBTagCSV;
  std::vector<Float_t> _daughters_jetBTagDeepCSV;
  std::vector<Float_t> _daughters_lepMVA_mvaId;

  std::vector<Float_t> _daughters_pca_x;
  std::vector<Float_t> _daughters_pca_y;
  std::vector<Float_t> _daughters_pca_z;

  std::vector<Float_t> _daughters_pcaRefitPV_x;
  std::vector<Float_t> _daughters_pcaRefitPV_y;
  std::vector<Float_t> _daughters_pcaRefitPV_z;

  std::vector<Float_t> _daughters_pcaGenPV_x;
  std::vector<Float_t> _daughters_pcaGenPV_y;
  std::vector<Float_t> _daughters_pcaGenPV_z;


  //Jets variables
  Int_t _numberOfJets;
  //std::vector<TLorentzVector> _jets;
  std::vector<std::vector<Long64_t>> _jets_VBFfirstTrigMatch; //FRA
  std::vector<std::vector<Long64_t>> _jets_VBFsecondTrigMatch; //FRA
  std::vector<Long64_t> _jets_VBFleadFilterMatch;    //FRA
  std::vector<Long64_t> _jets_VBFsubleadFilterMatch; //FRA
  std::vector<Float_t> _jets_px;
  std::vector<Float_t> _jets_py;
  std::vector<Float_t> _jets_pz;
  std::vector<Float_t> _jets_e;
  std::vector<Float_t> _jets_rawPt;
  std::vector<Float_t> _jets_area;
  std::vector<Float_t> _jets_mT;
  std::vector<Float_t> _jets_PUJetID;
  std::vector<Float_t> _jets_PUJetIDupdated;
  std::vector<Int_t>   _jets_PUJetIDupdated_WP;
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
  std::vector<Float_t> _jets_MUF;
  std::vector<Int_t>   _jets_neMult;
  std::vector<Int_t>   _jets_chMult;
  std::vector<Float_t> _jets_jecUnc;

  // JEC uncertainty sources
  std::vector<Float_t> _jets_jetUnc_AbsoluteFlavMap_up; // up variations
  std::vector<Float_t> _jets_jetUnc_AbsoluteMPFBias_up;
  std::vector<Float_t> _jets_jetUnc_AbsoluteScale_up;
  std::vector<Float_t> _jets_jetUnc_AbsoluteStat_up;
  std::vector<Float_t> _jets_jetUnc_FlavorQCD_up;
  std::vector<Float_t> _jets_jetUnc_Fragmentation_up;
  std::vector<Float_t> _jets_jetUnc_PileUpDataMC_up;
  std::vector<Float_t> _jets_jetUnc_PileUpPtBB_up;
  std::vector<Float_t> _jets_jetUnc_PileUpPtEC1_up;
  std::vector<Float_t> _jets_jetUnc_PileUpPtEC2_up;
  std::vector<Float_t> _jets_jetUnc_PileUpPtHF_up;
  std::vector<Float_t> _jets_jetUnc_PileUpPtRef_up;
  std::vector<Float_t> _jets_jetUnc_RelativeBal_up;
  std::vector<Float_t> _jets_jetUnc_RelativeFSR_up;
  std::vector<Float_t> _jets_jetUnc_RelativeJEREC1_up;
  std::vector<Float_t> _jets_jetUnc_RelativeJEREC2_up;
  std::vector<Float_t> _jets_jetUnc_RelativeJERHF_up;
  std::vector<Float_t> _jets_jetUnc_RelativePtBB_up;
  std::vector<Float_t> _jets_jetUnc_RelativePtEC1_up;
  std::vector<Float_t> _jets_jetUnc_RelativePtEC2_up;
  std::vector<Float_t> _jets_jetUnc_RelativePtHF_up;
  std::vector<Float_t> _jets_jetUnc_RelativeStatEC_up;
  std::vector<Float_t> _jets_jetUnc_RelativeStatFSR_up;
  std::vector<Float_t> _jets_jetUnc_RelativeStatHF_up;
  std::vector<Float_t> _jets_jetUnc_SinglePionECAL_up;
  std::vector<Float_t> _jets_jetUnc_SinglePionHCAL_up;
  std::vector<Float_t> _jets_jetUnc_TimePtEta_up;
  std::vector<Float_t> _jets_jetUnc_AbsoluteFlavMap_dw; // down variations
  std::vector<Float_t> _jets_jetUnc_AbsoluteMPFBias_dw;
  std::vector<Float_t> _jets_jetUnc_AbsoluteScale_dw;
  std::vector<Float_t> _jets_jetUnc_AbsoluteStat_dw;
  std::vector<Float_t> _jets_jetUnc_FlavorQCD_dw;
  std::vector<Float_t> _jets_jetUnc_Fragmentation_dw;
  std::vector<Float_t> _jets_jetUnc_PileUpDataMC_dw;
  std::vector<Float_t> _jets_jetUnc_PileUpPtBB_dw;
  std::vector<Float_t> _jets_jetUnc_PileUpPtEC1_dw;
  std::vector<Float_t> _jets_jetUnc_PileUpPtEC2_dw;
  std::vector<Float_t> _jets_jetUnc_PileUpPtHF_dw;
  std::vector<Float_t> _jets_jetUnc_PileUpPtRef_dw;
  std::vector<Float_t> _jets_jetUnc_RelativeBal_dw;
  std::vector<Float_t> _jets_jetUnc_RelativeFSR_dw;
  std::vector<Float_t> _jets_jetUnc_RelativeJEREC1_dw;
  std::vector<Float_t> _jets_jetUnc_RelativeJEREC2_dw;
  std::vector<Float_t> _jets_jetUnc_RelativeJERHF_dw;
  std::vector<Float_t> _jets_jetUnc_RelativePtBB_dw;
  std::vector<Float_t> _jets_jetUnc_RelativePtEC1_dw;
  std::vector<Float_t> _jets_jetUnc_RelativePtEC2_dw;
  std::vector<Float_t> _jets_jetUnc_RelativePtHF_dw;
  std::vector<Float_t> _jets_jetUnc_RelativeStatEC_dw;
  std::vector<Float_t> _jets_jetUnc_RelativeStatFSR_dw;
  std::vector<Float_t> _jets_jetUnc_RelativeStatHF_dw;
  std::vector<Float_t> _jets_jetUnc_SinglePionECAL_dw;
  std::vector<Float_t> _jets_jetUnc_SinglePionHCAL_dw;
  std::vector<Float_t> _jets_jetUnc_TimePtEta_dw;
  myJECMap jecSourceUncProviders;
  std::vector<std::string> m_jec_sources = {
    "AbsoluteFlavMap",
    "AbsoluteMPFBias",
    "AbsoluteScale",
    "AbsoluteStat",
    "FlavorQCD",
    "Fragmentation",
    "PileUpDataMC",
    "PileUpPtBB",
    "PileUpPtEC1",
    "PileUpPtEC2",
    "PileUpPtHF",
    "PileUpPtRef",
    "RelativeBal",
    "RelativeFSR",
    "RelativeJEREC1",
    "RelativeJEREC2",
    "RelativeJERHF",
    "RelativePtBB",
    "RelativePtEC1",
    "RelativePtEC2",
    "RelativePtHF",
    "RelativeStatEC",
    "RelativeStatFSR",
    "RelativeStatHF",
    "SinglePionECAL",
    "SinglePionHCAL",
    "TimePtEta" };
  std::map<std::string, std::vector<Float_t>> _SourceUncVal_up;
  std::map<std::string, std::vector<Float_t>> _SourceUncVal_dw;

  std::vector<Float_t> _jets_QGdiscr;

  std::vector<Float_t> _ak8jets_px;
  std::vector<Float_t> _ak8jets_py;
  std::vector<Float_t> _ak8jets_pz;
  std::vector<Float_t> _ak8jets_e;
  std::vector<Float_t> _ak8jets_SoftDropMass;
  std::vector<Float_t> _ak8jets_PrunedMass;
  std::vector<Float_t> _ak8jets_TrimmedMass;
  std::vector<Float_t> _ak8jets_FilteredMass;
  std::vector<Float_t> _ak8jets_tau1; // subjettiness
  std::vector<Float_t> _ak8jets_tau2; // subjettiness
  std::vector<Float_t> _ak8jets_tau3; // subjettiness
  std::vector<Float_t> _ak8jets_tau4; // subjettiness
  std::vector<Float_t> _ak8jets_CSV; // CSV score
  std::vector<Float_t> _ak8jets_deepCSV_probb; // CSV score
  std::vector<Float_t> _ak8jets_deepCSV_probbb; // CSV score
  std::vector<Int_t>   _ak8jets_nsubjets;

  // subjets of ak8 -- store ALL subjets, and link them with an idx to the ak8 jet vectors
  std::vector<Float_t> _subjets_px;
  std::vector<Float_t> _subjets_py;
  std::vector<Float_t> _subjets_pz;
  std::vector<Float_t> _subjets_e;
  std::vector<Float_t> _subjets_CSV;
  std::vector<Float_t> _subjets_deepCSV_probb;
  std::vector<Float_t> _subjets_deepCSV_probbb;
  std::vector<Int_t>   _subjets_ak8MotherIdx;

  std::vector<Int_t> _jets_Flavour; // parton flavour
  std::vector<Int_t> _jets_HadronFlavour; // hadron flavour
  std::vector<Int_t> _jets_genjetIndex; // index of matched gen jet in genjet vector

  std::vector<Float_t> _bdiscr;
  std::vector<Float_t> _bdiscr2; //CSVv2
  std::vector<Float_t> _bdiscr3;
  
  std::vector<Float_t> _bdiscr4; //DeepCSV_probb
  std::vector<Float_t> _bdiscr5; //DeepCSV_probbb
  std::vector<Float_t> _bdiscr6; //DeepCSV_probudsg
  std::vector<Float_t> _bdiscr7; //DeepCSV_probc
  std::vector<Float_t> _bdiscr8; //DeepCSV_probcc
  
  std::vector<Float_t> _bdiscr9;  //DeepFlavor_probb
  std::vector<Float_t> _bdiscr10; //DeepFlavor_probbb
  std::vector<Float_t> _bdiscr11; //DeepFlavor_problepb
  std::vector<Float_t> _bdiscr12; //DeepFlavor_probc
  std::vector<Float_t> _bdiscr13; //DeepFlavor_probuds
  std::vector<Float_t> _bdiscr14; //DeepFlavor_probg
  
  std::vector<Int_t> _jetID; //1=loose, 2=tight, 3=tightlepveto
  std::vector<Float_t> _jetrawf;

  std::vector<Float_t> _jets_JER; // Jet Energy Resolution

  // SUSY info
  TString _susyModel;

  //genH
  //std::vector<Float_t> _genH_px;
  //std::vector<Float_t> _genH_py;
  //std::vector<Float_t> _genH_pz;
  //std::vector<Float_t> _genH_e;

  // not a output tree branch, but used to assess object overlap in trg match
  // use as vTrgMatchedToDau_idx.at(idaughter).at(idxHLTPath)
  std::vector<std::vector<int>> vTrgMatchedToDau_idx;
};

const int HTauTauNtuplizer::ntauIds; // definition of static member

// ----Constructor and Destructor -----
HTauTauNtuplizer::HTauTauNtuplizer(const edm::ParameterSet& pset) : reweight(),
  triggerObjects_      (consumes<pat::TriggerObjectStandAloneCollection> (pset.getParameter<edm::InputTag>("triggerSet"))),
  //l1ExtraIsoTau_       (consumes<vector<l1extra::L1JetParticle>>         (pset.getParameter<edm::InputTag>("l1extraIsoTau"))) ,
  triggerBits_         (consumes<edm::TriggerResults>                    (pset.getParameter<edm::InputTag>("triggerResultsLabel"))),
  metFilterBits_       (consumes<edm::TriggerResults>                    (pset.getParameter<edm::InputTag>("metFilters"))),
  theVtxTag            (consumes<vector<Vertex>>                         (pset.getParameter<edm::InputTag>("vtxCollection"))),
  theSecVtxTag         (consumes<edm::View<reco::VertexCompositePtrCandidate>> (pset.getParameter<edm::InputTag>("secVtxCollection"))), //FRA
  theRhoTag            (consumes<double>                                 (pset.getParameter<edm::InputTag>("rhoCollection"))),
  theRhoMiniRelIsoTag  (consumes<double>                                 (pset.getParameter<edm::InputTag>("rhoMiniRelIsoCollection"))),
  theRhoForJERTag      (consumes<double>                                 (pset.getParameter<edm::InputTag>("rhoForJER"))), //FRA
  thePUTag             (consumes<vector<PileupSummaryInfo>>              (pset.getParameter<edm::InputTag>("puCollection"))),
  thePFCandTag         (consumes<edm::View<pat::PackedCandidate>>        (pset.getParameter<edm::InputTag>("PFCandCollection"))),
  theCandTag           (consumes<edm::View<pat::CompositeCandidate>>     (pset.getParameter<edm::InputTag>("candCollection"))),
  theJetTag            (consumes<edm::View<pat::Jet>>                    (pset.getParameter<edm::InputTag>("jetCollection"))),
  theFatJetTag         (consumes<edm::View<pat::Jet>>                    (pset.getParameter<edm::InputTag>("ak8jetCollection"))),
  theQGTaggerTag       (consumes<edm::ValueMap<float>>                   (pset.getParameter<edm::InputTag>("QGTagger"))),
  theLepTag            (consumes<edm::View<reco::Candidate>>             (pset.getParameter<edm::InputTag>("lepCollection"))),
  theLHETag            (consumes<LHEEventProduct>                        (pset.getParameter<edm::InputTag>("lheCollection"))),
  theGenTag            (consumes<GenEventInfoProduct>                    (pset.getParameter<edm::InputTag>("genCollection"))),
  theMetTag            (consumes<pat::METCollection>                     (pset.getParameter<edm::InputTag>("metCollection"))),
  theMetERTag            (consumes<pat::METCollection>                   (pset.getParameter<edm::InputTag>("metERCollection"))),
  thePUPPIMetTag       (consumes<pat::METCollection>                     (pset.getParameter<edm::InputTag>("PUPPImetCollection"))),
  thePFMETCovTag       (consumes<math::Error<2>::type>                   (pset.getParameter<edm::InputTag>("srcPFMETCov"))),
  thePFMETSignifTag    (consumes<double>                                 (pset.getParameter<edm::InputTag>("srcPFMETSignificance"))),
  theGenericTag        (consumes<edm::View<pat::GenericParticle>>        (pset.getParameter<edm::InputTag>("genericCollection"))),
  theGenJetTag         (consumes<edm::View<reco::GenJet>>                (pset.getParameter<edm::InputTag>("genjetCollection"))),
  theTotTag            (consumes<edm::MergeableCounter, edm::InLumi>     (pset.getParameter<edm::InputTag>("totCollection"))),
  thePassTag           (consumes<edm::MergeableCounter, edm::InLumi>     (pset.getParameter<edm::InputTag>("passCollection"))),
  theLHEPTag           (consumes<LHEEventProduct>                        (pset.getParameter<edm::InputTag>("lhepCollection"))),
  beamSpotTag          (consumes<reco::BeamSpot>                         (pset.getParameter<edm::InputTag>("beamSpot"))),
  theL1TauTag          (consumes<BXVector<l1t::Tau>>                     (pset.getParameter<edm::InputTag>("stage2TauCollection"))),
  theL1JetTag          (consumes<BXVector<l1t::Jet>>                     (pset.getParameter<edm::InputTag>("stage2JetCollection"))),
  theNBadMuTag         (consumes<int>                                    (pset.getParameter<edm::InputTag>("nBadMu"))),
  genLumiHeaderTag     (consumes<GenLumiInfoHeader, edm::InLumi>         (pset.getParameter<edm::InputTag>("genLumiHeaderTag"))),
  ecalBadCalibFilterUpdate_token  (consumes< bool >                      (pset.getParameter<edm::InputTag>("ecalBadCalibReducedMINIAODFilter")))

 {
  theFileName = pset.getUntrackedParameter<string>("fileName");
  theFSR = pset.getParameter<bool>("applyFSR");
  theisMC = pset.getParameter<bool>("IsMC");
  doCPVariables = pset.getParameter<bool>("doCPVariables");
  computeQGVar = pset.getParameter<bool>("computeQGVar");
  theJECName = pset.getUntrackedParameter<string>("JECset");
  // theUseNoHFPFMet = pset.getParameter<bool>("useNOHFMet");
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
  /*
  for (std::vector<edm::ParameterSet>::const_iterator iPSet = HLTList.begin();iPSet != HLTList.end(); ++iPSet) {
    const std::string& hlt = iPSet->getParameter<std::string>("HLT");
    const std::vector<std::string>& path1 = iPSet->getParameter<std::vector<std::string>>("path1");
    const std::vector<std::string>& path2 = iPSet->getParameter<std::vector<std::string>>("path2");
    const int& chan = iPSet->getParameter<int>("channel");
    // Build the mape
    myTriggerHelper->addTriggerMap(hlt,path1,path2,chan);
  }
  */
  if (HLTList.size() > 64) cout << endl << "** HTauTauNtuplizer : Warning : trigger list size exceeds 64, not enough bits in Long64_t type to store all" << endl << endl;
  
  for (std::vector<edm::ParameterSet>::const_iterator iPSet = HLTList.begin();iPSet != HLTList.end(); ++iPSet) {
    const std::string& hlt = iPSet->getParameter<std::string>("HLT");
    const std::vector<std::string>& path1 = iPSet->getParameter<std::vector<std::string>>("path1");
    const std::vector<std::string>& path2 = iPSet->getParameter<std::vector<std::string>>("path2");
    const std::vector<std::string>& path3 = iPSet->getParameter<std::vector<std::string>>("path3"); //FRA
    const std::vector<std::string>& path4 = iPSet->getParameter<std::vector<std::string>>("path4"); //FRA
    const int& leg1 = iPSet->getParameter<int>("leg1");
    const int& leg2 = iPSet->getParameter<int>("leg2");
    const double& pt1 = iPSet->getParameter<double>("pt1"); //FRA
    const double& pt2 = iPSet->getParameter<double>("pt2"); //FRA
    // Build the mape
    //myTriggerHelper->addTriggerMap(hlt,path1,path2,leg1,leg2);
    //myTriggerHelper->addTriggerMap(hlt,path1,path2,path3,path4,leg1,leg2); //FRA
    myTriggerHelper->addTriggerMap(hlt,path1,path2,path3,path4,leg1,leg2, pt1, pt2); //FRA
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
  _susyModel = ""; // not to initialize at every event, as updated ni new lumi block
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
  _mothers_trgSeparateMatch.clear();
  
  //_daughters.clear();
  if(DEBUG){
  _trigger_name.clear();
  _trigger_accept.clear();
  }
  _daughters_px.clear();
  _daughters_py.clear();
  _daughters_pz.clear();
  _daughters_e.clear();

  _daughters_charged_px.clear();
  _daughters_charged_py.clear();
  _daughters_charged_pz.clear();
  _daughters_charged_e.clear();

  _daughters_neutral_px.clear();
  _daughters_neutral_py.clear();
  _daughters_neutral_pz.clear();
  _daughters_neutral_e.clear();
   
  _daughters_TauUpExists.clear();
  _daughters_px_TauUp.clear();
  _daughters_py_TauUp.clear();
  _daughters_pz_TauUp.clear();
  _daughters_e_TauUp.clear();
  _daughters_TauDownExists.clear();
  _daughters_px_TauDown.clear();
  _daughters_py_TauDown.clear();
  _daughters_pz_TauDown.clear();
  _daughters_e_TauDown.clear();
  _daughters_charge.clear();
  _daughters_genindex.clear();
  _daughters_IetaIeta.clear();
  _daughters_full5x5_IetaIeta.clear();
  _daughters_hOverE.clear();
  _daughters_deltaEtaSuperClusterTrackAtVtx.clear();
  _daughters_deltaPhiSuperClusterTrackAtVtx.clear();
  _daughters_IoEmIoP.clear();
  _daughters_IoEmIoP_ttH.clear();
  _daughters_SCeta.clear();
  _daughters_depositR03_tracker.clear();  
  _daughters_depositR03_ecal.clear();  
  _daughters_depositR03_hcal.clear();  
  _daughters_decayModeFindingOldDMs.clear();
  _daughters_decayModeFindingNewDMs.clear();
  _daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits.clear();
  _daughters_byIsolationMVArun2v1DBoldDMwLTraw.clear();
  _daughters_byIsolationMVArun2017v2DBoldDMwLTraw2017.clear();      //FRA
  _daughters_byIsolationMVArun2017v1DBoldDMwLTraw2017.clear();      //FRA
  _daughters_byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017.clear(); //FRA
  _daughters_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017.clear(); //FRA
  _daughters_footprintCorrection.clear();
  _daughters_neutralIsoPtSumWeight.clear();
  _daughters_photonPtSumOutsideSignalCone.clear();
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
  _daughters_trgMatched.clear();
  _daughters_FilterFired.clear();
  _daughters_isGoodTriggerType.clear();
  _daughters_L3FilterFired.clear();
  _daughters_L3FilterFiredLast.clear();

  _daughters_jetNDauChargedMVASel.clear();
  _daughters_miniRelIsoCharged.clear();
  _daughters_miniRelIsoNeutral.clear();
  _daughters_jetPtRel.clear();
  _daughters_jetPtRatio.clear();
  _daughters_jetBTagCSV.clear();
  _daughters_jetBTagDeepCSV.clear();
  _daughters_lepMVA_mvaId.clear();

  _daughters_iseleBDT.clear();
  _daughters_iseleWPLoose.clear();
  _daughters_iseleWP80.clear();
  _daughters_iseleWP90.clear();
  _daughters_iseleNoIsoWPLoose.clear();
  _daughters_iseleNoIsoWP80.clear();
  _daughters_iseleNoIsoWP90.clear();
  _daughters_eleMVAnt.clear();
  _daughters_eleMVA_HZZ.clear();
  _daughters_passConversionVeto.clear();
  _daughters_eleMissingHits.clear();
  _daughters_eleMissingLostHits.clear();
  _daughters_iseleChargeConsistent.clear();
  _daughters_iseleCUT.clear();
  //_daughter2.clear();
  _softLeptons.clear();
  //_genDaughters.clear();

  _daughters_pca_x.clear();
  _daughters_pca_y.clear();
  _daughters_pca_z.clear();

  _daughters_pcaRefitPV_x.clear();
  _daughters_pcaRefitPV_y.clear();
  _daughters_pcaRefitPV_z.clear();

  _daughters_pcaGenPV_x.clear();
  _daughters_pcaGenPV_y.clear();
  _daughters_pcaGenPV_z.clear();
  
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

  _genpart_pca_x.clear();
  _genpart_pca_y.clear();
  _genpart_pca_z.clear();

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
  _genpart_TauGenDetailedDecayMode.clear();
  _genpart_flags.clear();

  _L1_tauEt.clear();
  _L1_tauEta.clear();
  _L1_tauPhi.clear();
  _L1_tauIso.clear();

  _L1_jetEt.clear();
  _L1_jetEta.clear();
  _L1_jetPhi.clear();
  
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
  _SVmassTauUp.clear();
  _SVmassTauDown.clear();
  _SVmassMETUp.clear();
  _SVmassMETDown.clear();

  _SVmassUnc.clear();
  _SVmassUncTauUp.clear();
  _SVmassUncTauDown.clear();
  _SVmassUncMETUp.clear();
  _SVmassUncMETDown.clear();

  _SVmassTransverse.clear();
  _SVmassTransverseTauUp.clear();
  _SVmassTransverseTauDown.clear();
  _SVmassTransverseMETUp.clear();
  _SVmassTransverseMETDown.clear();

  _SVmassTransverseUnc.clear();
  _SVmassTransverseUncTauUp.clear();
  _SVmassTransverseUncTauDown.clear();
  _SVmassTransverseUncMETUp.clear();
  _SVmassTransverseUncMETDown.clear();

  _SVpt.clear();
  _SVptTauUp.clear();
  _SVptTauDown.clear();
  _SVptMETUp.clear();
  _SVptMETDown.clear();

  _SVptUnc.clear();
  _SVptUncTauUp.clear();
  _SVptUncTauDown.clear();
  _SVptUncMETUp.clear();
  _SVptUncMETDown.clear();

  _SVeta.clear();
  _SVetaTauUp.clear();
  _SVetaTauDown.clear();
  _SVetaMETUp.clear();
  _SVetaMETDown.clear();

  _SVetaUnc.clear();
  _SVetaUncTauUp.clear();
  _SVetaUncTauDown.clear();
  _SVetaUncMETUp.clear();
  _SVetaUncMETDown.clear();

  _SVphi.clear();
  _SVphiTauUp.clear();
  _SVphiTauDown.clear();
  _SVphiMETUp.clear();
  _SVphiMETDown.clear();

  _SVphiUnc.clear();
  _SVphiUncTauUp.clear();
  _SVphiUncTauDown.clear();
  _SVphiUncMETUp.clear();
  _SVphiUncMETDown.clear();

  _SVMetRho.clear();
  _SVMetRhoTauUp.clear();
  _SVMetRhoTauDown.clear();
  _SVMetRhoMETUp.clear();
  _SVMetRhoMETDown.clear();

  _SVMetPhi.clear();
  _SVMetPhiTauUp.clear();
  _SVMetPhiTauDown.clear();
  _SVMetPhiMETUp.clear();
  _SVMetPhiMETDown.clear();

  _isOSCand.clear();
  _daughters_HLTpt.clear();
  _daughters_isL1IsoTau28Matched.clear();
  _daughters_highestEt_L1IsoTauMatched.clear();
  _metx.clear();
  _mety.clear();
  _metx_up.clear();
  _mety_up.clear();
  _metx_down.clear();
  _mety_down.clear();
  _metx_up_tes.clear();
  _mety_up_tes.clear();
  _metx_down_tes.clear();
  _mety_down_tes.clear();
  _uncorrmetx.clear();
  _uncorrmety.clear();
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
  _dxy_innerTrack.clear();
  _dz_innerTrack.clear();
  _daughters_rel_error_trackpt.clear();
  _SIP.clear();
  _decayType.clear();
  _daughters_tauID.clear();
  _combreliso.clear();
  _combreliso03.clear();
  _indexevents=0;
  _runNumber=0;
  _lumi=0;
  _NBadMu=0;
  _passecalBadCalibFilterUpdate=false;
  _triggerbit=0;
  _metfilterbit=0;
  _met=0;
  _met_er=0;
  _met_er_phi=0;
  _metphi=0.;
  _PUPPImet=0;
  _PUPPImetphi=0.;
  _PFMETCov00=0.;
  _PFMETCov01=0.;
  _PFMETCov10=0.;
  _PFMETCov11=0.;
  _PFMETsignif=0.;
  _MC_weight=0.;
  _npv=0;
  _lheHt=0;
  _lheNOutPartons=0;
  _lheNOutB=0;
  _lheNOutC=0;
  _npu=0.;
  _PUNumInteractions=0;
  _PUReweight=0.;
  _rho=0;
  _nup=-999;
  _MC_weight_scale_muF0p5=0.;
  _MC_weight_scale_muF2=0.;
  _MC_weight_scale_muR0p5=0.;
  _MC_weight_scale_muR2=0.;

//  _jets.clear();
  _jets_VBFfirstTrigMatch.clear(); //FRA
  _jets_VBFsecondTrigMatch.clear(); //FRA
  _jets_VBFleadFilterMatch.clear();    //FRA
  _jets_VBFsubleadFilterMatch.clear(); //FRA
  _jets_px.clear();
  _jets_py.clear();
  _jets_pz.clear();
  _jets_e.clear();
  _jets_rawPt.clear();
  _jets_area.clear();
  _jets_mT.clear();
  _jets_PUJetID.clear();
  _jets_PUJetIDupdated.clear();
  _jets_PUJetIDupdated_WP.clear();
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
  _jets_MUF.clear();
  _jets_neMult.clear();
  _jets_chMult.clear();
  _jets_Flavour.clear();
  _jets_HadronFlavour.clear();
  _jets_genjetIndex.clear();
  _jets_jecUnc.clear();
  // JEC uncertainty sources
  _jets_jetUnc_AbsoluteFlavMap_up.clear(); //up variations
  _jets_jetUnc_AbsoluteMPFBias_up.clear();
  _jets_jetUnc_AbsoluteScale_up.clear();
  _jets_jetUnc_AbsoluteStat_up.clear();
  _jets_jetUnc_FlavorQCD_up.clear();
  _jets_jetUnc_Fragmentation_up.clear();
  _jets_jetUnc_PileUpDataMC_up.clear();
  _jets_jetUnc_PileUpPtBB_up.clear();
  _jets_jetUnc_PileUpPtEC1_up.clear();
  _jets_jetUnc_PileUpPtEC2_up.clear();
  _jets_jetUnc_PileUpPtHF_up.clear();
  _jets_jetUnc_PileUpPtRef_up.clear();
  _jets_jetUnc_RelativeBal_up.clear();
  _jets_jetUnc_RelativeFSR_up.clear();
  _jets_jetUnc_RelativeJEREC1_up.clear();
  _jets_jetUnc_RelativeJEREC2_up.clear();
  _jets_jetUnc_RelativeJERHF_up.clear();
  _jets_jetUnc_RelativePtBB_up.clear();
  _jets_jetUnc_RelativePtEC1_up.clear();
  _jets_jetUnc_RelativePtEC2_up.clear();
  _jets_jetUnc_RelativePtHF_up.clear();
  _jets_jetUnc_RelativeStatEC_up.clear();
  _jets_jetUnc_RelativeStatFSR_up.clear();
  _jets_jetUnc_RelativeStatHF_up.clear();
  _jets_jetUnc_SinglePionECAL_up.clear();
  _jets_jetUnc_SinglePionHCAL_up.clear();
  _jets_jetUnc_TimePtEta_up.clear();
  _jets_jetUnc_AbsoluteFlavMap_dw.clear(); // down variations
  _jets_jetUnc_AbsoluteMPFBias_dw.clear();
  _jets_jetUnc_AbsoluteScale_dw.clear();
  _jets_jetUnc_AbsoluteStat_dw.clear();
  _jets_jetUnc_FlavorQCD_dw.clear();
  _jets_jetUnc_Fragmentation_dw.clear();
  _jets_jetUnc_PileUpDataMC_dw.clear();
  _jets_jetUnc_PileUpPtBB_dw.clear();
  _jets_jetUnc_PileUpPtEC1_dw.clear();
  _jets_jetUnc_PileUpPtEC2_dw.clear();
  _jets_jetUnc_PileUpPtHF_dw.clear();
  _jets_jetUnc_PileUpPtRef_dw.clear();
  _jets_jetUnc_RelativeBal_dw.clear();
  _jets_jetUnc_RelativeFSR_dw.clear();
  _jets_jetUnc_RelativeJEREC1_dw.clear();
  _jets_jetUnc_RelativeJEREC2_dw.clear();
  _jets_jetUnc_RelativeJERHF_dw.clear();
  _jets_jetUnc_RelativePtBB_dw.clear();
  _jets_jetUnc_RelativePtEC1_dw.clear();
  _jets_jetUnc_RelativePtEC2_dw.clear();
  _jets_jetUnc_RelativePtHF_dw.clear();
  _jets_jetUnc_RelativeStatEC_dw.clear();
  _jets_jetUnc_RelativeStatFSR_dw.clear();
  _jets_jetUnc_RelativeStatHF_dw.clear();
  _jets_jetUnc_SinglePionECAL_dw.clear();
  _jets_jetUnc_SinglePionHCAL_dw.clear();
  _jets_jetUnc_TimePtEta_dw.clear();
  for (std::map<std::string, std::vector<Float_t> >::iterator it=_SourceUncVal_up.begin(); it!=_SourceUncVal_up.end(); ++it)
    {
      it->second.clear();
    }
  for (std::map<std::string, std::vector<Float_t> >::iterator it=_SourceUncVal_dw.begin(); it!=_SourceUncVal_dw.end(); ++it)
    {
      it->second.clear();
    }

  _jets_QGdiscr.clear();
  _numberOfJets=0;
  _bdiscr.clear();
  _bdiscr2.clear();
  _bdiscr3.clear();
  _bdiscr4.clear();
  _bdiscr5.clear();
  _bdiscr6.clear();
  _bdiscr7.clear();
  _bdiscr8.clear();
  _bdiscr9.clear();
  _bdiscr10.clear();
  _bdiscr11.clear();
  _bdiscr12.clear();
  _bdiscr13.clear();
  _bdiscr14.clear();
  
  _jetID.clear();
  _jetrawf.clear();
  _jets_JER.clear(); // Jet Energy Resolution

  _ak8jets_px.clear();
  _ak8jets_py.clear();
  _ak8jets_pz.clear();
  _ak8jets_e.clear();
  _ak8jets_SoftDropMass.clear();
  _ak8jets_PrunedMass.clear();
  _ak8jets_TrimmedMass.clear();
  _ak8jets_FilteredMass.clear();
  _ak8jets_tau1.clear();
  _ak8jets_tau2.clear();
  _ak8jets_tau3.clear();
  _ak8jets_tau4.clear();
  _ak8jets_CSV.clear();
  _ak8jets_deepCSV_probb.clear();
  _ak8jets_deepCSV_probbb.clear();
  _ak8jets_nsubjets.clear();

  _subjets_px.clear();
  _subjets_py.clear();
  _subjets_pz.clear();
  _subjets_e.clear();
  _subjets_CSV.clear();
  _subjets_deepCSV_probb.clear();
  _subjets_deepCSV_probbb.clear();
  _subjets_ak8MotherIdx.clear();



  //_genH_px.clear();
  //_genH_py.clear();
  //_genH_pz.clear();
  //_genH_e.clear();

  // not a tree var, but has to be filled once per daughter - reset here
  vTrgMatchedToDau_idx.clear();
}

void HTauTauNtuplizer::beginJob(){
  edm::Service<TFileService> fs;
  myTree = fs->make<TTree>("HTauTauTree","HTauTauTree");
  int nbins=3+(myTriggerHelper->GetNTriggers());
  hCounter = fs->make<TH1F>("Counters","Counters",nbins,0,nbins);
  hTauIDs = fs->make<TH1F>("TauIDs","TauIDs",ntauIds,0,ntauIds);

  //Branches
  myTree->Branch("EventNumber",&_indexevents,"EventNumber/l");
  myTree->Branch("RunNumber",&_runNumber,"RunNumber/I");
  myTree->Branch("lumi",&_lumi,"lumi/I");
  myTree->Branch("NBadMu",&_NBadMu,"NBadMu/I");
  myTree->Branch("passecalBadCalibFilterUpdate",&_passecalBadCalibFilterUpdate,"passecalBadCalibFilterUpdate/O");
  myTree->Branch("triggerbit",&_triggerbit,"triggerbit/L");
  myTree->Branch("metfilterbit",&_metfilterbit,"metfilterbit/I");
  myTree->Branch("met",&_met,"met/F");
  myTree->Branch("met_er",&_met_er,"met_er/F");
  myTree->Branch("met_er_phi",&_met_er_phi,"met_er_phi/F");
  myTree->Branch("metphi",&_metphi,"metphi/F");
  myTree->Branch("PUPPImet",&_PUPPImet,"PUPPImet/F");
  myTree->Branch("PUPPImetphi",&_PUPPImetphi,"PUPPImetphi/F");
  if(DETAIL>=1){  
    myTree->Branch("daughters_IetaIeta",&_daughters_IetaIeta);
    myTree->Branch("daughters_full5x5_IetaIeta",&_daughters_full5x5_IetaIeta);
    myTree->Branch("daughters_hOverE",&_daughters_hOverE);
    myTree->Branch("daughters_deltaEtaSuperClusterTrackAtVtx",&_daughters_deltaEtaSuperClusterTrackAtVtx);
    myTree->Branch("daughters_deltaPhiSuperClusterTrackAtVtx",&_daughters_deltaPhiSuperClusterTrackAtVtx);
    myTree->Branch("daughters_IoEmIoP",&_daughters_IoEmIoP);
    myTree->Branch("daughters_IoEmIoP_ttH",&_daughters_IoEmIoP_ttH);
  }
  myTree->Branch("PFMETCov00",&_PFMETCov00,"PFMETCov00/F");
  myTree->Branch("PFMETCov01",&_PFMETCov01,"PFMETCov01/F");
  myTree->Branch("PFMETCov10",&_PFMETCov10,"PFMETCov10/F");
  myTree->Branch("PFMETCov11",&_PFMETCov11,"PFMETCov11/F");
  myTree->Branch("PFMETsignif", &_PFMETsignif, "PFMETsignif/F");
  myTree->Branch("npv",&_npv,"npv/I");  
  myTree->Branch("npu",&_npu,"npu/F"); 
  myTree->Branch("PUReweight",&_PUReweight,"PUReweight/F"); 
  myTree->Branch("rho",&_rho,"rho/F");  
  
  myTree->Branch("mothers_px",&_mothers_px);
  myTree->Branch("mothers_py",&_mothers_py);
  myTree->Branch("mothers_pz",&_mothers_pz);
  myTree->Branch("mothers_e",&_mothers_e);
  myTree->Branch("mothers_trgSeparateMatch", &_mothers_trgSeparateMatch);
  if(DEBUG){
    myTree->Branch("trigger_name",&_trigger_name);
    myTree->Branch("trigger_accept",&_trigger_accept);
  }
  myTree->Branch("daughters_px",&_daughters_px);
  myTree->Branch("daughters_py",&_daughters_py);
  myTree->Branch("daughters_pz",&_daughters_pz);
  myTree->Branch("daughters_e",&_daughters_e);
  myTree->Branch("daughters_charge",&_daughters_charge);



  
  if(doCPVariables){
    myTree->Branch("daughters_charged_px",&_daughters_charged_px);
    myTree->Branch("daughters_charged_py",&_daughters_charged_py);
    myTree->Branch("daughters_charged_pz",&_daughters_charged_pz);
    myTree->Branch("daughters_charged_e",&_daughters_charged_e);
    
    myTree->Branch("daughters_neutral_px",&_daughters_neutral_px);
    myTree->Branch("daughters_neutral_py",&_daughters_neutral_py);
    myTree->Branch("daughters_neutral_pz",&_daughters_neutral_pz);
    myTree->Branch("daughters_neutral_e",&_daughters_neutral_e);

  }


  if(writeSoftLep)myTree->Branch("softLeptons",&_softLeptons);
  if(theisMC){
    myTree->Branch("L1_tauEt",&_L1_tauEt);
    myTree->Branch("L1_tauEta",&_L1_tauEta);
    myTree->Branch("L1_tauPhi",&_L1_tauPhi);
    myTree->Branch("L1_tauIso",&_L1_tauIso);
    myTree->Branch("L1_jetEt",&_L1_jetEt);
    myTree->Branch("L1_jetEta",&_L1_jetEta);
    myTree->Branch("L1_jetPhi",&_L1_jetPhi);


     myTree->Branch("daughters_highestEt_L1IsoTauMatched", &_daughters_highestEt_L1IsoTauMatched);
      
    myTree->Branch("daughters_TauUpExists",&_daughters_TauUpExists);
      myTree->Branch("daughters_px_TauUp",&_daughters_px_TauUp);
      myTree->Branch("daughters_py_TauUp",&_daughters_py_TauUp);
      myTree->Branch("daughters_pz_TauUp",&_daughters_pz_TauUp);
      myTree->Branch("daughters_e_TauUp",&_daughters_e_TauUp);

      myTree->Branch("daughters_TauDownExists",&_daughters_TauDownExists);
      myTree->Branch("daughters_px_TauDown",&_daughters_px_TauDown);
      myTree->Branch("daughters_py_TauDown",&_daughters_py_TauDown);
      myTree->Branch("daughters_pz_TauDown",&_daughters_pz_TauDown);
      myTree->Branch("daughters_e_TauDown",&_daughters_e_TauDown);
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
    myTree->Branch("PUNumInteractions",&_PUNumInteractions,"PUNumInteractions/I");  
    myTree->Branch("daughters_genindex",&_daughters_genindex);
    myTree->Branch("MC_weight",&_MC_weight,"MC_weight/F");
    myTree->Branch("MC_weight_scale_muF0p5",&_MC_weight_scale_muF0p5,"MC_weight_scale_muF0p5/F");
    myTree->Branch("MC_weight_scale_muF2",&_MC_weight_scale_muF2,"MC_weight_scale_muF2/F");
    myTree->Branch("MC_weight_scale_muR0p5",&_MC_weight_scale_muR0p5,"MC_weight_scale_muR0p5/F");
    myTree->Branch("MC_weight_scale_muR2",&_MC_weight_scale_muR2,"MC_weight_scale_muR2/F");
    myTree->Branch("lheHt",&_lheHt,"lheHt/F");  
    myTree->Branch("lheNOutPartons", &_lheNOutPartons, "lheNOutPartons/I");
    myTree->Branch("lheNOutB", &_lheNOutB, "lheNOutB/I");
    myTree->Branch("lheNOutC", &_lheNOutC, "lheNOutC/I");
    myTree->Branch("aMCatNLOweight",&_aMCatNLOweight,"aMCatNLOweight/F");    
    myTree->Branch("genpart_px", &_genpart_px);
    myTree->Branch("genpart_py", &_genpart_py);
    myTree->Branch("genpart_pz", &_genpart_pz);
    myTree->Branch("genpart_e", &_genpart_e);
    if(doCPVariables){
      myTree->Branch("genpart_pca_x",&_genpart_pca_x);
      myTree->Branch("genpart_pca_y",&_genpart_pca_y);
      myTree->Branch("genpart_pca_z",&_genpart_pca_z);
    }
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
    myTree->Branch("genpart_TauGenDetailedDecayMode", &_genpart_TauGenDetailedDecayMode);
    myTree->Branch("genpart_flags", &_genpart_flags);

    myTree->Branch("genjet_px", &_genjet_px);
    myTree->Branch("genjet_py", &_genjet_py);
    myTree->Branch("genjet_pz", &_genjet_pz);
    myTree->Branch("genjet_e" , &_genjet_e);
    myTree->Branch("genjet_partonFlavour" , &_genjet_partonFlavour);
    myTree->Branch("genjet_hadronFlavour" , &_genjet_hadronFlavour);

    
    myTree->Branch("NUP", &_nup,"NUP/I");
    myTree->Branch("SVfit_fitMETPhiTauUp", &_SVMetPhiTauUp);
    myTree->Branch("SVfit_fitMETPhiTauDown", &_SVMetPhiTauDown);
    myTree->Branch("SVfit_fitMETPhiMETUp", &_SVMetPhiMETUp);
    myTree->Branch("SVfit_fitMETPhiMETDown", &_SVMetPhiMETDown);
    myTree->Branch("SVfit_fitMETRhoTauUp", &_SVMetRhoTauUp);
    myTree->Branch("SVfit_fitMETRhoTauDown", &_SVMetRhoTauDown);
    myTree->Branch("SVfit_fitMETRhoMETUp", &_SVMetRhoMETUp);
    myTree->Branch("SVfit_fitMETRhoMETDown", &_SVMetRhoMETDown);
    myTree->Branch("SVfit_phiUncTauUp", &_SVphiUncTauUp);
    myTree->Branch("SVfit_phiUncTauDown", &_SVphiUncTauDown);
    myTree->Branch("SVfit_phiUncMETUp", &_SVphiUncMETUp);
    myTree->Branch("SVfit_phiUncMETDown", &_SVphiUncMETDown);
    myTree->Branch("SVfit_phiTauUp", &_SVphiTauUp);
    myTree->Branch("SVfit_phiTauDown", &_SVphiTauDown);
    myTree->Branch("SVfit_phiMETUp", &_SVphiMETUp);
    myTree->Branch("SVfit_phiMETDown", &_SVphiMETDown);
    myTree->Branch("SVfit_etaUncTauUp", &_SVetaUncTauUp);
    myTree->Branch("SVfit_etaUncTauDown", &_SVetaUncTauDown);
    myTree->Branch("SVfit_etaUncMETUp", &_SVetaUncMETUp);
    myTree->Branch("SVfit_etaUncMETDown", &_SVetaUncMETDown);
    myTree->Branch("SVfit_etaTauUp", &_SVetaTauUp);
    myTree->Branch("SVfit_etaTauDown", &_SVetaTauDown);
    myTree->Branch("SVfit_etaMETUp", &_SVetaMETUp);
    myTree->Branch("SVfit_etaMETDown", &_SVetaMETDown);
    myTree->Branch("SVfit_ptUncTauUp", &_SVptUncTauUp);
    myTree->Branch("SVfit_ptUncTauDown", &_SVptUncTauDown);
    myTree->Branch("SVfit_ptUncMETUp", &_SVptUncMETUp);
    myTree->Branch("SVfit_ptUncMETDown", &_SVptUncMETDown);
    myTree->Branch("SVfit_ptTauUp", &_SVptTauUp);
    myTree->Branch("SVfit_ptTauDown", &_SVptTauDown);
    myTree->Branch("SVfit_ptMETUp", &_SVptMETUp);
    myTree->Branch("SVfit_ptMETDown", &_SVptMETDown);
    myTree->Branch("SVfitTransverseMassUncTauUp",&_SVmassTransverseUncTauUp);
    myTree->Branch("SVfitTransverseMassUncTauDown",&_SVmassTransverseUncTauDown);
    myTree->Branch("SVfitTransverseMassUncMETUp",&_SVmassTransverseUncMETUp);
    myTree->Branch("SVfitTransverseMassUncMETDown",&_SVmassTransverseUncMETDown);
    myTree->Branch("SVfitTransverseMassTauUp",&_SVmassTransverseTauUp);
    myTree->Branch("SVfitTransverseMassTauDown",&_SVmassTransverseTauDown);
    myTree->Branch("SVfitTransverseMassMETUp",&_SVmassTransverseMETUp);
    myTree->Branch("SVfitTransverseMassMETDown",&_SVmassTransverseMETDown);
    myTree->Branch("SVfitMassUncTauUp",&_SVmassUncTauUp);
    myTree->Branch("SVfitMassUncTauDown",&_SVmassUncTauDown);
    myTree->Branch("SVfitMassUncMETUp",&_SVmassUncMETUp);
    myTree->Branch("SVfitMassUncMETDown",&_SVmassUncMETDown);
    myTree->Branch("SVfitMassTauUp",&_SVmassTauUp);
    myTree->Branch("SVfitMassTauDown",&_SVmassTauDown);
    myTree->Branch("SVfitMassMETUp",&_SVmassMETUp);
    myTree->Branch("SVfitMassMETDown",&_SVmassMETDown);
  }// end if isMC
  //myTree->Branch("daughters2",&_daughter2);

  myTree->Branch("SVfitMass",&_SVmass);
  myTree->Branch("SVfitMassUnc",&_SVmassUnc);
  myTree->Branch("SVfitTransverseMass",&_SVmassTransverse);
  myTree->Branch("SVfitTransverseMassUnc",&_SVmassTransverseUnc);
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
  myTree->Branch("METx_UP",&_metx_up);
  myTree->Branch("METy_UP",&_mety_up);
  myTree->Branch("METx_DOWN",&_metx_down);
  myTree->Branch("METy_DOWN",&_mety_down);
  myTree->Branch("METx_UP_TES",&_metx_up_tes);
  myTree->Branch("METy_UP_TES",&_mety_up_tes);
  myTree->Branch("METx_DOWN_TES",&_metx_down_tes);
  myTree->Branch("METy_DOWN_TES",&_mety_down_tes);
  myTree->Branch("uncorrMETx",&_uncorrmetx);
  myTree->Branch("uncorrMETy",&_uncorrmety);
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
  myTree->Branch("dxy_innerTrack",&_dxy_innerTrack);
  myTree->Branch("dz_innerTrack",&_dz_innerTrack);
  myTree->Branch("daughters_rel_error_trackpt",&_daughters_rel_error_trackpt);
  myTree->Branch("SIP",&_SIP);
  myTree->Branch("daughters_iseleBDT",&_daughters_iseleBDT);
  myTree->Branch("daughters_iseleWPLoose",&_daughters_iseleWPLoose);
  myTree->Branch("daughters_iseleWP80",&_daughters_iseleWP80);
  myTree->Branch("daughters_iseleWP90",&_daughters_iseleWP90);
  myTree->Branch("daughters_iseleNoIsoWPLoose",&_daughters_iseleNoIsoWPLoose);
  myTree->Branch("daughters_iseleNoIsoWP80",&_daughters_iseleNoIsoWP80);
  myTree->Branch("daughters_iseleNoIsoWP90",&_daughters_iseleNoIsoWP90);
  myTree->Branch("daughters_eleMVAnt",&_daughters_eleMVAnt);
  myTree->Branch("daughters_eleMVA_HZZ",&_daughters_eleMVA_HZZ);
  myTree->Branch("daughters_passConversionVeto",&_daughters_passConversionVeto);
  myTree->Branch("daughters_eleMissingHits",&_daughters_eleMissingHits);
  myTree->Branch("daughters_eleMissingLostHits",&_daughters_eleMissingLostHits);
  myTree->Branch("daughters_iseleChargeConsistent",&_daughters_iseleChargeConsistent);
  myTree->Branch("daughters_eleCUTID",&_daughters_iseleCUT);
  myTree->Branch("decayMode",&_decayType);
  myTree->Branch("tauID",&_daughters_tauID);
  myTree->Branch("combreliso",& _combreliso);
  myTree->Branch("combreliso03",& _combreliso03);
  myTree->Branch("daughters_depositR03_tracker",&_daughters_depositR03_tracker);
  myTree->Branch("daughters_depositR03_ecal",&_daughters_depositR03_ecal);
  myTree->Branch("daughters_depositR03_hcal",&_daughters_depositR03_hcal);
  myTree->Branch("daughters_decayModeFindingOldDMs", &_daughters_decayModeFindingOldDMs);
  myTree->Branch("daughters_SCeta",&_daughters_SCeta);
  myTree->Branch("footprintCorrection",&_daughters_footprintCorrection);
  myTree->Branch("neutralIsoPtSumWeight",&_daughters_neutralIsoPtSumWeight);
  myTree->Branch("photonPtSumOutsideSignalCone",&_daughters_photonPtSumOutsideSignalCone);
  myTree->Branch("daughters_decayModeFindingNewDMs", &_daughters_decayModeFindingNewDMs);
  myTree->Branch("daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits", &_daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits);
  myTree->Branch("daughters_byIsolationMVArun2v1DBoldDMwLTraw",&_daughters_byIsolationMVArun2v1DBoldDMwLTraw);
  myTree->Branch("daughters_byIsolationMVArun2017v2DBoldDMwLTraw2017",&_daughters_byIsolationMVArun2017v2DBoldDMwLTraw2017); //FRA
  myTree->Branch("daughters_byIsolationMVArun2017v1DBoldDMwLTraw2017",&_daughters_byIsolationMVArun2017v1DBoldDMwLTraw2017); //FRA
  myTree->Branch("daughters_byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017",&_daughters_byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017); //FRA
  myTree->Branch("daughters_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017", &_daughters_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017); //FRA
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
  myTree->Branch("daughters_trgMatched", &_daughters_trgMatched);
  myTree->Branch("daughters_isTriggerObjectforPath", &_daughters_L3trigger);
  myTree->Branch("daughters_FilterFired",&_daughters_FilterFired);
  myTree->Branch("daughters_isGoodTriggerType",&_daughters_isGoodTriggerType);
  myTree->Branch("daughters_L3FilterFired",&_daughters_L3FilterFired);
  myTree->Branch("daughters_L3FilterFiredLast",&_daughters_L3FilterFiredLast);
  myTree->Branch("daughters_HLTpt",&_daughters_HLTpt);

  myTree->Branch("daughters_isL1IsoTau28Matched", &_daughters_isL1IsoTau28Matched);

  myTree->Branch("daughters_jetNDauChargedMVASel",&_daughters_jetNDauChargedMVASel);
  myTree->Branch("daughters_miniRelIsoCharged",&_daughters_miniRelIsoCharged);
  myTree->Branch("daughters_miniRelIsoNeutral",&_daughters_miniRelIsoNeutral);
  myTree->Branch("daughters_jetPtRel",&_daughters_jetPtRel);
  myTree->Branch("daughters_jetPtRatio",&_daughters_jetPtRatio);
  myTree->Branch("daughters_jetBTagCSV",&_daughters_jetBTagCSV);
  myTree->Branch("daughters_jetBTagDeepCSV",&_daughters_jetBTagDeepCSV);
  myTree->Branch("daughters_lepMVA_mvaId",&_daughters_lepMVA_mvaId);

  if(doCPVariables){
    myTree->Branch("daughters_pca_x",&_daughters_pca_x);
    myTree->Branch("daughters_pca_y",&_daughters_pca_y);
    myTree->Branch("daughters_pca_z",&_daughters_pca_z);
    myTree->Branch("daughters_pcaRefitPV_x",&_daughters_pcaRefitPV_x);
    myTree->Branch("daughters_pcaRefitPV_y",&_daughters_pcaRefitPV_y);
    myTree->Branch("daughters_pcaRefitPV_z",&_daughters_pcaRefitPV_z);
    myTree->Branch("daughters_pcaGenPV_x",&_daughters_pcaGenPV_x);
    myTree->Branch("daughters_pcaGenPV_y",&_daughters_pcaGenPV_y);
    myTree->Branch("daughters_pcaGenPV_z",&_daughters_pcaGenPV_z);
  }

  myTree->Branch("JetsNumber",&_numberOfJets,"JetsNumber/I");
  myTree->Branch("jets_VBFfirstTrigMatch",&_jets_VBFfirstTrigMatch); //FRA
  myTree->Branch("jets_VBFsecondTrigMatch",&_jets_VBFsecondTrigMatch); //FRA
  myTree->Branch("jets_VBFleadFilterMatch"   , &_jets_VBFleadFilterMatch   ); //FRA
  myTree->Branch("jets_VBFsubleadFilterMatch", &_jets_VBFsubleadFilterMatch); //FRA
  myTree->Branch("jets_px",&_jets_px);
  myTree->Branch("jets_py",&_jets_py);
  myTree->Branch("jets_pz",&_jets_pz);
  myTree->Branch("jets_e",&_jets_e);
  myTree->Branch("jets_rawPt", &_jets_rawPt);
  myTree->Branch("jets_area", &_jets_area);
  myTree->Branch("jets_mT", &_jets_mT);
  myTree->Branch("jets_Flavour",&_jets_Flavour);
  myTree->Branch("jets_HadronFlavour",&_jets_HadronFlavour);
  myTree->Branch("jets_genjetIndex", &_jets_genjetIndex);
  myTree->Branch("jets_PUJetID",&_jets_PUJetID);
  myTree->Branch("jets_PUJetIDupdated",&_jets_PUJetIDupdated);
  myTree->Branch("jets_PUJetIDupdated_WP",&_jets_PUJetIDupdated_WP);
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
  myTree->Branch("jets_MUF"    , &_jets_MUF);
  myTree->Branch("jets_neMult" , &_jets_neMult);
  myTree->Branch("jets_chMult" , &_jets_chMult);
  myTree->Branch("jets_jecUnc" , &_jets_jecUnc);
  // JEC Uncertainty sources
  myTree->Branch("jets_jetUnc_AbsoluteFlavMap_up" , &_SourceUncVal_up["AbsoluteFlavMap"]); // up variations
  myTree->Branch("jets_jetUnc_AbsoluteMPFBias_up"  , &_SourceUncVal_up["AbsoluteMPFBias"]);
  myTree->Branch("jets_jetUnc_AbsoluteScale_up"    , &_SourceUncVal_up["AbsoluteScale"]);
  myTree->Branch("jets_jetUnc_AbsoluteStat_up"     , &_SourceUncVal_up["AbsoluteStat"]);
  myTree->Branch("jets_jetUnc_FlavorQCD_up"        , &_SourceUncVal_up["FlavorQCD"]);
  myTree->Branch("jets_jetUnc_Fragmentation_up"    , &_SourceUncVal_up["Fragmentation"]);
  myTree->Branch("jets_jetUnc_PileUpDataMC_up"     , &_SourceUncVal_up["PileUpDataMC"]);
  myTree->Branch("jets_jetUnc_PileUpPtBB_up"       , &_SourceUncVal_up["PileUpPtBB"]);
  myTree->Branch("jets_jetUnc_PileUpPtEC1_up"      , &_SourceUncVal_up["PileUpPtEC1"]);
  myTree->Branch("jets_jetUnc_PileUpPtEC2_up"      , &_SourceUncVal_up["PileUpPtEC2"]);
  myTree->Branch("jets_jetUnc_PileUpPtHF_up"       , &_SourceUncVal_up["PileUpPtHF"]);
  myTree->Branch("jets_jetUnc_PileUpPtRef_up"      , &_SourceUncVal_up["PileUpPtRef"]);
  myTree->Branch("jets_jetUnc_RelativeBal_up"      , &_SourceUncVal_up["RelativeBal"]);
  myTree->Branch("jets_jetUnc_RelativeFSR_up"      , &_SourceUncVal_up["RelativeFSR"]);
  myTree->Branch("jets_jetUnc_RelativeJEREC1_up"   , &_SourceUncVal_up["RelativeJEREC1"]);
  myTree->Branch("jets_jetUnc_RelativeJEREC2_up"   , &_SourceUncVal_up["RelativeJEREC2"]);
  myTree->Branch("jets_jetUnc_RelativeJERHF_up"    , &_SourceUncVal_up["RelativeJERHF"]);
  myTree->Branch("jets_jetUnc_RelativePtBB_up"     , &_SourceUncVal_up["RelativePtBB"]);
  myTree->Branch("jets_jetUnc_RelativePtEC1_up"    , &_SourceUncVal_up["RelativePtEC1"]);
  myTree->Branch("jets_jetUnc_RelativePtEC2_up"    , &_SourceUncVal_up["RelativePtEC2"]);
  myTree->Branch("jets_jetUnc_RelativePtHF_up"     , &_SourceUncVal_up["RelativePtHF"]);
  myTree->Branch("jets_jetUnc_RelativeStatEC_up"   , &_SourceUncVal_up["RelativeStatEC"]);
  myTree->Branch("jets_jetUnc_RelativeStatFSR_up"  , &_SourceUncVal_up["RelativeStatFSR"]);
  myTree->Branch("jets_jetUnc_RelativeStatHF_up"   , &_SourceUncVal_up["RelativeStatHF"]);
  myTree->Branch("jets_jetUnc_SinglePionECAL_up"   , &_SourceUncVal_up["SinglePionECAL"]);
  myTree->Branch("jets_jetUnc_SinglePionHCAL_up"   , &_SourceUncVal_up["SinglePionHCAL"]);
  myTree->Branch("jets_jetUnc_TimePtEta_up"        , &_SourceUncVal_up["TimePtEta"]);
  myTree->Branch("jets_jetUnc_AbsoluteFlavMap_dw" , &_SourceUncVal_dw["AbsoluteFlavMap"]); // down variations
  myTree->Branch("jets_jetUnc_AbsoluteMPFBias_dw"  , &_SourceUncVal_dw["AbsoluteMPFBias"]);
  myTree->Branch("jets_jetUnc_AbsoluteScale_dw"    , &_SourceUncVal_dw["AbsoluteScale"]);
  myTree->Branch("jets_jetUnc_AbsoluteStat_dw"     , &_SourceUncVal_dw["AbsoluteStat"]);
  myTree->Branch("jets_jetUnc_FlavorQCD_dw"        , &_SourceUncVal_dw["FlavorQCD"]);
  myTree->Branch("jets_jetUnc_Fragmentation_dw"    , &_SourceUncVal_dw["Fragmentation"]);
  myTree->Branch("jets_jetUnc_PileUpDataMC_dw"     , &_SourceUncVal_dw["PileUpDataMC"]);
  myTree->Branch("jets_jetUnc_PileUpPtBB_dw"       , &_SourceUncVal_dw["PileUpPtBB"]);
  myTree->Branch("jets_jetUnc_PileUpPtEC1_dw"      , &_SourceUncVal_dw["PileUpPtEC1"]);
  myTree->Branch("jets_jetUnc_PileUpPtEC2_dw"      , &_SourceUncVal_dw["PileUpPtEC2"]);
  myTree->Branch("jets_jetUnc_PileUpPtHF_dw"       , &_SourceUncVal_dw["PileUpPtHF"]);
  myTree->Branch("jets_jetUnc_PileUpPtRef_dw"      , &_SourceUncVal_dw["PileUpPtRef"]);
  myTree->Branch("jets_jetUnc_RelativeBal_dw"      , &_SourceUncVal_dw["RelativeBal"]);
  myTree->Branch("jets_jetUnc_RelativeFSR_dw"      , &_SourceUncVal_dw["RelativeFSR"]);
  myTree->Branch("jets_jetUnc_RelativeJEREC1_dw"   , &_SourceUncVal_dw["RelativeJEREC1"]);
  myTree->Branch("jets_jetUnc_RelativeJEREC2_dw"   , &_SourceUncVal_dw["RelativeJEREC2"]);
  myTree->Branch("jets_jetUnc_RelativeJERHF_dw"    , &_SourceUncVal_dw["RelativeJERHF"]);
  myTree->Branch("jets_jetUnc_RelativePtBB_dw"     , &_SourceUncVal_dw["RelativePtBB"]);
  myTree->Branch("jets_jetUnc_RelativePtEC1_dw"    , &_SourceUncVal_dw["RelativePtEC1"]);
  myTree->Branch("jets_jetUnc_RelativePtEC2_dw"    , &_SourceUncVal_dw["RelativePtEC2"]);
  myTree->Branch("jets_jetUnc_RelativePtHF_dw"     , &_SourceUncVal_dw["RelativePtHF"]);
  myTree->Branch("jets_jetUnc_RelativeStatEC_dw"   , &_SourceUncVal_dw["RelativeStatEC"]);
  myTree->Branch("jets_jetUnc_RelativeStatFSR_dw"  , &_SourceUncVal_dw["RelativeStatFSR"]);
  myTree->Branch("jets_jetUnc_RelativeStatHF_dw"   , &_SourceUncVal_dw["RelativeStatHF"]);
  myTree->Branch("jets_jetUnc_SinglePionECAL_dw"   , &_SourceUncVal_dw["SinglePionECAL"]);
  myTree->Branch("jets_jetUnc_SinglePionHCAL_dw"   , &_SourceUncVal_dw["SinglePionHCAL"]);
  myTree->Branch("jets_jetUnc_TimePtEta_dw"        , &_SourceUncVal_dw["TimePtEta"]);

  myTree->Branch("bDiscriminator",&_bdiscr);
  myTree->Branch("bCSVscore",&_bdiscr2);
  myTree->Branch("pfCombinedMVAV2BJetTags",&_bdiscr3);
  myTree->Branch("bDeepCSV_probb",&_bdiscr4);
  myTree->Branch("bDeepCSV_probbb",&_bdiscr5);
  myTree->Branch("bDeepCSV_probudsg",&_bdiscr6);
  myTree->Branch("bDeepCSV_probc",&_bdiscr7);
  myTree->Branch("bDeepCSV_probcc",&_bdiscr8);
  myTree->Branch("bDeepFlavor_probb",&_bdiscr9);
  myTree->Branch("bDeepFlavor_probbb",&_bdiscr10);
  myTree->Branch("bDeepFlavor_problepb",&_bdiscr11);
  myTree->Branch("bDeepFlavor_probc",&_bdiscr12);
  myTree->Branch("bDeepFlavor_probuds",&_bdiscr13);
  myTree->Branch("bDeepFlavor_probg",&_bdiscr14);
  
  myTree->Branch("PFjetID",&_jetID);
  myTree->Branch("jetRawf",&_jetrawf);
  myTree->Branch("jets_JER",&_jets_JER);

  myTree->Branch("susyModel",&_susyModel);
  if (computeQGVar) myTree->Branch("jets_QGdiscr" , &_jets_QGdiscr);
  myTree->Branch("ak8jets_px", &_ak8jets_px);
  myTree->Branch("ak8jets_py", &_ak8jets_py);
  myTree->Branch("ak8jets_pz", &_ak8jets_pz);
  myTree->Branch("ak8jets_e", &_ak8jets_e);
  myTree->Branch("ak8jets_SoftDropMass", &_ak8jets_SoftDropMass);
  myTree->Branch("ak8jets_PrunedMass", &_ak8jets_PrunedMass);
  myTree->Branch("ak8jets_TrimmedMass", &_ak8jets_TrimmedMass);
  myTree->Branch("ak8jets_FilteredMass", &_ak8jets_FilteredMass);
  myTree->Branch("ak8jets_tau1", &_ak8jets_tau1);
  myTree->Branch("ak8jets_tau2", &_ak8jets_tau2);
  myTree->Branch("ak8jets_tau3", &_ak8jets_tau3);
  myTree->Branch("ak8jets_tau4", &_ak8jets_tau4);
  myTree->Branch("ak8jets_CSV", &_ak8jets_CSV);
  myTree->Branch("ak8jets_deepCSV_probb", &_ak8jets_deepCSV_probb);
  myTree->Branch("ak8jets_deepCSV_probbb", &_ak8jets_deepCSV_probbb);
  myTree->Branch("ak8jets_nsubjets", &_ak8jets_nsubjets);

  myTree->Branch("subjets_px", &_subjets_px);
  myTree->Branch("subjets_py", &_subjets_py);
  myTree->Branch("subjets_pz", &_subjets_pz);
  myTree->Branch("subjets_e", &_subjets_e);
  myTree->Branch("subjets_CSV", &_subjets_CSV);
  myTree->Branch("subjets_deepCSV_probb", &_subjets_deepCSV_probb);
  myTree->Branch("subjets_deepCSV_probbb", &_subjets_deepCSV_probbb);
  myTree->Branch("subjets_ak8MotherIdx", &_subjets_ak8MotherIdx);


  myTree->Branch("pv_x", &_pv_x);
  myTree->Branch("pv_y", &_pv_y);
  myTree->Branch("pv_z", &_pv_z);

  myTree->Branch("pvRefit_x", &_pvRefit_x);
  myTree->Branch("pvRefit_y", &_pvRefit_y);
  myTree->Branch("pvRefit_z", &_pvRefit_z);

  myTree->Branch("pvGen_x", &_pvGen_x);
  myTree->Branch("pvGen_y", &_pvGen_y);
  myTree->Branch("pvGen_z", &_pvGen_z);

  myTree->Branch("isRefitPV", &_isRefitPV);
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

  //for debugging purposes //FRA
  //if ( event.id().event() != 235 ) return; //FRA
  //cout << "################# Run number: " << event.id().event() << endl; //FRA

  if (doCPVariables) findPrimaryVertices(event, eSetup);
    
  Handle<vector<reco::Vertex> >  vertexs;
  //event.getByLabel("offlineSlimmedPrimaryVertices",vertex);
  event.getByToken(theVtxTag,vertexs);
  
  // JEC fill map with "sourceName, var"
  // Up variations
  _SourceUncVal_up.emplace("AbsoluteFlavMap" ,_jets_jetUnc_AbsoluteFlavMap_up);
  _SourceUncVal_up.emplace("AbsoluteMPFBias" ,_jets_jetUnc_AbsoluteMPFBias_up);
  _SourceUncVal_up.emplace("AbsoluteScale"   ,_jets_jetUnc_AbsoluteScale_up);
  _SourceUncVal_up.emplace("AbsoluteStat"    ,_jets_jetUnc_AbsoluteStat_up);
  _SourceUncVal_up.emplace("FlavorQCD"       ,_jets_jetUnc_FlavorQCD_up);
  _SourceUncVal_up.emplace("Fragmentation"   ,_jets_jetUnc_Fragmentation_up);
  _SourceUncVal_up.emplace("PileUpDataMC"    ,_jets_jetUnc_PileUpDataMC_up);
  _SourceUncVal_up.emplace("PileUpPtBB"      ,_jets_jetUnc_PileUpPtBB_up);
  _SourceUncVal_up.emplace("PileUpPtEC1"     ,_jets_jetUnc_PileUpPtEC1_up);
  _SourceUncVal_up.emplace("PileUpPtEC2"     ,_jets_jetUnc_PileUpPtEC2_up);
  _SourceUncVal_up.emplace("PileUpPtHF"      ,_jets_jetUnc_PileUpPtHF_up);
  _SourceUncVal_up.emplace("PileUpPtRef"     ,_jets_jetUnc_PileUpPtRef_up);
  _SourceUncVal_up.emplace("RelativeBal"     ,_jets_jetUnc_RelativeBal_up);
  _SourceUncVal_up.emplace("RelativeFSR"     ,_jets_jetUnc_RelativeFSR_up);
  _SourceUncVal_up.emplace("RelativeJEREC1"  ,_jets_jetUnc_RelativeJEREC1_up);
  _SourceUncVal_up.emplace("RelativeJEREC2"  ,_jets_jetUnc_RelativeJEREC2_up);
  _SourceUncVal_up.emplace("RelativeJERHF"   ,_jets_jetUnc_RelativeJERHF_up);
  _SourceUncVal_up.emplace("RelativePtBB"    ,_jets_jetUnc_RelativePtBB_up);
  _SourceUncVal_up.emplace("RelativePtEC1"   ,_jets_jetUnc_RelativePtEC1_up);
  _SourceUncVal_up.emplace("RelativePtEC2"   ,_jets_jetUnc_RelativePtEC2_up);
  _SourceUncVal_up.emplace("RelativePtHF"    ,_jets_jetUnc_RelativePtHF_up);
  _SourceUncVal_up.emplace("RelativeStatEC"  ,_jets_jetUnc_RelativeStatEC_up);
  _SourceUncVal_up.emplace("RelativeStatFSR" ,_jets_jetUnc_RelativeStatFSR_up);
  _SourceUncVal_up.emplace("RelativeStatHF"  ,_jets_jetUnc_RelativeStatHF_up);
  _SourceUncVal_up.emplace("SinglePionECAL"  ,_jets_jetUnc_SinglePionECAL_up);
  _SourceUncVal_up.emplace("SinglePionHCAL"  ,_jets_jetUnc_SinglePionHCAL_up);
  _SourceUncVal_up.emplace("TimePtEta"       ,_jets_jetUnc_TimePtEta_up);
  // Down variations
  _SourceUncVal_dw.emplace("AbsoluteFlavMap" ,_jets_jetUnc_AbsoluteFlavMap_dw);
  _SourceUncVal_dw.emplace("AbsoluteMPFBias" ,_jets_jetUnc_AbsoluteMPFBias_dw);
  _SourceUncVal_dw.emplace("AbsoluteScale"   ,_jets_jetUnc_AbsoluteScale_dw);
  _SourceUncVal_dw.emplace("AbsoluteStat"    ,_jets_jetUnc_AbsoluteStat_dw);
  _SourceUncVal_dw.emplace("FlavorQCD"       ,_jets_jetUnc_FlavorQCD_dw);
  _SourceUncVal_dw.emplace("Fragmentation"   ,_jets_jetUnc_Fragmentation_dw);
  _SourceUncVal_dw.emplace("PileUpDataMC"    ,_jets_jetUnc_PileUpDataMC_dw);
  _SourceUncVal_dw.emplace("PileUpPtBB"      ,_jets_jetUnc_PileUpPtBB_dw);
  _SourceUncVal_dw.emplace("PileUpPtEC1"     ,_jets_jetUnc_PileUpPtEC1_dw);
  _SourceUncVal_dw.emplace("PileUpPtEC2"     ,_jets_jetUnc_PileUpPtEC2_dw);
  _SourceUncVal_dw.emplace("PileUpPtHF"      ,_jets_jetUnc_PileUpPtHF_dw);
  _SourceUncVal_dw.emplace("PileUpPtRef"     ,_jets_jetUnc_PileUpPtRef_dw);
  _SourceUncVal_dw.emplace("RelativeBal"     ,_jets_jetUnc_RelativeBal_dw);
  _SourceUncVal_dw.emplace("RelativeFSR"     ,_jets_jetUnc_RelativeFSR_dw);
  _SourceUncVal_dw.emplace("RelativeJEREC1"  ,_jets_jetUnc_RelativeJEREC1_dw);
  _SourceUncVal_dw.emplace("RelativeJEREC2"  ,_jets_jetUnc_RelativeJEREC2_dw);
  _SourceUncVal_dw.emplace("RelativeJERHF"   ,_jets_jetUnc_RelativeJERHF_dw);
  _SourceUncVal_dw.emplace("RelativePtBB"    ,_jets_jetUnc_RelativePtBB_dw);
  _SourceUncVal_dw.emplace("RelativePtEC1"   ,_jets_jetUnc_RelativePtEC1_dw);
  _SourceUncVal_dw.emplace("RelativePtEC2"   ,_jets_jetUnc_RelativePtEC2_dw);
  _SourceUncVal_dw.emplace("RelativePtHF"    ,_jets_jetUnc_RelativePtHF_dw);
  _SourceUncVal_dw.emplace("RelativeStatEC"  ,_jets_jetUnc_RelativeStatEC_dw);
  _SourceUncVal_dw.emplace("RelativeStatFSR" ,_jets_jetUnc_RelativeStatFSR_dw);
  _SourceUncVal_dw.emplace("RelativeStatHF"  ,_jets_jetUnc_RelativeStatHF_dw);
  _SourceUncVal_dw.emplace("SinglePionECAL"  ,_jets_jetUnc_SinglePionECAL_dw);
  _SourceUncVal_dw.emplace("SinglePionHCAL"  ,_jets_jetUnc_SinglePionHCAL_dw);
  _SourceUncVal_dw.emplace("TimePtEta"       ,_jets_jetUnc_TimePtEta_dw);

  //----------------------------------------------------------------------
  // Analyze MC history. THIS HAS TO BE DONE BEFORE ANY RETURN STATEMENT
  // (eg skim or trigger), in order to update the gen counters correctly!!!
  
  //  std::vector<const reco::Candidate *> genZs;
  // std::vector<const reco::Candidate *> genZLeps;
  
  _npv = vertexs->size();
   if (theisMC) {
    Handle<std::vector< PileupSummaryInfo > >  PupInfo;
    //event.getByLabel(edm::InputTag("slimmedAddPileupInfo"), PupInfo);    
    event.getByToken(thePUTag, PupInfo);
    std::vector<PileupSummaryInfo>::const_iterator PVI;
    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      if(PVI->getBunchCrossing() == 0) { 
        _PUNumInteractions  = PVI->getPU_NumInteractions();
        float nTrueInt = PVI->getTrueNumInteractions();
        _npu = nTrueInt;        
        _PUReweight = reweight.weight(2012,2012,nTrueInt);
        break;
      } 
    }
  }


  // pile up information -- rho
  edm::Handle<double> rhoHandle;
  //event.getByLabel("fixedGridRhoFastjetAll", rhoHandle);
  event.getByToken(theRhoTag, rhoHandle);
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
     
    //try {event.getByLabel("externalLHEProducer", lheEventProduct);} catch (...) {;}
    try {event.getByToken(theLHEPTag, lheEventProduct);} catch (...) {;}
    if (lheEventProduct.isValid())
    {
      const lhef::HEPEUP& lheEvent = lheEventProduct->hepeup();
      std::vector<lhef::HEPEUP::FiveVector> lheParticles = lheEvent.PUP;
      double lheHt = 0.;
      int lheNOutPartons = 0;
      int lheNOutB = 0;
      int lheNOutC = 0;
      size_t numParticles = lheParticles.size();
      for ( size_t idxParticle = 0; idxParticle < numParticles; ++idxParticle ) {
        int absPdgId = TMath::Abs(lheEvent.IDUP[idxParticle]);
        int status = lheEvent.ISTUP[idxParticle];
        if ( status == 1 && ((absPdgId >= 1 &&  absPdgId<= 6) ||  absPdgId== 21) ) { // quarks and gluons
            // cout << "DEBUG: APDGID: " << absPdgId << endl;
            lheHt += TMath::Sqrt((lheParticles[idxParticle][0])*(lheParticles[idxParticle][0]) + (lheParticles[idxParticle][1])*(lheParticles[idxParticle][1])); // first entry is px, second py
            ++lheNOutPartons;
            if (absPdgId == 5) ++lheNOutB ;
            if (absPdgId == 4) ++lheNOutC ;
        }
      }
       _lheHt = lheHt;
       _lheNOutPartons = lheNOutPartons;
       _lheNOutB = lheNOutB;
       _lheNOutC = lheNOutC;
     //cout<<"lheHt = "<<lheHt<<endl;
    }
    //else cout << "LHE product not found" << endl;
  }
  
  _triggerbit = myTriggerHelper->FindTriggerBit(event,foundPaths,indexOfPath,triggerBits);
  _metfilterbit = myTriggerHelper->FindMETBit(event, metFilterBits_);
  Long64_t tbit = _triggerbit;
  //std::cout << " -- tbit: " << std::bitset<64>(tbit) << std::endl;
  for(int itr=0;itr<myTriggerHelper->GetNTriggers();itr++) {
    if(myTriggerHelper->IsTriggerFired(tbit,itr)) hCounter->Fill(itr+3);
  }

  //Get candidate collection
  edm::Handle<edm::View<pat::CompositeCandidate>>candHandle;
  edm::Handle<edm::View<reco::Candidate>>dauHandle;
  edm::Handle<edm::View<pat::Jet>>jetHandle;
  edm::Handle<edm::View<pat::Jet>>fatjetHandle;
  edm::Handle<edm::ValueMap<float>>qgTaggerHandle;
  edm::Handle<BXVector<l1t::Tau>>L1TauHandle;
  edm::Handle<BXVector<l1t::Jet>>L1JetHandle;
  edm::Handle<pat::METCollection> metHandle;
  edm::Handle<pat::METCollection> metERHandle;
  edm::Handle<pat::METCollection> PUPPImetHandle;
  edm::Handle<math::Error<2>::type> covHandle;
  edm::Handle<double> METsignficanceHandle;
  edm::Handle<GenFilterInfo> embeddingWeightHandle;
  edm::Handle<edm::TriggerResults> triggerResults;
  edm::Handle<int> NBadMuHandle;
  edm::Handle< bool > passecalBadCalibFilterUpdate ;

  // protect in case of events where trigger hasn't fired --> no collection created 
  event.getByToken(theCandTag,candHandle);
  if (!candHandle.isValid()) return;

  event.getByToken(theCandTag,candHandle);
  event.getByToken(theJetTag,jetHandle);
  if(computeQGVar)
    event.getByToken(theQGTaggerTag,qgTaggerHandle);
  event.getByToken(theL1TauTag,L1TauHandle);
  event.getByToken(theL1JetTag,L1JetHandle);
  event.getByToken(theFatJetTag,fatjetHandle);
  event.getByToken(theLepTag,dauHandle);
  event.getByToken(theMetTag,metHandle);
  event.getByToken(theMetERTag,metERHandle);
  event.getByToken(thePUPPIMetTag,PUPPImetHandle);
  event.getByToken(thePFMETCovTag,covHandle);
  event.getByToken(thePFMETSignifTag,METsignficanceHandle);
  event.getByToken(theNBadMuTag,NBadMuHandle);
  event.getByToken(ecalBadCalibFilterUpdate_token,passecalBadCalibFilterUpdate);

  if(theisMC){
    edm::Handle<LHEEventProduct> lheeventinfo;
    event.getByToken(theLHEPTag,lheeventinfo);

    edm::Handle<GenEventInfoProduct> genEvt;
    event.getByToken(theGenTag,genEvt);
    _aMCatNLOweight=genEvt->weight();
    _MC_weight = _aMCatNLOweight; // duplicated

    if (lheeventinfo.isValid()) {
      _nup=lheeventinfo->hepeup().NUP;
      if (lheeventinfo->weights().size() > 6) // access weights only if weights() is filled
      {
        _MC_weight_scale_muF0p5 = _aMCatNLOweight*(lheeventinfo->weights()[2].wgt)/(lheeventinfo->originalXWGTUP()); // muF = 0.5 | muR = 1
        _MC_weight_scale_muF2 = _aMCatNLOweight*(lheeventinfo->weights()[1].wgt)/(lheeventinfo->originalXWGTUP()); // muF = 2 | muR = 1
        _MC_weight_scale_muR0p5 = _aMCatNLOweight*(lheeventinfo->weights()[6].wgt)/(lheeventinfo->originalXWGTUP()); // muF = 1 | muR = 0.5
        _MC_weight_scale_muR2 = _aMCatNLOweight*(lheeventinfo->weights()[3].wgt)/(lheeventinfo->originalXWGTUP()); // muF = 1 | muR = 2
      }
    }

  }

  /*if (theUseNoHFPFMet) event.getByLabel("slimmedMETsNoHF",metHandle);
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
  }*/


  const edm::View<pat::CompositeCandidate>* cands = candHandle.product();
  const edm::View<reco::Candidate>* daus = dauHandle.product();
  const edm::View<pat::Jet>* jets = jetHandle.product();
  const edm::View<pat::Jet>* fatjets = fatjetHandle.product();
  const pat::MET &met = metHandle->front();
  const pat::MET &met_er = metERHandle->front();
  const pat::MET &PUPPImet = PUPPImetHandle->front();
  const BXVector<l1t::Tau>* L1Tau = L1TauHandle.product();
  const BXVector<l1t::Jet>* L1Jet = L1JetHandle.product();
  //myNtuple->InitializeVariables();
    
  _indexevents = event.id().event();
  _runNumber = event.id().run();
  _lumi=event.luminosityBlock();
  // _met = met.sumEt(); // scalar sum of the pf candidates
  _met = met.pt();
  _met_er = met_er.pt();
  _met_er_phi = met_er.phi();
  _metphi = met.phi();
  _PUPPImet = PUPPImet.pt();
  _PUPPImetphi = PUPPImet.phi();
  _PFMETCov00 = (*covHandle)(0,0); 
  _PFMETCov10 = (*covHandle)(1,0);
  _PFMETCov01 = _PFMETCov10; // (1,0) is the only one saved
  _PFMETCov11 = (*covHandle)(1,1);
  _PFMETsignif = (*METsignficanceHandle);
  _NBadMu = (*NBadMuHandle);
  _passecalBadCalibFilterUpdate =  (*passecalBadCalibFilterUpdate );
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
  FillSoftLeptons(daus,event,eSetup,theFSR,jets,L1Tau);
  
  //Loop on Jets

  // Accessing the JEC uncertainties 
  //ak4  
  edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  eSetup.get<JetCorrectionsRecord>().get("AK4PFchs",JetCorParColl);
  JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  JetCorrectionUncertainty jecUnc (JetCorPar);

  // Accessing the JEC uncertainties sources - !! FIXME !! - seems like uncertainty sources are all the same for MC and DATA (all eras)
  if(theisMC)
  {
    for (const auto& source: m_jec_sources) {
      JetCorrectorParameters source_parameters("JECUncertaintySources/Fall17_17Nov2017_V6_MC_UncertaintySources_AK4PFchs.txt", source);
      std::unique_ptr<JetCorrectionUncertainty> source_uncertainty(new JetCorrectionUncertainty(source_parameters));
      jecSourceUncProviders.emplace(source, std::move(source_uncertainty));
    }
  }
  else
  {
    for (const auto& source: m_jec_sources) {
      JetCorrectorParameters source_parameters("JECUncertaintySources/Fall17_17Nov2017B_V6_DATA_UncertaintySources_AK4PFchs.txt", source);
      std::unique_ptr<JetCorrectionUncertainty> source_uncertainty(new JetCorrectionUncertainty(source_parameters));
      jecSourceUncProviders.emplace(source, std::move(source_uncertainty));
    }
  }

  _numberOfJets = 0;
  if(writeJets){
    //_numberOfJets = FillJet(jets, event, &jecUnc);
    _numberOfJets = FillJet(jets, event, eSetup, &jecUnc, &jecSourceUncProviders);

    if(computeQGVar){ //Needs jetHandle + qgTaggerHandle
      for(auto jet = jetHandle->begin(); jet != jetHandle->end(); ++jet){
        edm::RefToBase<pat::Jet> jetRef(edm::Ref<edm::View<pat::Jet> >(jetHandle, jet - jetHandle->begin()));
        float qgLikelihood = (*qgTaggerHandle)[jetRef];
        _jets_QGdiscr.push_back(qgLikelihood);
      }
    }
    
  }
  if(writeFatJets) FillFatJet(fatjets, event);
  if(writeL1 && theisMC) FillL1Obj(L1Tau, L1Jet, event);
  
  //FRA: Matching jets with trigger filters for VBF
  VBFtrigMatch(jets, event);
    
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
    float thisMETpx_up = cand.userFloat("MEt_px_UP");
    float thisMETpy_up = cand.userFloat("MEt_py_UP");
    float thisMETpx_down = cand.userFloat("MEt_px_DOWN");
    float thisMETpy_down = cand.userFloat("MEt_py_DOWN");
    float thisMETpx_up_tes = cand.userFloat("MEt_px_UP_TES");
    float thisMETpy_up_tes = cand.userFloat("MEt_py_UP_TES");
    float thisMETpx_down_tes = cand.userFloat("MEt_px_DOWN_TES");
    float thisMETpy_down_tes = cand.userFloat("MEt_py_DOWN_TES");
    float thisMETpx_uncorr = ( cand.hasUserFloat("uncorrMEt_px") ) ? cand.userFloat("uncorrMEt_px") : -999.;
    float thisMETpy_uncorr = ( cand.hasUserFloat("uncorrMEt_py") ) ? cand.userFloat("uncorrMEt_py") : -999.;
    
    bool hasUp   = cand.hasUserFloat ("SVfitMassTauUp");
    bool hasDown = cand.hasUserFloat ("SVfitMassTauDown");
    bool hasMETUp   = cand.hasUserFloat ("SVfitMassMETUp");
    bool hasMETDown = cand.hasUserFloat ("SVfitMassMETDown");

    _SVmass.push_back(cand.userFloat("SVfitMass"));
    _SVmassTauUp.push_back  ( (hasUp   ? cand.userFloat("SVfitMassTauUp")   : -999. ));
    _SVmassTauDown.push_back( (hasDown ? cand.userFloat("SVfitMassTauDown") : -999. ));
    _SVmassMETUp.push_back  ( (hasMETUp   ? cand.userFloat("SVfitMassMETUp")   : -999. ));
    _SVmassMETDown.push_back( (hasMETDown ? cand.userFloat("SVfitMassMETDown") : -999. ));

    _SVmassUnc.push_back(cand.userFloat("SVfitMassUnc"));
    _SVmassUncTauUp.push_back  ( (hasUp   ? cand.userFloat("SVfitMassUncTauUp")   : -999. ));
    _SVmassUncTauDown.push_back( (hasDown ? cand.userFloat("SVfitMassUncTauDown") : -999. ));
    _SVmassUncMETUp.push_back  ( (hasMETUp   ? cand.userFloat("SVfitMassUncMETUp")   : -999. ));
    _SVmassUncMETDown.push_back( (hasMETDown ? cand.userFloat("SVfitMassUncMETDown") : -999. ));

    _SVmassTransverse.push_back(cand.userFloat("SVfitTransverseMass"));
    _SVmassTransverseTauUp.push_back  ( (hasUp   ? cand.userFloat("SVfitTransverseMassTauUp")  : -999. ));
    _SVmassTransverseTauDown.push_back( (hasDown ? cand.userFloat("SVfitTransverseMassTauDown"): -999. ));
    _SVmassTransverseMETUp.push_back  ( (hasMETUp   ? cand.userFloat("SVfitTransverseMassMETUp")  : -999. ));
    _SVmassTransverseMETDown.push_back( (hasMETDown ? cand.userFloat("SVfitTransverseMassMETDown"): -999. ));

    _SVmassTransverseUnc.push_back(cand.userFloat("SVfitTransverseMassUnc"));
    _SVmassTransverseUncTauUp.push_back  ( (hasUp   ? cand.userFloat("SVfitTransverseMassUncTauUp")  : -999. ));
    _SVmassTransverseUncTauDown.push_back( (hasDown ? cand.userFloat("SVfitTransverseMassUncTauDown"): -999. ));
    _SVmassTransverseUncMETUp.push_back  ( (hasMETUp   ? cand.userFloat("SVfitTransverseMassUncMETUp")  : -999. ));
    _SVmassTransverseUncMETDown.push_back( (hasMETDown ? cand.userFloat("SVfitTransverseMassUncMETDown"): -999. ));

    _SVpt.push_back(cand.userFloat("SVfit_pt"));
    _SVptTauUp.push_back  ( (hasUp   ? cand.userFloat("SVfit_ptTauUp")  : -999. ));
    _SVptTauDown.push_back( (hasDown ? cand.userFloat("SVfit_ptTauDown"): -999. ));
    _SVptMETUp.push_back  ( (hasMETUp   ? cand.userFloat("SVfit_ptMETUp")  : -999. ));
    _SVptMETDown.push_back( (hasMETDown ? cand.userFloat("SVfit_ptMETDown"): -999. ));

    _SVptUnc.push_back(cand.userFloat("SVfit_ptUnc"));
    _SVptUncTauUp.push_back  ( (hasUp   ? cand.userFloat("SVfit_ptUncTauUp")  : -999. ));
    _SVptUncTauDown.push_back( (hasDown ? cand.userFloat("SVfit_ptUncTauDown"): -999. ));
    _SVptUncMETUp.push_back  ( (hasMETUp   ? cand.userFloat("SVfit_ptUncMETUp")  : -999. ));
    _SVptUncMETDown.push_back( (hasMETDown ? cand.userFloat("SVfit_ptUncMETDown"): -999. ));

    _SVeta.push_back(cand.userFloat("SVfit_eta"));
    _SVetaTauUp.push_back  ( (hasUp   ? cand.userFloat("SVfit_etaTauUp")  : -999. ));
    _SVetaTauDown.push_back( (hasDown ? cand.userFloat("SVfit_etaTauDown"): -999. ));
    _SVetaMETUp.push_back  ( (hasMETUp   ? cand.userFloat("SVfit_etaMETUp")  : -999. ));
    _SVetaMETDown.push_back( (hasMETDown ? cand.userFloat("SVfit_etaMETDown"): -999. ));

    _SVetaUnc.push_back(cand.userFloat("SVfit_etaUnc"));
    _SVetaUncTauUp.push_back  ( (hasUp   ? cand.userFloat("SVfit_etaUncTauUp")  : -999. ));
    _SVetaUncTauDown.push_back( (hasDown ? cand.userFloat("SVfit_etaUncTauDown"): -999. ));
    _SVetaUncMETUp.push_back  ( (hasMETUp   ? cand.userFloat("SVfit_etaUncMETUp")  : -999. ));
    _SVetaUncMETDown.push_back( (hasMETDown ? cand.userFloat("SVfit_etaUncMETDown"): -999. ));

    _SVphi.push_back(cand.userFloat("SVfit_phi"));
    _SVphiTauUp.push_back  ( (hasUp   ? cand.userFloat("SVfit_phiTauUp")  : -999. ));
    _SVphiTauDown.push_back( (hasDown ? cand.userFloat("SVfit_phiTauDown"): -999. ));
    _SVphiMETUp.push_back  ( (hasMETUp   ? cand.userFloat("SVfit_phiMETUp")  : -999. ));
    _SVphiMETDown.push_back( (hasMETDown ? cand.userFloat("SVfit_phiMETDown"): -999. ));

    _SVphiUnc.push_back(cand.userFloat("SVfit_phiUnc"));
    _SVphiUncTauUp.push_back  ( (hasUp   ? cand.userFloat("SVfit_phiUncTauUp")  : -999. ));
    _SVphiUncTauDown.push_back( (hasDown ? cand.userFloat("SVfit_phiUncTauDown"): -999. ));
    _SVphiUncMETUp.push_back  ( (hasMETUp   ? cand.userFloat("SVfit_phiUncMETUp")  : -999. ));
    _SVphiUncMETDown.push_back( (hasMETDown ? cand.userFloat("SVfit_phiUncMETDown"): -999. ));

    _SVMetRho.push_back(cand.userFloat("SVfit_METRho"));
    _SVMetRhoTauUp.push_back  ( (hasUp   ? cand.userFloat("SVfit_METRhoTauUp")  : -999. ));
    _SVMetRhoTauDown.push_back( (hasDown ? cand.userFloat("SVfit_METRhoTauDown"): -999. ));
    _SVMetRhoMETUp.push_back  ( (hasMETUp   ? cand.userFloat("SVfit_METRhoMETUp")  : -999. ));
    _SVMetRhoMETDown.push_back( (hasMETDown ? cand.userFloat("SVfit_METRhoMETDown"): -999. ));

    _SVMetPhi.push_back(cand.userFloat("SVfit_METPhi"));
    _SVMetPhiTauUp.push_back  ( (hasUp   ? cand.userFloat("SVfit_METPhiTauUp")  : -999. ));
    _SVMetPhiTauDown.push_back( (hasDown ? cand.userFloat("SVfit_METPhiTauDown"): -999. ));
    _SVMetPhiMETUp.push_back  ( (hasMETUp   ? cand.userFloat("SVfit_METPhiMETUp")  : -999. ));
    _SVMetPhiMETDown.push_back( (hasMETDown ? cand.userFloat("SVfit_METPhiMETDown"): -999. ));

    _metx.push_back(thisMETpx);
    _mety.push_back(thisMETpy);
    _metx_up.push_back(thisMETpx_up);
    _mety_up.push_back(thisMETpy_up);
    _metx_down.push_back(thisMETpx_down);
    _mety_down.push_back(thisMETpy_down);
    _metx_up_tes.push_back(thisMETpx_up_tes);
    _mety_up_tes.push_back(thisMETpy_up_tes);
    _metx_down_tes.push_back(thisMETpx_down_tes);
    _mety_down_tes.push_back(thisMETpy_down_tes);
    _uncorrmetx.push_back(thisMETpx_uncorr);
    _uncorrmety.push_back(thisMETpy_uncorr);
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

    // use info computed in FillSoftLeptons to check if legs are matched to separate trigger objects for each trigger
    Long64_t trgSeparateMatch = 0;
    vector<int> vMatchesDau1 = vTrgMatchedToDau_idx.at(_indexDau1.back());
    vector<int> vMatchesDau2 = vTrgMatchedToDau_idx.at(_indexDau2.back());

    for (uint trgidx = 0; trgidx < vMatchesDau1.size(); ++trgidx)
    {
      int match1 = vMatchesDau1.at(trgidx) ;
      int match2 = vMatchesDau2.at(trgidx) ;

      if (match1 != match2 || match1 == -1) // if all good (different match) store 1 in the corresp hlt bit
        trgSeparateMatch |= ((Long64_t) 1 << trgidx);
    }
    _mothers_trgSeparateMatch.push_back(trgSeparateMatch);
    // cout << "Mother: " << _mothers_trgSeparateMatch.size() -1 << " " << trgSeparateMatch << endl;



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

//FRA: VBF trigger matching for jets in the VBF topology
void HTauTauNtuplizer::VBFtrigMatch (const edm::View<pat::Jet> *jets, const edm::Event& event){

  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;

  event.getByToken(triggerObjects_, triggerObjects);
  event.getByToken(triggerBits_, triggerBits);
  const edm::TriggerNames &names = event.triggerNames(*triggerBits);

  // Variables for storing the result
  std::vector<Long64_t> VBFfirstTrigMatched;
  std::vector<Long64_t> VBFsecondTrigMatched;
  
  //int numero = 0;
  
  // Loop on the jets in the event
  for(edm::View<pat::Jet>::const_iterator ijet = jets->begin(); ijet!=jets->end();++ijet)
  {
    //numero = numero+1; //FRA
    //std::cout << "Jet " << numero << std::endl; //FRA
    //std::cout << " -- pt:  " << ijet->pt() << std::endl; //FRA
    //std::cout << " -- phi: " << ijet->phi() << std::endl; //FRA
    
    // For each jet, check all the trigger objects, both 'path3' and 'path4'
    Long64_t firstJetMatched  = 0;
    Long64_t secondJetMatched = 0;
  
    // Loop on the Trigger Objects in the event
    for (size_t idxto = 0; idxto < triggerObjects->size(); ++idxto)
    {
      pat::TriggerObjectStandAlone obj = triggerObjects->at(idxto);
      
      // Check if the TO type is a jet, otherwise continue with next TO
      if (obj.type(85) != true ) continue;
    
      // DeltaR2 match between jet and trigger object
      if(deltaR2(obj,*ijet)<0.25)
      {
        // Unpacking Filter Labels
        obj.unpackFilterLabels(event,*triggerBits);
        
        // Unpacking Path Names
        obj.unpackPathNames(names);
        std::vector<std::string> pathNamesAll  = obj.pathNames(false);
        std::vector<std::string> pathNamesLast = obj.pathNames(true);
        
        // debug: checking TO filter labels //FRA
        //if ( obj.type(85))
        //{
          //const std::vector<std::string>& VLabels = obj.filterLabels(); //FRA
          //printing TO labels //FRA
          //std::cout << " -- VLabels for TO "<< idxto << " - pt: " << obj.pt() << " - phi: " << obj.phi() << std::endl; //FRA
          //for (uint ll = 0; ll < VLabels.size(); ++ll) cout << "    -- " << VLabels.at(ll) << endl; //FRA
        //}
        
        // Loop on the HLT path names in the Trigger Object
        for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h)
        {
          int triggerbit = myTriggerHelper->FindTriggerNumber(pathNamesAll[h],true);
          if (triggerbit < 0) continue ; // not a path I want to save
          
          triggerMapper trgmap = myTriggerHelper->GetTriggerMap(pathNamesAll[h]);
          
          /*if (trgmap.GetNfiltersleg3() > 0)
          {
            std::cout << "HLTPath: " << trgmap.GetHLTPath() << std::endl;
            std::cout << "trgmap.GetNfiltersleg3(): " << trgmap.GetNfiltersleg3() << std::endl;
            std::cout << "trgmap.GetNfiltersleg4(): " << trgmap.GetNfiltersleg4() << std::endl;
            for (int ifilt=0;ifilt<trgmap.GetNfiltersleg1();ifilt++) std::cout << " - Leg 1 : " << trgmap.Getfilter(true,ifilt) << std::endl;
            for (int ifilt=0;ifilt<trgmap.GetNfiltersleg2();ifilt++) std::cout << " - Leg 2 : " << trgmap.Getfilter(false,ifilt) << std::endl;
            for (int ifilt=0;ifilt<trgmap.GetNfiltersleg3();ifilt++) std::cout << " - Leg 3 : " << trgmap.GetfilterVBF(true,ifilt) << std::endl;
            for (int ifilt=0;ifilt<trgmap.GetNfiltersleg4();ifilt++) std::cout << " - Leg 4 : " << trgmap.GetfilterVBF(false,ifilt) << std::endl;
          }*/
          
          bool VBFfirstMatch = true;
          bool VBFsecondMatch = true;
          
          // TO filter labels
          const std::vector<std::string>& vLabels = obj.filterLabels();
        
          // Loop on Leg3 filter (2 jets with pt40)
          if (trgmap.GetNfiltersleg3()>0)
          {
            for(int ifilt=0;ifilt<trgmap.GetNfiltersleg3();ifilt++) //change to leg3
            {
              string label = trgmap.GetfilterVBF(true,ifilt);
              //std::cout << " ---------------------------------------------- @@ leg 3 looking for " << label << std::endl;
              if (label.empty()) {VBFfirstMatch=false; continue;}
              if (find(vLabels.begin(), vLabels.end(), label) == vLabels.end()) VBFfirstMatch=false;
            }
          }
          else VBFfirstMatch = false;
          //std::cout << "VBFfirstMatch: " << VBFfirstMatch << std::endl; //FRA
          
          // Loop on Leg4 filter (1 jet with pt115)
          if (trgmap.GetNfiltersleg4()>0)
          {
            for(int ifilt=0;ifilt<trgmap.GetNfiltersleg4();ifilt++) //change to leg4
            {
              string label = trgmap.GetfilterVBF(false,ifilt);
              //std::cout << " ---------------------------------------------- @@ leg 4 looking for " << label << std::endl;
              if (label.empty()) {VBFsecondMatch=false; continue;}
              if (find(vLabels.begin(), vLabels.end(), label) == vLabels.end()) VBFsecondMatch=false;
            }
          }
          else VBFsecondMatch = false;
          //std::cout << "VBFsecondMatch: " << VBFsecondMatch << std::endl; //FRA

          
          // If two jets matching the VBF filters are found, the events has
          // the event has fired the VBF HLT path
          //std::cout << "first: " << VBFfirstMatch << " - second: " << VBFsecondMatch << std::endl;
          if(VBFfirstMatch)
          {
            //std::cout << " ###### GOOD FIRST VBF MATCH ######" << std::endl; //FRA
            //std::cout << " ***** searching trigger : " << myTriggerHelper -> printTriggerName(triggerbit) << " " << trgmap.GetHLTPath() << std::endl; //FRA
            //std::cout << " --> triggerbit: " << triggerbit << std::endl; //FRA
            //std::cout << "TO filters: " << std::endl;
            //for (unsigned int kj=0;kj<vLabels.size();kj++) std::cout << " - filter : " << vLabels[kj] << std::endl;
            firstJetMatched |= (long(1) <<triggerbit);
          }
          
          if(VBFsecondMatch)
          {
            //std::cout << " ###### GOOD SECOND VBF MATCH ######" << std::endl; //FRA
            //std::cout << " ***** searching trigger : " << myTriggerHelper -> printTriggerName(triggerbit) << " " << trgmap.GetHLTPath() << std::endl; //FRA
            //std::cout << " --> triggerbit: " << triggerbit << std::endl; //FRA
            //std::cout << "TO filters: " << std::endl;
            //for (unsigned int kj=0;kj<vLabels.size();kj++) std::cout << " - filter : " << vLabels[kj] << std::endl;
            secondJetMatched |= (long(1) <<triggerbit);
          }
        } // Loop on HLT paths
        
      } // DeltaR2 match
      
    } // Loop on Trigger Objects
    
    VBFfirstTrigMatched.push_back(firstJetMatched);
    VBFsecondTrigMatched.push_back(secondJetMatched);
    
  } // Loop on jets
  
  // Fill the tree variable
  _jets_VBFfirstTrigMatch.push_back(VBFfirstTrigMatched);
  _jets_VBFsecondTrigMatch.push_back(VBFsecondTrigMatched);

}


//Fill jets quantities
//int HTauTauNtuplizer::FillJet(const edm::View<pat::Jet> *jets, const edm::Event& event, JetCorrectionUncertainty* jecUnc){
int HTauTauNtuplizer::FillJet(const edm::View<pat::Jet> *jets, const edm::Event& event, edm::EventSetup const& iSetup, JetCorrectionUncertainty* jecUnc, myJECMap* jecSourceUncProviders){

  // TriggerBits and TriggerObjets (for VBF trigger matching)
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;

  event.getByToken(triggerObjects_, triggerObjects);
  event.getByToken(triggerBits_, triggerBits);
  const edm::TriggerNames &names = event.triggerNames(*triggerBits);

  // Save a subsample of TOs with (type==85) && (hasLeadFilter || hasSubLeadFilter)
  std::vector<pat::TriggerObjectStandAlone> goodTOs;
  for (size_t idxto = 0; idxto < triggerObjects->size(); ++idxto)
  {
    // Get the TO and unpack the filter labels
    pat::TriggerObjectStandAlone TO = triggerObjects->at(idxto);
    TO.unpackFilterLabels(event,*triggerBits);

    if ( (TO.type(85)==true) && (TO.hasFilterLabel("hltMatchedVBFOnePFJet2CrossCleanedFromDoubleLooseChargedIsoPFTau20") || TO.hasFilterLabel("hltMatchedVBFTwoPFJets2CrossCleanedFromDoubleLooseChargedIsoPFTau20") ) )
      goodTOs.push_back(TO);
  }

  // Getting the primary Vertex (FRA 2017)
  Handle<vector<reco::Vertex> >  vertexs;
  event.getByToken(theVtxTag,vertexs);
  const auto & pv = (*vertexs)[0]; // get the firstPV
  
  // Getting the secondaryVertexCollection (FRA 2017)
  Handle<edm::View<reco::VertexCompositePtrCandidate>> secVtxsHandle;
  event.getByToken(theSecVtxTag,secVtxsHandle);
  if (!secVtxsHandle.isValid())
  {
    throw cms::Exception("ProductNotValid") << "slimmedSecondaryVertices product not valid";
  }
  const edm::View<reco::VertexCompositePtrCandidate> * secVtxs = secVtxsHandle.product();

  // Getting the Rho for Jet Enegy Resolution
  edm::Handle<double> rhoJERHandle;
  event.getByToken(theRhoForJERTag, rhoJERHandle);

  // Initialize Jet Energy Resolution (FRA 2017)
  JME::JetResolution JERresolution = JME::JetResolution::get(iSetup, "AK4PFchs_pt");

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
    _jets_PUJetIDupdated.push_back(ijet->hasUserFloat("pileupJetIdUpdated:fullDiscriminant") ? ijet->userFloat("pileupJetIdUpdated:fullDiscriminant") : -999);
    _jets_PUJetIDupdated_WP.push_back(ijet->hasUserInt("pileupJetIdUpdated:fullId") ? ijet->userInt("pileupJetIdUpdated:fullId") : -999);

    
    //float vtxPx = ijet->userFloat ("vtxPx");                      //FRA: not anymore available in 2017
    //float vtxPy = ijet->userFloat ("vtxPy");                      //FRA: not anymore available in 2017
    //_jets_vtxMass.push_back(ijet->userFloat("vtxMass"));          //FRA: not anymore available in 2017
    //_jets_vtx3dL. push_back(ijet->userFloat("vtx3DVal"));         //FRA: not anymore available in 2017
    //_jets_vtxNtrk.push_back(ijet->userFloat("vtxNtracks"));       //FRA: not anymore available in 2017
    //_jets_vtx3deL.push_back(ijet->userFloat("vtx3DSig"));         //FRA: not anymore available in 2017

    //FRA: new way to access these variables (in 2017)
    Float_t vtxPt   = 0.0;
    Float_t vtxMass = 0.0;
    Float_t vtx3dL  = 0.0;
    Float_t vtxNtrk = 0.0;
    Float_t vtx3deL = 0.0;
    
    VertexDistance3D vdist;
    float maxFoundSignificance=0;
    for(edm::View<reco::VertexCompositePtrCandidate>::const_iterator isv = secVtxs->begin(); isv!=secVtxs->end(); ++isv)
    {
      GlobalVector flightDir(isv->vertex().x() - pv.x(), isv->vertex().y() - pv.y(), isv->vertex().z() - pv.z());
	  GlobalVector jetDir(ijet->px(),ijet->py(),ijet->pz());
      if( deltaR2( flightDir, jetDir ) < 0.09 )
      {
        Measurement1D dl= vdist.distance(pv,VertexState(RecoVertex::convertPos(isv->position()),RecoVertex::convertError(isv->error())));
		if(dl.significance() > maxFoundSignificance)
        {
		  maxFoundSignificance=dl.significance();
		  vtxPt   = isv->pt();
		  vtxMass = isv->p4().M();
		  vtx3dL  = dl.value();
		  vtx3deL = dl.error();
		  vtxNtrk = isv->numberOfSourceCandidatePtrs();
		}
      }
    }
    
    _jets_vtxPt.  push_back(vtxPt);
    _jets_vtxMass.push_back(vtxMass);
    _jets_vtx3dL. push_back(vtx3dL);
    _jets_vtxNtrk.push_back(vtxNtrk);
    _jets_vtx3deL.push_back(vtx3deL);
    

    _bdiscr.push_back(ijet->bDiscriminator("pfJetProbabilityBJetTags"));
    _bdiscr2.push_back(ijet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    _bdiscr3.push_back(ijet->bDiscriminator("pfCombinedMVAV2BJetTags"));
    
    // DeepCSV
    _bdiscr4.push_back(ijet->bDiscriminator("pfDeepCSVJetTags:probb"));
    _bdiscr5.push_back(ijet->bDiscriminator("pfDeepCSVJetTags:probbb"));
    _bdiscr6.push_back(ijet->bDiscriminator("pfDeepCSVJetTags:probudsg"));
    _bdiscr7.push_back(ijet->bDiscriminator("pfDeepCSVJetTags:probc"));
    _bdiscr8.push_back(ijet->bDiscriminator("pfDeepCSVJetTags:probcc"));
    
    // DeepFlavor
    _bdiscr9.push_back(ijet->bDiscriminator("pfDeepFlavourJetTags:probb"));
    _bdiscr10.push_back(ijet->bDiscriminator("pfDeepFlavourJetTags:probbb"));
    _bdiscr11.push_back(ijet->bDiscriminator("pfDeepFlavourJetTags:problepb"));
    _bdiscr12.push_back(ijet->bDiscriminator("pfDeepFlavourJetTags:probc"));
    _bdiscr13.push_back(ijet->bDiscriminator("pfDeepFlavourJetTags:probuds"));
    _bdiscr14.push_back(ijet->bDiscriminator("pfDeepFlavourJetTags:probg"));

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
    _jets_neMult .push_back(NumNeutralParticles);  
    _jets_MUF    .push_back(MUF);

    int jetid=0; 
    //PHYS14
    /*
    if((NHF<0.99 && NEMF<0.99 && NumConst>1 && MUF<0.8) && ((absjeta<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || absjeta>2.4)){
      jetid++;
      if( (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((absjeta<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || absjeta>2.4)  ) jetid++;
    }
    */
    //Spring15
    // if(absjeta<=3.0){
    //   if((NHF<0.99 && NEMF<0.99 && NumConst>1) && ((absjeta<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || absjeta>2.4) ){
    //     jetid++;
    //     if( (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((absjeta<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || absjeta>2.4) ) {
    //       jetid++;
    //       if( (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((absjeta<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || absjeta>2.4)) jetid++;
    //     }
    //   }
    // }else{
    //   if(NEMF<0.90 && NumNeutralParticles>10 ){
    //     jetid++;
    //     jetid++; //TIGHT and LOOSE are the same in this eta region
    //   }
    // }  
    
    // https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Recommendations_for_13_TeV_data
    // bool looseJetID = false;
    // bool tightJetID = false;
    // bool tightLepVetoJetID = false;
    // if (absjeta <= 2.7)
    // {
    //   looseJetID = ( (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((absjeta<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || absjeta>2.4) );
    //   tightJetID = ( (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((absjeta<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || absjeta>2.4) );
    //   tightLepVetoJetID = ( (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((absjeta<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || absjeta>2.4) );
    // }
    // else if (absjeta <= 3.0)
    // {
    //   looseJetID = (NEMF<0.90 && NumNeutralParticles>2 ) ;
    //   tightJetID = looseJetID;
    // }
    // else
    // {
    //   looseJetID = (NEMF<0.90 && NumNeutralParticles>10 );
    //   tightJetID = looseJetID;
    // }
    // if (looseJetID) ++jetid;
    // if (tightJetID) ++jetid;
    // if (tightLepVetoJetID) ++jetid;

    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID13TeVRun2017
    // Since the tight JetID efficiency is > 99% everywhere, loose is not recommended anymore
    bool tightJetID = false;
    bool tightLepVetoJetID = false;
    if (absjeta <= 2.7)
    {
      tightJetID = ( (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((absjeta<=2.4 && CHF>0 && CHM>0) || absjeta>2.4) );
      tightLepVetoJetID = ( (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((absjeta<=2.4 && CHF>0 && CHM>0 && CEMF<0.80) || absjeta>2.4) );
    }
    else if (absjeta <= 3.0)
    {
      tightJetID = (NEMF>0.02 && NEMF<0.99 && NumNeutralParticles>2 );
    }
    else
    {
      tightJetID = (NEMF<0.90 && NHF>0.02 && NumNeutralParticles>10 );
    }
    if (tightJetID) ++jetid;
    if (tightLepVetoJetID) ++jetid;

    _jetID.push_back(jetid);
    float jecFactor = ijet->jecFactor("Uncorrected") ;
    float jetRawPt = jecFactor * ijet->pt();
    //float jetRawPt2 = ijet->pt() / jecFactor; // this is wrong
    _jets_rawPt.push_back ( jetRawPt );
    _jets_area.push_back (ijet->jetArea());
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

    // JEC uncertainties sources
    for (myJECMap::iterator it=jecSourceUncProviders->begin(); it!=jecSourceUncProviders->end(); ++it)
    {
      // up variations
      it->second->setJetEta(ijet->eta());
      it->second->setJetPt(ijet->pt());
      float uncertainty_up = it->second->getUncertainty(true);
      _SourceUncVal_up[it->first].push_back(uncertainty_up);

      // down variations
      it->second->setJetEta(ijet->eta());
      it->second->setJetPt(ijet->pt());
      float uncertainty_dw = it->second->getUncertainty(false);
      _SourceUncVal_dw[it->first].push_back(uncertainty_dw);
    }

    // Jet Energy Resolution
    JME::JetParameters JER_parameters;
    JER_parameters.setJetPt(ijet->pt());
    JER_parameters.setJetEta(ijet->eta());
    JER_parameters.setRho(*rhoJERHandle);
    _jets_JER.push_back( JERresolution.getResolution(JER_parameters) * ijet->energy() ); // JER*energy beacuse cmssw gives the % of JER, while KinFit wants resolution in GeV


    // VBF trigger matching
    Long64_t VBFleadFilterMatch    = 0;
    Long64_t VBFsubleadFilterMatch = 0;
    // loop on the TOs that are already of type==85 and passed one of the VBF jet filters
    for (size_t idxto = 0; idxto < goodTOs.size(); ++idxto)
    {
      // Get the TO
      pat::TriggerObjectStandAlone obj = goodTOs.at(idxto);

      // bools to save if it has the filter or not
      bool isVBFsubleadFilterMatched = false;
      bool isVBFleadFilterMatched = false;

      // DeltaR2 match between jet and trigger object
      if(deltaR2(obj,*ijet)<0.25)
      {
        // Unpacking Filter Labels and Path Names
        obj.unpackFilterLabels(event,*triggerBits);
        obj.unpackPathNames(names);

        //Get HLT path names
        std::vector<std::string> pathNamesAll  = obj.pathNames(false);

        // Loop on the HLT path names in the Trigger Object
        for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h)
        {
          int triggerbit = myTriggerHelper->FindTriggerNumber(pathNamesAll[h],true);
          if (triggerbit < 0) continue ; // not a path I want to save

          triggerMapper trgmap = myTriggerHelper->GetTriggerMap(pathNamesAll[h]);

          // TO filter labels
          const std::vector<std::string>& vLabels = obj.filterLabels();

          // Loop on Leg3 filter (2 jets with pt40) ---> subleadFilter
          if (trgmap.GetNfiltersleg3()>0)
          {
            for(int ifilt=0;ifilt<trgmap.GetNfiltersleg3();ifilt++)
            {
              string label = trgmap.GetfilterVBF(true,ifilt);
              if (label.empty()) continue;
              if (find(vLabels.begin(), vLabels.end(), label) != vLabels.end()) isVBFsubleadFilterMatched=true;
            }
          }
          else isVBFsubleadFilterMatched = false;

          // Loop on Leg4 filter (1 jet with pt115)  ---> leadFilter
          if (trgmap.GetNfiltersleg4()>0)
          {
            for(int ifilt=0;ifilt<trgmap.GetNfiltersleg4();ifilt++) //change to leg4
            {
              string label = trgmap.GetfilterVBF(false,ifilt);
              if (label.empty()) continue;
              if (find(vLabels.begin(), vLabels.end(), label) != vLabels.end()) isVBFleadFilterMatched=true;
            }
          }
          else isVBFleadFilterMatched = false;

        } // end loop on HLT paths

      } // end if dR2<0.25

      if(isVBFsubleadFilterMatched) VBFsubleadFilterMatch |= (Long64_t(1) <<idxto);
      if(isVBFleadFilterMatched)    VBFleadFilterMatch    |= (Long64_t(1) <<idxto);

    } //end loop on goodTOs

    // Fill branches with result of VBF trig matching
    _jets_VBFleadFilterMatch.push_back(VBFleadFilterMatch);
    _jets_VBFsubleadFilterMatch.push_back(VBFsubleadFilterMatch);

  } // end loop on jets

  if ( theisMC )
  {
    edm::Handle<edm::View<reco::GenJet>> genJetHandle;
    //event.getByLabel ("slimmedGenJets", genJetHandle);
    event.getByToken (theGenJetTag, genJetHandle);
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

void HTauTauNtuplizer::FillFatJet(const edm::View<pat::Jet>* fatjets, const edm::Event&)
{
    for(edm::View<pat::Jet>::const_iterator ijet = fatjets->begin(); ijet!=fatjets->end();++ijet)
    {
      _ak8jets_px.push_back( (float) ijet->px());
      _ak8jets_py.push_back( (float) ijet->py());
      _ak8jets_pz.push_back( (float) ijet->pz());
      _ak8jets_e.push_back( (float) ijet->energy());    
      //_ak8jets_SoftDropMass.push_back (ijet->hasUserFloat("ak8PFJetsCHSSoftDropMass") ? ijet->userFloat("ak8PFJetsCHSSoftDropMass") : -999 );
      //_ak8jets_PrunedMass.push_back   (ijet->hasUserFloat("ak8PFJetsCHSPrunedMass")   ? ijet->userFloat("ak8PFJetsCHSPrunedMass")   : -999 );
      //_ak8jets_TrimmedMass.push_back  (ijet->hasUserFloat("ak8PFJetsCHSTrimmedMass")  ? ijet->userFloat("ak8PFJetsCHSTrimmedMass")  : -999 );
      //_ak8jets_FilteredMass.push_back (ijet->hasUserFloat("ak8PFJetsCHSFilteredMass") ? ijet->userFloat("ak8PFJetsCHSFilteredMass") : -999 );
      //_ak8jets_tau1.push_back         (ijet->hasUserFloat("NjettinessAK8:tau1")       ? ijet->userFloat("NjettinessAK8:tau1")       : -999 );
      //_ak8jets_tau2.push_back         (ijet->hasUserFloat("NjettinessAK8:tau2")       ? ijet->userFloat("NjettinessAK8:tau2")       : -999 );
      //_ak8jets_tau3.push_back         (ijet->hasUserFloat("NjettinessAK8:tau3")       ? ijet->userFloat("NjettinessAK8:tau3")       : -999 );
      _ak8jets_SoftDropMass.push_back (ijet->hasUserFloat("ak8PFJetsPuppiSoftDropMass") ? ijet->userFloat("ak8PFJetsPuppiSoftDropMass") : -999 );
      _ak8jets_PrunedMass.push_back   (ijet->hasUserFloat("ak8PFJetsPuppiPrunedMass")   ? ijet->userFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSPrunedMass")   : -999 );
      _ak8jets_TrimmedMass.push_back  (ijet->hasUserFloat("ak8PFJetsPuppiTrimmedMass")  ? ijet->userFloat("ak8PFJetsPuppiTrimmedMass")  : -999 );
      _ak8jets_FilteredMass.push_back (ijet->hasUserFloat("ak8PFJetsPuppiFilteredMass") ? ijet->userFloat("ak8PFJetsPuppiFilteredMass") : -999 );
      _ak8jets_tau1.push_back         (ijet->hasUserFloat("NjettinessAK8Puppi:tau1")    ? ijet->userFloat("NjettinessAK8Puppi:tau1")    : -999 );
      _ak8jets_tau2.push_back         (ijet->hasUserFloat("NjettinessAK8Puppi:tau2")    ? ijet->userFloat("NjettinessAK8Puppi:tau2")    : -999 );
      _ak8jets_tau3.push_back         (ijet->hasUserFloat("NjettinessAK8Puppi:tau3")    ? ijet->userFloat("NjettinessAK8Puppi:tau3")    : -999 );
      _ak8jets_tau4.push_back         (ijet->hasUserFloat("NjettinessAK8Puppi:tau4")    ? ijet->userFloat("NjettinessAK8Puppi:tau4")    : -999 );
      _ak8jets_CSV.push_back(ijet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
      _ak8jets_deepCSV_probb.push_back(ijet->bDiscriminator("pfDeepCSVJetTags:probb"));
      _ak8jets_deepCSV_probbb.push_back(ijet->bDiscriminator("pfDeepCSVJetTags:probbb"));

      // store subjets for soft drop
      int nsubj = 0;
      //if (ijet->hasSubjets("SoftDrop"))
      if (ijet->hasSubjets("SoftDropPuppi"))
      {
        //pat::JetPtrCollection const & subj = ijet->subjets("SoftDrop");
        pat::JetPtrCollection const & subj = ijet->subjets("SoftDropPuppi");
        // cout << "============= IJET " << ijet - fatjets->begin() << " ============= " << endl;
        for (auto isubj = subj.begin(); isubj!=subj.end(); isubj++)
        {
          nsubj += 1;
          _subjets_px.push_back((*isubj)->px());
          _subjets_py.push_back((*isubj)->py());
          _subjets_pz.push_back((*isubj)->pz());
          _subjets_e.push_back((*isubj)->energy());
          _subjets_CSV.push_back((*isubj)->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
          _subjets_deepCSV_probb.push_back((*isubj)->bDiscriminator("pfDeepCSVJetTags:probb"));
          _subjets_deepCSV_probbb.push_back((*isubj)->bDiscriminator("pfDeepCSVJetTags:probbb"));
          _subjets_ak8MotherIdx.push_back(ijet - fatjets->begin()); // idx of fatjet in fatjet vector
          // cout << " * " << (*isubj)->pt() << " " << isubj - subj.begin() << endl;
        }
      }
      _ak8jets_nsubjets.push_back(nsubj);
    }
}


//Fill all leptons (we keep them all for veto purposes
void HTauTauNtuplizer::FillSoftLeptons(const edm::View<reco::Candidate> *daus,
				       const edm::Event& event, const edm::EventSetup& setup,
				       bool theFSR, const edm::View<pat::Jet> *jets,const BXVector<l1t::Tau>* l1taus ){
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;

  event.getByToken(triggerObjects_, triggerObjects);
  event.getByToken(triggerBits_, triggerBits);
  const edm::TriggerNames &names = event.triggerNames(*triggerBits);
  
  //edm::Handle<vector<l1extra::L1JetParticle>> L1ExtraIsoTau;
  //event.getByToken(l1ExtraIsoTau_, L1ExtraIsoTau);


  edm::Handle<edm::View<pat::PackedCandidate> >pfCandHandle;
  event.getByToken(thePFCandTag,pfCandHandle);
  const edm::View<pat::PackedCandidate>* pfCands = pfCandHandle.product();

  std::vector<const pat::PackedCandidate *> pfCands_charged;
  std::vector<const pat::PackedCandidate *> pfCands_neutral;
  LeptonIsoHelper::PFIso_particles(pfCands, pfCands_charged, pfCands_neutral);

  edm::Handle<double> rhoHandle_miniRelIso;
  event.getByToken(theRhoMiniRelIsoTag, rhoHandle_miniRelIso);
  float rho_miniRelIso = *rhoHandle_miniRelIso;

  for(edm::View<reco::Candidate>::const_iterator daui = daus->begin(); daui!=daus->end();++daui){

    const reco::Candidate* cand = &(*daui);
    math::XYZTLorentzVector pfour = cand->p4();
    math::XYZTLorentzVector chargedP4 = cand->p4();
    math::XYZTLorentzVector neutralP4;
    
    TLorentzVector pfourTauUp;
    TLorentzVector pfourTauDown;

    bool existUp   = userdatahelpers::hasUserInt(cand,"TauUpExists"); // simply to check if the userfloat exists, i.e. if it is a tau
    bool existDown = userdatahelpers::hasUserInt(cand,"TauDownExists"); // simply to check if the userfloat exists, i.e. if it is a tau

    int hasUp   = ( existUp   ? userdatahelpers::getUserInt(cand,"TauUpExists")   : false) ;   // actual check of the value of the userfloat
    int hasDown = ( existDown ? userdatahelpers::getUserInt(cand,"TauDownExists") : false) ; // actual check of the value of the userfloat
   
    
    if(hasUp)
    {
      pfourTauUp.SetPxPyPzE(userdatahelpers::getUserFloat(cand,"px_TauUp"),userdatahelpers::getUserFloat(cand,"py_TauUp"),userdatahelpers::getUserFloat(cand,"pz_TauUp"),userdatahelpers::getUserFloat(cand,"e_TauUp"));
    }
    else
    {
      pfourTauUp.SetPxPyPzE(-999.,-999.,-999.,-999.);
    }
    if(hasDown)
    {
      pfourTauDown.SetPxPyPzE(userdatahelpers::getUserFloat(cand,"px_TauDown"),userdatahelpers::getUserFloat(cand,"py_TauDown"),userdatahelpers::getUserFloat(cand,"pz_TauDown"),userdatahelpers::getUserFloat(cand,"e_TauDown"));
    }
    else
    {
      pfourTauDown.SetPxPyPzE(-999.,-999.,-999.,-999.);
    }
    
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

    _daughters_TauUpExists.push_back( (existUp ? hasUp : 0) );
    _daughters_px_TauUp.push_back((float)pfourTauUp.Px());
    _daughters_py_TauUp.push_back((float)pfourTauUp.Py());
    _daughters_pz_TauUp.push_back((float)pfourTauUp.Pz());
    _daughters_e_TauUp.push_back((float)pfourTauUp.E());

    _daughters_TauDownExists.push_back( (existDown ? hasDown : 0) );
    _daughters_px_TauDown.push_back((float)pfourTauDown.Px());
    _daughters_py_TauDown.push_back((float)pfourTauDown.Py());
    _daughters_pz_TauDown.push_back((float)pfourTauDown.Pz());
    _daughters_e_TauDown.push_back((float)pfourTauDown.E());

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
    _combreliso03.push_back( userdatahelpers::hasUserFloat(cand,"combRelIsoPF03") ? userdatahelpers::getUserFloat(cand,"combRelIsoPF03") : -1 );
    _dxy.push_back(userdatahelpers::getUserFloat(cand,"dxy"));
    _dz.push_back(userdatahelpers::getUserFloat(cand,"dz"));
    //_SIP.push_back(userdatahelpers::getUserFloat(cand,"SIP"));
    //int type = -1; 
    //if( userdatahelpers::hasUserInt(cand,"isTESShifted") ) type = ParticleType::TAU;
    //else if (userdatahelpers::hasUserInt(cand,"isPFMuon")) type = ParticleType::MUON;
    //else if (userdatahelpers::hasUserInt(cand,"isEleID80")) type = ParticleType::ELECTRON;
    //else printf("===ERROR!!! UNRECOGNIZED PARTICLE\n");
     int type = ParticleType::TAU;
     if(cand->isMuon()) type = ParticleType::MUON;
     else if(cand->isElectron()) type = ParticleType::ELECTRON;
    _particleType.push_back(type);
    

    //Find closest jet for lepton MVA
    float dRmin_cand_jet = 0.4;
    pat::Jet closest_jet;

    for(edm::View<pat::Jet>::const_iterator jeti = jets->begin(); jeti!=jets->end();++jeti){

      float dR_cand_jet = deltaR(*cand,*jeti);
      if(dR_cand_jet<dRmin_cand_jet){
	closest_jet = (*jeti);
	dRmin_cand_jet = dR_cand_jet;
      }

    }


    // variables
    float discr=-1.;
    int muIDflag = 0;
    bool isgood = false;
    bool iseleLoose = false;
    bool isele80=false;
    bool isele90=false;
    bool iselenoisoLoose = false;
    bool iselenoiso80=false;
    bool iselenoiso90=false;
    float elemva=-2;
    float elemva_HZZ=-2;
    bool isconversionveto=false;
    int elemissinghits = 999;
    int elemissinglosthits = 999;
    bool iselechargeconsistent=false;

    int decay=-1;
    float ieta=-1,full5x5_ieta=-1,hOverE=-1,etasuperatvtx=-1,phisuperatvtx=-1,IoEmIoP=-999.,IoEmIoP_ttH=-999.,depositTracker=-1,depositEcal=-1,depositHcal=-1,SCeta=-999.;
    int decayModeFindingOldDMs=-1, decayModeFindingNewDMs=-1; // tau 13 TeV ID
    float byCombinedIsolationDeltaBetaCorrRaw3Hits=-1., chargedIsoPtSum=-1., neutralIsoPtSum=-1., puCorrPtSum=-1.; // tau 13 TeV RAW iso info
    int numChargedParticlesSignalCone=-1, numNeutralHadronsSignalCone=-1, numPhotonsSignalCone=-1, numParticlesSignalCone=-1, numChargedParticlesIsoCone=-1, numNeutralHadronsIsoCone=-1, numPhotonsIsoCone=-1, numParticlesIsoCone=-1;
    float leadChargedParticlePt=-1., trackRefPt=-1.;
    int typeOfMuon=0;
    float byIsolationMVArun2v1DBoldDMwLTraw=-1, byIsolationMVArun2017v2DBoldDMwLTraw2017=-1, byIsolationMVArun2017v1DBoldDMwLTraw2017=-1, byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017=-1; //FRA
    int  byVVLooseIsolationMVArun2017v2DBoldDMwLT2017=-1; //FRA
    Long64_t tauIDflag = 0;
    float footprintCorrection, neutralIsoPtSumWeight, photonPtSumOutsideSignalCone;

    float dxy_innerTrack = -1., dz_innerTrack = -1., sip = -1., error_trackpt=-1.;
    int jetNDauChargedMVASel = -1;
    float miniRelIsoCharged = -1., miniRelIsoNeutral = -1.;
    float jetPtRel = -1., jetPtRatio = -1., jetBTagCSV=-1., jetBTagDeepCSV=-1.;
    float lepMVA_mvaId = -1.;

    //
    GlobalPoint aPVPoint(_pv_x, _pv_y, _pv_z);
    GlobalPoint aPVRefitPoint(_pvRefit_x, _pvRefit_y, _pvRefit_z);    
    GlobalPoint aPVGenPoint(_pvGen_x, _pvGen_y, _pvGen_z);

    TVector3 pcaPV = getPCA(event, setup, cand->bestTrack(), aPVPoint);
    TVector3 pcaRefitPV = getPCA(event, setup, cand->bestTrack(), aPVRefitPoint);
    TVector3 pcaGenPV;
    if(theisMC) pcaGenPV = getPCA(event, setup, cand->bestTrack(), aPVGenPoint);

    if(type==ParticleType::MUON){	
      muIDflag=userdatahelpers::getUserInt(cand,"muonID");
      discr = (float) muIDflag; // not really needed, will use the muonID branch in ntuples...
      if(userdatahelpers::getUserFloat(cand,"isPFMuon"))typeOfMuon |= 1 << 0;
      if(userdatahelpers::getUserFloat(cand,"isGlobalMuon"))typeOfMuon |= 1 << 1;
      if(userdatahelpers::getUserFloat(cand,"isTrackerMuon"))typeOfMuon |= 1 << 2;
      depositTracker=userdatahelpers::getUserFloat(cand,"DepositR03TrackerOfficial");
      depositEcal=userdatahelpers::getUserFloat(cand,"DepositR03ECal");
      depositHcal=userdatahelpers::getUserFloat(cand,"DepositR03Hcal");

      dxy_innerTrack = userdatahelpers::getUserFloat(cand,"dxy_innerTrack");
      dz_innerTrack = userdatahelpers::getUserFloat(cand,"dz_innerTrack");
      error_trackpt = userdatahelpers::getUserFloat(cand,"rel_error_trackpt");
      sip = userdatahelpers::getUserFloat(cand,"SIP");

      jetNDauChargedMVASel= LeptonIsoHelper::jetNDauChargedMVASel(cand, closest_jet);
      std::pair<float,float> miniRelIso = LeptonIsoHelper::miniRelIso_ChargedNeutral(cand, pfCands_charged, pfCands_neutral, rho_miniRelIso);
      miniRelIsoCharged = miniRelIso.first;
      miniRelIsoNeutral = miniRelIso.second;

      jetPtRel = LeptonIsoHelper::jetPtRel(*cand, closest_jet,theJECName);
      jetPtRatio = LeptonIsoHelper::jetPtRatio(*cand, closest_jet,theJECName);
      jetBTagCSV = closest_jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
      jetBTagDeepCSV = closest_jet.bDiscriminator("pfDeepCSVJetTags:probb") + closest_jet.bDiscriminator("pfDeepCSVJetTags:probbb");

      lepMVA_mvaId  = userdatahelpers::getUserFloat(cand,"segmentCompatibility");

    }else if(type==ParticleType::ELECTRON){
      discr=userdatahelpers::getUserFloat(cand,"BDT");
      ieta=userdatahelpers::getUserFloat(cand,"sigmaIetaIeta");
      full5x5_ieta=userdatahelpers::getUserFloat(cand,"full5x5_sigmaIetaIeta");
      hOverE=userdatahelpers::getUserFloat(cand,"hOverE");
      etasuperatvtx=userdatahelpers::getUserFloat(cand,"deltaEtaSuperClusterTrackAtVtx");
      phisuperatvtx=userdatahelpers::getUserFloat(cand,"deltaPhiSuperClusterTrackAtVtx");
      IoEmIoP=userdatahelpers::getUserFloat(cand,"IoEmIoP");
      IoEmIoP_ttH=userdatahelpers::getUserFloat(cand,"IoEmIoP_ttH");
      SCeta = userdatahelpers::getUserFloat(cand,"SCeta");
      if(userdatahelpers::getUserInt(cand,"isBDT") == 1)isgood=true;
      if(userdatahelpers::getUserInt(cand,"isEleIDLoose") == 1) iseleLoose=true;
      if(userdatahelpers::getUserInt(cand,"isEleID80") == 1) isele80=true;
      if(userdatahelpers::getUserInt(cand,"isEleID90") == 1) isele90=true;
      if(userdatahelpers::getUserInt(cand,"isEleNoIsoIDLoose") == 1) iselenoisoLoose=true;
      if(userdatahelpers::getUserInt(cand,"isEleNoIsoID80") == 1) iselenoiso80=true;
      if(userdatahelpers::getUserInt(cand,"isEleNoIsoID90") == 1) iselenoiso90=true;
      elemva=(userdatahelpers::getUserFloat(cand,"eleMVAvalue"));
      elemva_HZZ=(userdatahelpers::getUserFloat(cand,"HZZeleMVAvalue"));
      if(userdatahelpers::getUserInt(cand,"isConversionVeto") == 1)isconversionveto=true;
      error_trackpt = userdatahelpers::getUserFloat(cand,"rel_error_trackpt");
      elemissinghits = userdatahelpers::getUserInt(cand,"missingHit");
      elemissinglosthits = userdatahelpers::getUserInt(cand,"missingLostHit");
      if((userdatahelpers::getUserInt(cand,"isGsfCtfScPixChargeConsistent") + userdatahelpers::getUserInt(cand,"isGsfScPixChargeConsistent"))>1)iselechargeconsistent=true;

      //if(userdatahelpers::getUserInt(cand,"isCUT"))isgoodcut=true;

      sip = userdatahelpers::getUserFloat(cand,"SIP");
      jetNDauChargedMVASel= LeptonIsoHelper::jetNDauChargedMVASel(cand, closest_jet);
      std::pair<float,float> miniRelIso = LeptonIsoHelper::miniRelIso_ChargedNeutral(cand, pfCands_charged, pfCands_neutral, rho_miniRelIso);
      miniRelIsoCharged = miniRelIso.first;
      miniRelIsoNeutral = miniRelIso.second;
      
      jetPtRel = LeptonIsoHelper::jetPtRel(*cand, closest_jet,theJECName);
      jetPtRatio = LeptonIsoHelper::jetPtRatio(*cand, closest_jet,theJECName);
      jetBTagCSV = closest_jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
      jetBTagDeepCSV = closest_jet.bDiscriminator("pfDeepCSVJetTags:probb") + closest_jet.bDiscriminator("pfDeepCSVJetTags:probbb");

      lepMVA_mvaId = elemva_HZZ;

    }else if(type==ParticleType::TAU){
      discr=userdatahelpers::getUserFloat(cand,"HPSDiscriminator");
      decay = userdatahelpers::getUserFloat(cand,"decayMode");
      decayModeFindingOldDMs = userdatahelpers::getUserInt (cand, "decayModeFinding");
      decayModeFindingNewDMs = userdatahelpers::getUserInt (cand, "decayModeFindingNewDMs");
      for (uint itau =0; itau<ntauIds; itau++){
        int id = userdatahelpers::getUserInt (cand,  tauIDStrings[itau]);
        if(id>0){
          tauIDflag |= (Long64_t(1) << itau);
          hTauIDs->Fill(id);
        }
      }
      footprintCorrection = userdatahelpers::getUserFloat (cand, "footprintCorrection");
      neutralIsoPtSumWeight = userdatahelpers::getUserFloat (cand, "neutralIsoPtSumWeight");
      photonPtSumOutsideSignalCone = userdatahelpers::getUserFloat (cand, "photonPtSumOutsideSignalCone");

      byCombinedIsolationDeltaBetaCorrRaw3Hits = userdatahelpers::getUserFloat (cand, "byCombinedIsolationDeltaBetaCorrRaw3Hits");
      byIsolationMVArun2v1DBoldDMwLTraw=userdatahelpers::getUserFloat (cand, "byIsolationMVArun2v1DBoldDMwLTraw");
      byIsolationMVArun2017v2DBoldDMwLTraw2017=userdatahelpers::getUserFloat (cand, "byIsolationMVArun2017v2DBoldDMwLTraw2017"); //FRA
      byIsolationMVArun2017v1DBoldDMwLTraw2017=userdatahelpers::getUserFloat (cand, "byIsolationMVArun2017v1DBoldDMwLTraw2017"); //FRA
      byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017=userdatahelpers::getUserFloat (cand, "byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017"); //FRA
      byVVLooseIsolationMVArun2017v2DBoldDMwLT2017= userdatahelpers::getUserInt (cand, "byVVLooseIsolationMVArun2017v2DBoldDMwLT2017"); //FRA
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

      const pat::Tau *taon  = dynamic_cast<const pat::Tau*>(cand);
      if(taon){
	pcaPV = getPCA(event, setup, taon->leadChargedHadrCand()->bestTrack(), aPVPoint);
	pcaRefitPV = getPCA(event, setup, taon->leadChargedHadrCand()->bestTrack(), aPVRefitPoint);
	if(theisMC) pcaGenPV = getPCA(event, setup, taon->leadChargedHadrCand()->bestTrack(), aPVGenPoint);

	reco::CandidatePtrVector chCands = taon->signalChargedHadrCands();
	reco::CandidatePtrVector neCands = taon->signalGammaCands();
	chargedP4 = math::XYZTLorentzVector();
	neutralP4 = math::XYZTLorentzVector();
	for(reco::CandidatePtrVector::const_iterator id=chCands.begin();id!=chCands.end(); ++id) chargedP4 += (*id)->p4();
	for(reco::CandidatePtrVector::const_iterator id=neCands.begin();id!=neCands.end(); ++id) neutralP4 += (*id)->p4();
      }
    }
    _discriminator.push_back(discr);
    _daughters_typeOfMuon.push_back(typeOfMuon);
    _daughters_muonID.push_back(muIDflag);
    _daughters_tauID.push_back(tauIDflag);
    _daughters_footprintCorrection.push_back(footprintCorrection);
    _daughters_neutralIsoPtSumWeight.push_back(neutralIsoPtSumWeight);
    _daughters_photonPtSumOutsideSignalCone.push_back(photonPtSumOutsideSignalCone);

    _daughters_charge.push_back(cand->charge());
    _daughters_iseleBDT.push_back(isgood);
    _daughters_iseleWPLoose.push_back(iseleLoose);
    _daughters_iseleWP80.push_back(isele80);
    _daughters_iseleWP90.push_back(isele90);
    _daughters_iseleNoIsoWPLoose.push_back(iselenoisoLoose);
    _daughters_iseleNoIsoWP80.push_back(iselenoiso80);
    _daughters_iseleNoIsoWP90.push_back(iselenoiso90);
    _daughters_eleMVAnt.push_back(elemva); 
    _daughters_eleMVA_HZZ.push_back(elemva_HZZ);     
    _daughters_passConversionVeto.push_back(isconversionveto);
    _daughters_eleMissingHits.push_back(elemissinghits);
    _daughters_eleMissingLostHits.push_back(elemissinglosthits);
    _daughters_iseleChargeConsistent.push_back(iselechargeconsistent);

    //_daughters_iseleCUT.push_back(userdatahelpers::getUserInt(cand,"isCUT"));
    _decayType.push_back(decay);
    _daughters_IetaIeta.push_back(ieta);
    _daughters_full5x5_IetaIeta.push_back(full5x5_ieta);
    _daughters_hOverE.push_back(hOverE);
    _daughters_deltaEtaSuperClusterTrackAtVtx.push_back(etasuperatvtx);
    _daughters_deltaPhiSuperClusterTrackAtVtx.push_back(phisuperatvtx);
    _daughters_IoEmIoP.push_back(IoEmIoP);
    _daughters_IoEmIoP_ttH.push_back(IoEmIoP_ttH);
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
    _daughters_byIsolationMVArun2v1DBoldDMwLTraw.push_back(byIsolationMVArun2v1DBoldDMwLTraw);
    _daughters_byIsolationMVArun2017v2DBoldDMwLTraw2017.push_back(byIsolationMVArun2017v2DBoldDMwLTraw2017); //FRA
    _daughters_byIsolationMVArun2017v1DBoldDMwLTraw2017.push_back(byIsolationMVArun2017v1DBoldDMwLTraw2017); //FRA
    _daughters_byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017.push_back(byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017); //FRA
    _daughters_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017.push_back(byVVLooseIsolationMVArun2017v2DBoldDMwLT2017); //FRA
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

    _dxy_innerTrack.push_back(dxy_innerTrack);
    _dz_innerTrack.push_back(dz_innerTrack);
    _daughters_rel_error_trackpt.push_back(error_trackpt);
    _SIP.push_back(sip);

    _daughters_jetNDauChargedMVASel.push_back(jetNDauChargedMVASel);
    _daughters_miniRelIsoCharged.push_back(miniRelIsoCharged);
    _daughters_miniRelIsoNeutral.push_back(miniRelIsoNeutral);
    _daughters_jetPtRel.push_back(jetPtRel);
    _daughters_jetPtRatio.push_back(jetPtRatio);
    _daughters_jetBTagCSV.push_back(jetBTagCSV);
    _daughters_jetBTagDeepCSV.push_back(jetBTagDeepCSV);
    _daughters_lepMVA_mvaId.push_back(lepMVA_mvaId);

    _daughters_pca_x.push_back(pcaPV.X());
    _daughters_pca_y.push_back(pcaPV.Y());
    _daughters_pca_z.push_back(pcaPV.Z());

    _daughters_pcaRefitPV_x.push_back(pcaRefitPV.X());
    _daughters_pcaRefitPV_y.push_back(pcaRefitPV.Y());
    _daughters_pcaRefitPV_z.push_back(pcaRefitPV.Z());

    _daughters_pcaGenPV_x.push_back(pcaGenPV.X());
    _daughters_pcaGenPV_y.push_back(pcaGenPV.Y());
    _daughters_pcaGenPV_z.push_back(pcaGenPV.Z());

    _daughters_charged_px.push_back(chargedP4.X());
    _daughters_charged_py.push_back(chargedP4.Y());
    _daughters_charged_pz.push_back(chargedP4.Z());
    _daughters_charged_e.push_back(chargedP4.T());

    _daughters_neutral_px.push_back(neutralP4.X());
    _daughters_neutral_py.push_back(neutralP4.Y());
    _daughters_neutral_pz.push_back(neutralP4.Z());
    _daughters_neutral_e.push_back(neutralP4.T());

    //TRIGGER MATCHING
    Long64_t LFtriggerbit=0,L3triggerbit=0,filterFired=0;
    Long64_t trgMatched = 0;
    Long64_t triggertypeIsGood = 0;
    float hltpt=0;
    
    // list of indexes of all TO standalone that are matched to this specific daughter and pass HLT filter(s)
    // use as: toStandaloneMatched.at(idxHLT).at(idx tostadalone)
    vector<vector<int>> toStandaloneMatched (myTriggerHelper->GetNTriggers(), vector<int>(0));

    // for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
    for (size_t idxto = 0; idxto < triggerObjects->size(); ++idxto)
    {
      pat::TriggerObjectStandAlone obj = triggerObjects->at(idxto);

      // unpacking Filter Labels
      obj.unpackFilterLabels(event,*triggerBits); //FRA
      
      //check if the trigger object matches cand
      bool triggerType=false;

      if(deltaR2(obj,*cand)<0.25){
      
	    //cout << "######### NEW OBJECT MATCHED to offline " << cand->pdgId() << " of pt = " << cand->pt() << " HLT obj pt " << obj.pt() << endl;

        if (type==ParticleType::TAU && (obj.hasTriggerObjectType(trigger::TriggerTau)|| obj.hasTriggerObjectType(trigger::TriggerL1TauJet)))triggerType=true;
        if (type==ParticleType::ELECTRON && (obj.hasTriggerObjectType(trigger::TriggerElectron) || obj.hasTriggerObjectType(trigger::TriggerPhoton)))triggerType=true;
        if (type==ParticleType::MUON && (obj.hasTriggerObjectType(trigger::TriggerMuon)))triggerType=true;
        
        //check fired paths
        obj.unpackPathNames(names);
        std::vector<std::string> pathNamesAll  = obj.pathNames(false);
        std::vector<std::string> pathNamesLast = obj.pathNames(true);

        for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {

          int triggerbit = myTriggerHelper->FindTriggerNumber(pathNamesAll[h],true);
          if (triggerbit < 0) continue ; // not a path I want to save
          bool isLF   = obj.hasPathName( pathNamesAll[h], true, false );
          bool isL3   = obj.hasPathName( pathNamesAll[h], false, true );

          triggerMapper trgmap = myTriggerHelper->GetTriggerMap(pathNamesAll[h]);
          bool isfilterGood = true;
          int IDsearch = 0;
          if (type==ParticleType::ELECTRON) IDsearch = 11;
          else if (type==ParticleType::MUON) IDsearch = 13;
          else if(type==ParticleType::TAU) IDsearch = 15;
          int legPosition = trgmap.GetLegFromID(IDsearch);

          if (legPosition == 1)
          {
            for(int ifilt=0;ifilt<trgmap.GetNfiltersleg1();ifilt++)
            {
              string label = trgmap.Getfilter(true,ifilt);
              if (label.empty()) continue;
              if(! obj.hasFilterLabel(label.c_str()))isfilterGood=false;
            }
          }
          else if (legPosition == 2)
          {
            for(int ifilt=0;ifilt<trgmap.GetNfiltersleg2();ifilt++)
            {
              string label = trgmap.Getfilter(false,ifilt);
              if (label.empty()) continue;
              if(! obj.hasFilterLabel(label.c_str()))isfilterGood=false;
            }
          }
          else isfilterGood = false;
          
          //_isFilterFiredLast;
          if(isfilterGood)filterFired |= long(1) <<triggerbit;
          if(triggerType) triggertypeIsGood |= long(1) << triggerbit;
          if(isLF)LFtriggerbit |= long(1) <<triggerbit;
          if(isL3)L3triggerbit |= long(1) <<triggerbit;
        } // loop on all trigger paths

        // -------------- now do matching "filter-wise" to do x-check
        // trigger matching checking labels
        const std::vector<std::string>& vLabels = obj.filterLabels();
        for (int triggerbit = 0; triggerbit < myTriggerHelper->GetNTriggers(); ++triggerbit)
        {
          triggerMapper trgmap = myTriggerHelper->GetTriggerMap(triggerbit);
            
          bool istrgMatched = true;
          int IDsearch = 0;
          if (type==ParticleType::ELECTRON)  IDsearch = 11;
          else if (type==ParticleType::MUON) IDsearch = 13;
          else if(type==ParticleType::TAU)   IDsearch = 15;
          int legPosition = trgmap.GetLegFromID(IDsearch);

          // debug
          // cout << "***** searching trigger : " << myTriggerHelper -> printTriggerName(triggerbit) << " " << trgmap.GetHLTPath() << endl;
          // cout << "all this object labels: ID " << IDsearch << " --> leg position : " << legPosition << endl;
          // cout << "Nfilters . 1 : " << trgmap.GetNfiltersleg1() << " || 2 : " << trgmap.GetNfiltersleg2() << endl;
          // for (uint ll = 0; ll < vLabels.size(); ++ll) cout << "   -- " << vLabels.at(ll) << endl;

          if (legPosition == 1)
          {
            for(int ifilt=0;ifilt<trgmap.GetNfiltersleg1();ifilt++)
            {
              string label = trgmap.Getfilter(true,ifilt);
              // cout << " @@ leg 1 looking for " << label << endl;
              if (label.empty()) continue;
              if (find(vLabels.begin(), vLabels.end(), label) == vLabels.end()) istrgMatched=false;
            }
          }
          else if (legPosition == 2)
          {
            for(int ifilt=0;ifilt<trgmap.GetNfiltersleg2();ifilt++)
            {
              string label = trgmap.Getfilter(false,ifilt);
              // cout << " @@ leg 2 looking for " << label << endl;
              if (label.empty()) continue;
              if (find(vLabels.begin(), vLabels.end(), label) == vLabels.end()) istrgMatched=false;
            }
          }
          else istrgMatched = false;
          
          // // Check the pT of the candidate for leg1 and 2 //FRA
          //if (legPosition == 1)
          //{
          //  if ( cand->pt() < trgmap.GetPtCut1() ) istrgMatched=false;
          //}
          //else if (legPosition == 2)
          //{
          //  if ( cand->pt() < trgmap.GetPtCut2() ) istrgMatched=false;
          //}
          //else
          //  istrgMatched=false;

          // FIXME: should I check type? --> no, multiple filters should be enough
          if(istrgMatched)
          {
            //std::cout << " ###### GOOD MATCH ######" << std::endl; //FRA
            //std::cout << "***** searching trigger : " << myTriggerHelper -> printTriggerName(triggerbit) << " " << trgmap.GetHLTPath() << std::endl; //FRA
            //std::cout << "--> triggerbit: " << triggerbit << std::endl; //FRA
            //std::cout << "--> label: " << label << std::endl; //FRA
            //std::cout << "--> BEFORE trgMatched: " << std::bitset<64>(trgMatched) << std::endl; //FRA
            
            //printing TO labels //FRA
            //for (uint ll = 0; ll < vLabels.size(); ++ll) cout << "    -- " << vLabels.at(ll) << endl; //FRA
            
            trgMatched |= (long(1) <<triggerbit);
            toStandaloneMatched.at(triggerbit).push_back(idxto);
            //std::cout << "--> AFTER  trgMatched: " << std::bitset<64>(trgMatched) << std::endl; //FRA
          }
          // cout << "istrgMatched ? " << istrgMatched << endl;

        } // loop on triggerbit from 0 to GetNTriggers()

      } // if dR < 0.25
    } // loop on all trigger candidates
    _daughters_isGoodTriggerType.push_back(triggertypeIsGood);
    _daughters_FilterFired.push_back(filterFired);
    _daughters_L3FilterFired.push_back(LFtriggerbit);
    _daughters_L3FilterFiredLast.push_back(L3triggerbit);    
    _daughters_trgMatched.push_back(trgMatched);    
    _daughters_HLTpt.push_back(hltpt);

    vector<int> vTrgMatchedIdx;
    for (int idxHLT = 0; idxHLT < myTriggerHelper->GetNTriggers(); ++idxHLT)
    {
      // if I have more than 1 match, I will be sure that different hlt objects are matched in a pair,
      // so I don't care and I put a value of -1
      // if I have no match, I don't care so I put -1 as well
      if (toStandaloneMatched.at(idxHLT).size() != 1) 
        vTrgMatchedIdx.push_back(-1);
      // if I have exactly 1 match, I store the Trigger Object index.
      // I will compare it later for the two legs of the pair and see if I matched separate objects
      else
        vTrgMatchedIdx.push_back(toStandaloneMatched.at(idxHLT).at(0));
    }
    vTrgMatchedToDau_idx.push_back(vTrgMatchedIdx);


    // L1 candidate matching -- to correct for the missing seed
    bool isL1IsoTauMatched = false;
    /*if(L1ExtraIsoTau.isValid()){
      for (unsigned int iL1IsoTau = 0; iL1IsoTau < L1ExtraIsoTau->size(); iL1IsoTau++)
      {
        const l1extra::L1JetParticle& L1IsoTau = (*L1ExtraIsoTau).at(iL1IsoTau);
        // cout << "IL1TauL: " << iL1IsoTau << " - " << L1IsoTau.pt() << " " << L1IsoTau.eta() << " " << L1IsoTau.phi() << " " << isL1IsoTauMatched << endl;
        // 0.5 cone match + pT requirement as in data taking
	cout<<"iso tau pt "<<L1IsoTau.pt()<<" eta "<<L1IsoTau.eta()<<" "<<L1IsoTau.phi()<<endl;
        if(L1IsoTau.pt() > 32 && deltaR2(L1IsoTau,*cand)<0.25)
          {

            isL1IsoTauMatched = true;
            break;
          }
      }
      }*/
    if (theisMC)
      {
	std::vector<Float_t> L1IsoTau_et; 
	for (int ibx = l1taus->getFirstBX(); ibx <= l1taus->getLastBX(); ++ibx)
	  {
	    for (BXVector<l1t::Tau>::const_iterator it=l1taus->begin(ibx); it!=l1taus->end(ibx); it++)
	      {
		if (it->et() > 0 && ibx ==0){
		  if (it->hwIso() > 0.5){
		    TLorentzVector tlv_L1Tau;
		    TLorentzVector tlv_Tau;
		    tlv_L1Tau.SetPtEtaPhiM(it->et(),
					   it->eta(),
					   it->phi(),
					   0.);
		    
		    tlv_Tau.SetPtEtaPhiM(cand->pt(),
					 cand->eta(),
					 cand->phi(),
					 0.);
		    
		    if ((tlv_L1Tau.DeltaR(tlv_Tau)*tlv_L1Tau.DeltaR(tlv_Tau)) < 0.25) {
		      isL1IsoTauMatched = true;
		      L1IsoTau_et.push_back(it->et());
		    }
		    
		  }
		  
		}
	      }
	    
	  }
	
	if(isL1IsoTauMatched) {
	  std::sort(L1IsoTau_et.begin(), L1IsoTau_et.end());
	  Float_t L1IsoTau_etMax = *L1IsoTau_et.rbegin();
	  _daughters_highestEt_L1IsoTauMatched.push_back(L1IsoTau_etMax) ;
	}
      }
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
    //event.getByLabel ("genInfo", candHandle);
    event.getByToken (theGenericTag, candHandle);
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
	int TauGenDetailedDecayMode = -1;
	TVector3 pca(99,99,99);
        
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
	if (igen->hasUserInt("tauGenDetailedDecayMode"))   TauGenDetailedDecayMode = igen->userInt("tauGenDetailedDecayMode");
	if (igen->hasUserFloat("pca_x")) pca = TVector3(igen->userFloat("pca_x"),igen->userFloat("pca_y"),igen->userFloat("pca_z"));
	  
        
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
	_genpart_TauGenDetailedDecayMode.push_back(TauGenDetailedDecayMode);
	_genpart_pca_x.push_back(pca.X());
	_genpart_pca_y.push_back(pca.Y());
	_genpart_pca_z.push_back(pca.Z());
        
        //const pat::GenericParticle* genClone = &(*igen);
        //int flags = CreateFlagsWord (genClone);
        int flags = igen -> userInt ("generalGenFlags");
        _genpart_flags.push_back(flags);

	if(igen->hasUserInt("HZDecayMode") || igen->hasUserInt("WDecayMode") || igen->hasUserInt("TopDecayMode")){
	  _pvGen_x = igen->vx();
	  _pvGen_y = igen->vy();
	  _pvGen_z = igen->vz();
	}
    }
}


void HTauTauNtuplizer::FillGenJetInfo(const edm::Event& event)
{
    edm::Handle<edm::View<reco::GenJet>> genJetHandle;
    //event.getByLabel ("slimmedGenJets", genJetHandle);
    event.getByToken (theGenJetTag, genJetHandle);
    unsigned int genJetSize = genJetHandle->size();

    // to retrieve gen jet flavour from matched gen jet
    edm::Handle<edm::View<pat::Jet>> patjetHandle;
    //event.getByLabel("jets", patjetHandle);
    event.getByToken(theJetTag, patjetHandle);
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

      // cout << igj << " " << genJet.pdgId() << endl;
      // only mesons and barions in the list, but no B to infer the jet flavour    
      // for (uint ic = 0; ic < genJet.numberOfDaughters(); ++ic)
      // {
      //   const reco::Candidate* cand = genJet.daughter(ic);
      //   cout << ic << " " << cand->pdgId() << " " << cand->pt() << " " << endl;
      // }
    }

    return;

}

// return index of gen matched to reco lepton, and -1 if not existing or not found
int HTauTauNtuplizer::GetMatchedGen (const reco::Candidate* genL, const edm::Event& event)
{
    //cout.precision(15); // just to check real precision
    
    edm::Handle<edm::View<pat::GenericParticle>>candHandle;
    //event.getByLabel ("genInfo", candHandle);
    event.getByToken (theGenericTag, candHandle);
        
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

//Fill L1 objects
void HTauTauNtuplizer::FillL1Obj(const BXVector<l1t::Tau>* taus, const BXVector<l1t::Jet>* jets, const edm::Event& event){

  for (int ibx = taus->getFirstBX(); ibx <= taus->getLastBX(); ++ibx)
    {
      for (BXVector<l1t::Tau>::const_iterator it=taus->begin(ibx); it!=taus->end(ibx); it++)
	{
	  if (it->et() > 0 && ibx ==0){

	    _L1_tauEt .push_back(it->et());
	    _L1_tauEta.push_back(it->eta());
	    _L1_tauPhi.push_back(it->phi());
	    _L1_tauIso.push_back(it->hwIso());
	  }
	}
    }

  for (int ibx = jets->getFirstBX(); ibx <= jets->getLastBX(); ++ibx)
    {
      for (BXVector<l1t::Jet>::const_iterator it=jets->begin(ibx); it!=jets->end(ibx); it++)
	{
	  if (it->et() > 0&& ibx ==0){
	    
	    _L1_jetEt .push_back(it->et());
	    _L1_jetEta.push_back(it->eta());
	    _L1_jetPhi.push_back(it->phi());
	   
	  }
	}
    }
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
void HTauTauNtuplizer::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const& iSetup){
  if (theisMC)
  {
    // try {event.getByToken(theLHEPTag, lheEventProduct);} catch (...) {;}
    // if (lheEventProduct.isValid())
    edm::Handle<GenLumiInfoHeader> gen_header;
    iLumi.getByToken(genLumiHeaderTag, gen_header);
    if (gen_header.isValid())
    { 
      string model = gen_header->configDescription(); 
      // cout<<model<<endl;  // prints, e.g. T1tttt_1500_100
      _susyModel = model;
    }
  }
}
void HTauTauNtuplizer::endLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const& iSetup)
{
  // Total number of events is the sum of the events in each of these luminosity blocks
  edm::Handle<edm::MergeableCounter> nEventsTotalCounter;
  iLumi.getByToken(theTotTag, nEventsTotalCounter);
  Nevt_Gen += nEventsTotalCounter->value;

  edm::Handle<edm::MergeableCounter> nEventsPassTrigCounter;
  iLumi.getByToken(thePassTag, nEventsPassTrigCounter);
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
  //if (!i.daughter(cand1i)->isMuon() && !i.daughter(cand1i)->isElectron()) isoi= -userdatahelpers::getUserFloat(i.daughter(cand1i),"byIsolationMVArun2v1DBoldDMwLTraw");
  //if (!j.daughter(cand1j)->isMuon() && !j.daughter(cand1j)->isElectron()) isoj= -userdatahelpers::getUserFloat(j.daughter(cand1j),"byIsolationMVArun2v1DBoldDMwLTraw");
  if (!i.daughter(cand1i)->isMuon() && !i.daughter(cand1i)->isElectron()) isoi= -userdatahelpers::getUserFloat(i.daughter(cand1i),"byIsolationMVArun2017v2DBoldDMwLTraw2017");
  if (!j.daughter(cand1j)->isMuon() && !j.daughter(cand1j)->isElectron()) isoj= -userdatahelpers::getUserFloat(j.daughter(cand1j),"byIsolationMVArun2017v2DBoldDMwLTraw2017");

  if (isoi<isoj)return true;
  else if(isoi>isoj)return false;

  //step 2, leg 1 Pt
  if(i.daughter(cand1i)->pt()>j.daughter(cand1j)->pt()) return true;
  else if(i.daughter(cand1i)->pt()<j.daughter(cand1j)->pt()) return false;

  //step 3, leg 2 ISO
  isoi=userdatahelpers::getUserFloat(i.daughter(1-cand1i),"combRelIsoPF");
  isoj=userdatahelpers::getUserFloat(j.daughter(1-cand1j),"combRelIsoPF");
  //if (!i.daughter(1-cand1i)->isMuon() && !i.daughter(1-cand1i)->isElectron()) isoi= -userdatahelpers::getUserFloat(i.daughter(1-cand1i),"byIsolationMVArun2v1DBoldDMwLTraw");
  //if (!j.daughter(1-cand1j)->isMuon() && !j.daughter(1-cand1j)->isElectron()) isoj= -userdatahelpers::getUserFloat(j.daughter(1-cand1j),"byIsolationMVArun2v1DBoldDMwLTraw");
  if (!i.daughter(1-cand1i)->isMuon() && !i.daughter(1-cand1i)->isElectron()) isoi= -userdatahelpers::getUserFloat(i.daughter(1-cand1i),"byIsolationMVArun2017v2DBoldDMwLTraw2017");
  if (!j.daughter(1-cand1j)->isMuon() && !j.daughter(1-cand1j)->isElectron()) isoj= -userdatahelpers::getUserFloat(j.daughter(1-cand1j),"byIsolationMVArun2017v2DBoldDMwLTraw2017");

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
    math::XYZTLorentzVector METP4 (METx, METy, 0, 0); // I only care about transverse plane
    double dphi = deltaPhi(visP4.phi(), METP4.phi());
    return sqrt(2.*visP4.Pt()*METP4.Pt()*(1.-cos(dphi)));

    // this code is equivalent, but can give NaN and rounding errors when met and lepton are very collinear
    // float MET = TMath::Sqrt (METx*METx + METy*METy);
    // math::XYZTLorentzVector METP4 (METx, METy, 0, MET);
    
    // float scalSum = MET + visP4.pt();
    // math::XYZTLorentzVector vecSum (visP4);
    // vecSum += METP4;
    // float vecSumPt = vecSum.pt();
    
    // return TMath::Sqrt (scalSum*scalSum - vecSumPt*vecSumPt);
}

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
bool HTauTauNtuplizer::refitPV(const edm::Event & iEvent, const edm::EventSetup & iSetup){

  edm::ESHandle<TransientTrackBuilder> transTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transTrackBuilder);

  edm::Handle<edm::View<pat::PackedCandidate> >pfCandHandle;
  iEvent.getByToken(thePFCandTag,pfCandHandle);
  const edm::View<pat::PackedCandidate>* cands = pfCandHandle.product();

  Handle<vector<reco::Vertex> >  vertices;
  iEvent.getByToken(theVtxTag,vertices);

  edm::Handle<reco::BeamSpot> beamSpot;
  iEvent.getByToken(beamSpotTag, beamSpot);
 
  TransientVertex transVtx;

  //Get tracks associated wiht pfPV
  reco::TrackCollection pvTracks;
  TLorentzVector aTrack;
  for(size_t i=0; i<cands->size(); ++i){
    if((*cands)[i].charge()==0 || (*cands)[i].vertexRef().isNull()) continue;
    if(!(*cands)[i].bestTrack()) continue;
    //if(!(*cands)[i].trackHighPurity()) continue;  // FRA
    //if((*cands)[i].pt()<=0.5) continue;           // FRA
    
    unsigned int key = (*cands)[i].vertexRef().key();
    int quality = (*cands)[i].pvAssociationQuality();

    if(key!=0 ||
       (quality!=pat::PackedCandidate::UsedInFitTight
	&& quality!=pat::PackedCandidate::UsedInFitLoose)) continue;

    pvTracks.push_back(*((*cands)[i].bestTrack()));
  }
  ///Built transient tracks from tracks.
  std::vector<reco::TransientTrack> transTracks;  
  for(auto iter: pvTracks) transTracks.push_back(transTrackBuilder->build(iter));

  bool fitOk = false;  
  if(transTracks.size() >= 2 ) {
    AdaptiveVertexFitter avf;
    avf.setWeightThreshold(0.001); 
    try {
      transVtx = avf.vertex(transTracks, *beamSpot);
      fitOk = true; 
    } catch (...) {
      fitOk = false; 
      std::cout<<"Vtx fit failed!"<<std::endl;
    }
  }

  fitOk = fitOk && transVtx.isValid() && fabs(transVtx.position().x())<1 && fabs(transVtx.position().y())<1;
  
  if(fitOk) {
    ///NOTE: we take original vertex z position, as this gives the best reults on CP
    ///variables. To be understood; probable reason are missing tracks with Pt<0.95GeV
    _pvRefit_x = transVtx.position().x();
    _pvRefit_y = transVtx.position().y();
    //_pvRefit_z = transVtx.position().z();
    _pvRefit_z = (*vertices)[0].z();
  }
  else {
    _pvRefit_x = (*vertices)[0].x();
    _pvRefit_y = (*vertices)[0].y();
    _pvRefit_z = (*vertices)[0].z();
  }

  return fitOk;
}
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
bool HTauTauNtuplizer::findPrimaryVertices(const edm::Event & iEvent, const edm::EventSetup & iSetup){

  Handle<vector<reco::Vertex> >  vertices;
  iEvent.getByToken(theVtxTag,vertices);
  
  if(vertices->size()==0) return false;   //at least one vertex

  _pv_x = (*vertices)[0].x();
  _pv_y = (*vertices)[0].y();
  _pv_z = (*vertices)[0].z();

  _isRefitPV = refitPV(iEvent, iSetup);

  return true;
}
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
TVector3 HTauTauNtuplizer::getPCA(const edm::Event & iEvent, const edm::EventSetup & iSetup,
				  const reco::Track *aTrack,	   
				  const GlobalPoint & aPoint){
  TVector3 aPCA;
  if(!doCPVariables || !aTrack ||  _npv==0 || aTrack->pt()<2) return aPCA;

  edm::ESHandle<TransientTrackBuilder> transTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transTrackBuilder);
  if(!transTrackBuilder.isValid()){
    std::cout<<"Problem with TransientTrackBuilder"<<std::endl;
    return aPCA;
  }

  reco::TransientTrack transTrk=transTrackBuilder->build(aTrack);
  
  AnalyticalImpactPointExtrapolator extrapolator(transTrk.field());
  GlobalPoint pos  = extrapolator.extrapolate(transTrk.impactPointState(),aPoint).globalPosition();

  aPCA.SetX(pos.x() - aPoint.x());
  aPCA.SetY(pos.y() - aPoint.y());
  aPCA.SetZ(pos.z() - aPoint.z());

  return aPCA;
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////


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
