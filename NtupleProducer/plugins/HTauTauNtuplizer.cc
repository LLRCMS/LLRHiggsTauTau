/** \class HTauTauNtuplizer
 *
 *  No description available.
 *
 *  $Date: 2014/10/31 10:08:19 $
 *  $Revision: 1.00 $
 *  \author G. Ortona (LLR)
 */

// system include files
#include <memory>
#include <cmath>
#include <vector>
#include <string>
#include <TNtuple.h>
//#include <XYZTLorentzVector.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/Framework/interface/LuminosityBlock.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Common/interface/TriggerNames.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>

#include <DataFormats/Common/interface/TriggerResults.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
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

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include <CommonTools/UtilAlgos/interface/TFileService.h>

#include <Muon/MuonAnalysisTools/interface/MuonEffectiveArea.h>

#include <LLRHiggsTauTau/NtupleProducer/interface/CutSet.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/LeptonIsoHelper.h>
//#include <LLRHiggsTauTau/NtupleProducer/interface/DaughterDataHelpers.h>
//#include <LLRHiggsTauTau/NtupleProducer/interface/FinalStates.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/MCHistoryTools.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/PUReweight.h>
//#include <LLRHiggsTauTau/NtupleProducer/interface/VBFCandidateJetSelector.h>
//#include <LLRHiggsTauTau/NtupleProducer/interface/bitops.h>
//#include <LLRHiggsTauTau/NtupleProducer/interface/Fisher.h>
//#include <LLRHiggsTauTau/NtupleProducer/interface/HTauTauConfigHelper.h>
//#include "HZZ4lNtupleFactory.h"

#include "TLorentzVector.h"

// namespace {
//   bool writePhotons = false;  // Write photons in the tree. 
//   bool writeJets = true;     // Write jets in the tree. 
// }


using namespace std;
using namespace edm;
using namespace reco;

// class declaration


class HTauTauNtuplizer : public edm::EDAnalyzer {
 public:
  /// Constructor
  explicit HTauTauNtuplizer(const edm::ParameterSet&);
    
  /// Destructor
  ~HTauTauNtuplizer();  

  TString makeVarString();

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
  //----To implement here-----
  //virtual void FillCandidate(const pat::CompositeCandidate& higgs, bool evtPass, const edm::Event&, const Int_t CRflag);
  //virtual void FillPhoton(const pat::Photon& photon);
  //virtual void FillJet(const cmg::PFJet& jet);

  // ----------member data ---------------------------
  //HTauTauConfigHelper myHelper;
  int theChannel;
  std::string theCandLabel;
  TString theFileName;

  //Output Objects
  TTree *myTree;//->See from ntuplefactory in zz4l
  TH1F *hCounter;

  //flags
  static const int nOutVars =14;
  Bool_t isMC;
  bool applyTrigger;    // Only events passing trigger
  bool applySkim;       //  "     "      "     skim
  bool skipEmptyEvents; // Skip events whith no candidate in the collection
  //PUReweight reweight;

  //counters
  Int_t Nevt_Gen;
  Int_t Npairs;

  //Output variables
  std::vector<math::XYZTLorentzVector> _mothers;
  std::vector<math::XYZTLorentzVector> _daughter1;
  std::vector<math::XYZTLorentzVector> _daughter2;
  std::vector<Int_t> _indexmot;
  Int_t _indexevents;
  Int_t _runNumber;
};

// ----Constructor and Desctructor -----
HTauTauNtuplizer::HTauTauNtuplizer(const edm::ParameterSet& pset) {
  theCandLabel = pset.getUntrackedParameter<string>("CandCollection");
  //theChannel = myHelper.channel();
  theFileName = pset.getUntrackedParameter<string>("fileName");
  skipEmptyEvents = pset.getParameter<bool>("skipEmptyEvents");
  //writeBestCandOnly = pset.getParameter<bool>("onlyBestCandidate");
  //sampleName = pset.getParameter<string>("sampleName");
  Nevt_Gen=0;
  Npairs=0;

  Initialize();
}

HTauTauNtuplizer::~HTauTauNtuplizer(){}
//

void HTauTauNtuplizer::Initialize(){
  _mothers.clear();
  _daughter1.clear();
  _daughter2.clear();
  _indexmot.clear();
  _indexevents=0;
  _runNumber=0;
}

void HTauTauNtuplizer::beginJob(){
  edm::Service<TFileService> fs;
  //TString varString = makeVarString();
  myTree = fs->make<TTree>("HTauTauTree","HTauTauTree");
  hCounter = fs->make<TH1F>("Counters","Counters",2,0,2);

  //Branches
  myTree->Branch("EventNumber",&_indexevents,"EventNumber/I");
  myTree->Branch("RunNumber",&_runNumber,"RunNumber/I");
  myTree->Branch("mothers",&_mothers);
  myTree->Branch("daughters1",&_daughter1);
  myTree->Branch("daughters2",&_daughter2);
  myTree->Branch("indexMothers",&_indexmot);

}

//For semplicity and flexibility, list output variables here. Possibly, amke flags for different
//outputs, but take care of modifying the filling acocrdingly
//This function is obsolete, TBR
TString HTauTauNtuplizer::makeVarString(){
  TString varlist[nOutVars] = {
    "run",
    "iev",
    "P0Pair",
    "P1Pair",
    "P2Pair",
    "P3Pair",
    "P01",
    "P11",
    "P21",
    "P31",
    "P02",
    "P12",
    "P22",
    "P32"
  };
  TString out = varlist[0].Data();
  for(int i=1;i<nOutVars;i++){out.Append(":");out.Append(varlist[i].Data());}
  return out;
}
// ----Analyzer (main) ----

// ------------ method called for each event  ------------
void HTauTauNtuplizer::analyze(const edm::Event& event, const edm::EventSetup& eSetup)
{
  Nevt_Gen++; 
  //Initialize();

  Handle<vector<reco::Vertex> >  vertexs;
  event.getByLabel("goodPrimaryVertices",vertexs);
  
  //----------------------------------------------------------------------
  // Analyze MC history. THIS HAS TO BE DONE BEFORE ANY RETURN STATEMENT
  // (eg skim or trigger), in order to update the gen counters correctly!!!
  
  //  std::vector<const reco::Candidate *> genZs;
  // std::vector<const reco::Candidate *> genZLeps;
  // if (isMC) {
    
  //   Handle<std::vector< PileupSummaryInfo > >  PupInfo;
  //   event.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);
    
  //   std::vector<PileupSummaryInfo>::const_iterator PVI;
  //   for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
  //     if(PVI->getBunchCrossing() == 0) { 
  // 	nObsInt  = PVI->getPU_NumInteractions();
  // 	nTrueInt = PVI->getTrueNumInteractions();
  // 	break;
  //     } 
  //   }
  // }


  //Get candidate collection
  //edm::Handle<edm::View<pat::CompositeCandidate>>candHandle;
  //event.getByLabel("barellCand",candHandle);
  //const edm::View<pat::CompositeCandidate>* cands = candHandle.product();
  edm::Handle<edm::View<reco::CompositeCandidate>>candHandle;
  event.getByLabel("barellCand",candHandle);
  const edm::View<reco::CompositeCandidate>* cands = candHandle.product();

  //myNtuple->InitializeVariables();
  _indexevents = event.id().event();
  _runNumber = event.id().run();

  //Do all the stuff here
  //Compute the variables needed for the output and store them in the ntuple
  int iMot=0;
  for(edm::View<reco::CompositeCandidate>::const_iterator candi = cands->begin(); candi!=cands->end();++candi){
    Npairs++;
    const reco::CompositeCandidate& cand = (*candi);
    _mothers.push_back(cand.p4());
    _daughter1.push_back(cand.daughter(0)->p4());
    _daughter2.push_back(cand.daughter(1)->p4());
    //We need to find a way to avoid repetitions of daughter in order to save space
    _indexmot.push_back(iMot);
  
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
  }
  myTree->Fill();
  //return;
}

void HTauTauNtuplizer::endJob(){
  hCounter->SetBinContent(1,Nevt_Gen);
  hCounter->SetBinContent(2,Npairs);
}

void HTauTauNtuplizer::beginRun(edm::Run const&, edm::EventSetup const&){
}
void HTauTauNtuplizer::endRun(edm::Run const&, edm::EventSetup const&){
}
void HTauTauNtuplizer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){
}
void HTauTauNtuplizer::endLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const& iSetup)
{
  Float_t Nevt_preskim = -1.;
  edm::Handle<edm::MergeableCounter> preSkimCounter;
  if (iLumi.getByLabel("preSkimCounter", preSkimCounter)) { // Counter before skim. Does not exist for non-skimmed samples.
    Nevt_preskim = preSkimCounter->value;
  }  
  if (Nevt_preskim>=0.) {
    Nevt_Gen = Nevt_Gen + Nevt_preskim; 
  }
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
