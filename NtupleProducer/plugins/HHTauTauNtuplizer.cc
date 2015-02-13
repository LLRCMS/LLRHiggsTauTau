/** \class HHTauTauNtuplizer
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

#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
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
//#include "LLRHiggsTauTau/NtupleProducer/Utils/OfflineProducerHelper.h"

#include "TLorentzVector.h"

 namespace {
//   bool writePhotons = false;  // Write photons in the tree. 
   bool writeJets = true;     // Write jets in the tree. 
   //   bool writeSoftLep = false;
   bool DEBUG = false;
 }

using namespace std;
using namespace edm;
using namespace reco;

// class declaration

class HHTauTauNtuplizer : public edm::EDAnalyzer {
 public:
  /// Constructor
  explicit HHTauTauNtuplizer(const edm::ParameterSet&);
    
  /// Destructor
  ~HHTauTauNtuplizer();  

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
  //Int_t FindCandIndex(const reco::Candidate&, Int_t iCand);
  //----To implement here-----
  //virtual void FillCandidate(const pat::CompositeCandidate& higgs, bool evtPass, const edm::Event&, const Int_t CRflag);
  //virtual void FillPhoton(const pat::Photon& photon);
  int FillJet(const edm::View<pat::Jet>* jet);
  void FillbQuarks(const edm::Event&);

  // ----------member data ---------------------------
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
  HLTConfigProvider hltConfig_;
  //Output Objects
  TTree *myTree;//->See from ntuplefactory in zz4l
  TH1F *hCounter;

  //flags
  static const int nOutVars =14;
  bool applyTrigger;    // Only events passing trigger
  bool applySkim;       //  "     "      "     skim
  bool skipEmptyEvents; // Skip events whith no candidate in the collection
  //PUReweight reweight;

  //counters
  Int_t Nevt_Gen;
  Int_t Nevt_HLT;
  
  //Event Output variables
  Int_t _indexevents;
  Int_t _runNumber;
  Int_t _triggerbit;
  Float_t _met;
  Float_t _metphi;
  
  //Mothers output variables
  Float_t _SVmass;
  Float_t _metx;
  Float_t _mety;
  Float_t _mothers_PX;
  Float_t _mothers_PY;
  Float_t _mothers_PZ;
  Float_t _mothers_E;

  //Leptons variables
  //  std::vector<Int_t> _pdgdau;
  Float_t _dau0_PX;
  Float_t _dau0_PY;
  Float_t _dau0_PZ;
  Float_t _dau0_E;
  Float_t _dau0_ISO;
  Float_t _dau0_Discriminator;//BDT for ele, discriminator for tau
  Float_t _dau0_dxy;
  Float_t _dau0_dz;
  Int_t _dau0_decayMode;//for taus only

  Float_t _dau1_PX;
  Float_t _dau1_PY;
  Float_t _dau1_PZ;
  Float_t _dau1_E;
  Float_t _dau1_ISO;
  Float_t _dau1_Discriminator;//BDT for ele, discriminator for tau
  Float_t _dau1_dxy;
  Float_t _dau1_dz;
  Int_t _dau1_decayMode;//for taus only

  //Jets variables
  Int_t _numberOfJets;
  std::vector<math::XYZTLorentzVector> _jets;
  std::vector<Int_t> _jetFlavour;
  std::vector<Float_t> _bdiscr;
  std::vector<Float_t> _bdiscr2;

  //b quarks variables
  std::vector<math::XYZTLorentzVector> _bquarks;
  std::vector<Float_t> _bmotmass;
};

// ----Constructor and Destructor -----
HHTauTauNtuplizer::HHTauTauNtuplizer(const edm::ParameterSet& pset) {
  theCandLabel = pset.getUntrackedParameter<string>("CandCollection");
  //theChannel = myHelper.channel();
  theFileName = pset.getUntrackedParameter<string>("fileName");
  skipEmptyEvents = pset.getParameter<bool>("skipEmptyEvents");
  theFSR = pset.getParameter<bool>("applyFSR");
  theisMC = pset.getParameter<bool>("IsMC");
  //writeBestCandOnly = pset.getParameter<bool>("onlyBestCandidate");
  //sampleName = pset.getParameter<string>("sampleName");
  Nevt_Gen=0;
 
  //triggerResultsLabel = InputTag("TriggerResults","","HLT");
  processName= pset.getParameter<edm::InputTag>("triggerResultsLabel");
  Initialize();
}

HHTauTauNtuplizer::~HHTauTauNtuplizer(){}
//

void HHTauTauNtuplizer::Initialize(){
  _bquarks.clear();
  _bmotmass.clear();
  _SVmass=0;
  _metx=0;
  _mety=0;
  _mothers_PX=0;
  _mothers_PY=0;
  _mothers_PZ=0;
  _mothers_E=0;

  //Leptons variables
  //  std::vector<Int_t> _pdgdau;
  _dau0_PX=0;
  _dau0_PY=0;
  _dau0_PZ=0;
  _dau0_E=0;
  _dau0_ISO=0;
  _dau0_Discriminator=0;//BDT for ele, discriminator for tau
  _dau0_dxy=0;
  _dau0_dz=0;
  _dau0_decayMode=0;//for taus only

  _dau1_PX=0;
  _dau1_PY=0;
  _dau1_PZ=0;
  _dau1_E=0;
  _dau1_ISO=0;
  _dau1_Discriminator=0;//BDT for ele, discriminator for tau
  _dau1_dxy=0;
  _dau1_dz=0;
  _dau1_decayMode=0;//for taus only
  _indexevents=0;
  _runNumber=0;
  _triggerbit=0;
  _met=0;
  _metphi=0.;
  _jets.clear();
  _jetFlavour.clear();
  _numberOfJets=0;
  _bdiscr.clear();
  _bdiscr2.clear();
}

void HHTauTauNtuplizer::beginJob(){
  edm::Service<TFileService> fs;
  myTree = fs->make<TTree>("HTauTauTree","HTauTauTree");
  hCounter = fs->make<TH1F>("Counters","Counters",2,0,2);

  //Branches
  myTree->Branch("EventNumber",&_indexevents,"EventNumber/I");
  myTree->Branch("RunNumber",&_runNumber,"RunNumber/I");
  myTree->Branch("triggerbit",&_triggerbit,"triggerbit/I");
  myTree->Branch("met",&_met,"met/F");
  myTree->Branch("metphi",&_metphi,"metphi/F");  
  if(theisMC){
    //    myTree->Branch("genDaughters",&_genDaughters);
    myTree->Branch("bquarks",&_bquarks);
    myTree->Branch("bmotmass",&_bmotmass);
  }

  myTree->Branch("SVfitMass",&_SVmass,"SVfitMass/F");
  myTree->Branch("METx",&_metx,"MET_x/F");
  myTree->Branch("METy",&_mety,"MET_y/F");
  myTree->Branch("Px",&_mothers_PX,"Px/F");
  myTree->Branch("Py",&_mothers_PY,"Py/F");
  myTree->Branch("Pz",&_mothers_PZ,"Pz/F");
  myTree->Branch("E",&_mothers_E,"E/F");

  myTree->Branch("dau0_Px",&_dau0_PX,"dau0_Px/F");
  myTree->Branch("dau0_Py",&_dau0_PY,"dau0_Py/F");
  myTree->Branch("dau0_Pz",&_dau0_PZ,"dau0_Pz/F");
  myTree->Branch("dau0_E",&_dau0_E,"dau0_E/F");
  myTree->Branch("dau0_Discriminator",&_dau0_Discriminator,"dau0_Discriminator/F");
  myTree->Branch("dau0_dxy",&_dau0_dxy,"dau0_dxy/F");
  myTree->Branch("dau0_dz",&_dau0_dz,"dau0_dz/F");
  myTree->Branch("dau0_decayMode",&_dau0_decayMode,"dau0_decayMode/F");
  myTree->Branch("dau0_combreliso",& _dau0_ISO,"dau0_combreliso/F");

  myTree->Branch("dau1_Px",&_dau1_PX,"dau1_Px/F");
  myTree->Branch("dau1_Py",&_dau1_PY,"dau1_Py/F");
  myTree->Branch("dau1_Pz",&_dau1_PZ,"dau1_Pz/F");
  myTree->Branch("dau1_E",&_dau1_E,"dau1_E/F");
  myTree->Branch("dau1_Discriminator",&_dau1_Discriminator,"dau1_Discriminator/F");
  myTree->Branch("dau1_dxy",&_dau1_dxy,"dau1_dxy/F");
  myTree->Branch("dau1_dz",&_dau1_dz,"dau1_dz/F");
  myTree->Branch("dau1_decayMode",&_dau1_decayMode,"dau1_decayMode/F");
  myTree->Branch("dau1_combreliso",& _dau1_ISO,"dau1_combreliso/F");

  myTree->Branch("JetsNumber",&_numberOfJets,"JetsNumber/I");
  myTree->Branch("jets",&_jets);
  myTree->Branch("jetFlavour",&_jetFlavour);
  myTree->Branch("bDiscriminator",&_bdiscr);
  myTree->Branch("bCSVscore",&_bdiscr2);
}

/*
Int_t HHTauTauNtuplizer::FindCandIndex(const reco::Candidate& cand,Int_t iCand=0){
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
*/
// ----Analyzer (main) ----
// ------------ method called for each event  ------------
void HHTauTauNtuplizer::analyze(const edm::Event& event, const edm::EventSetup& eSetup)
{
  Nevt_Gen++; 
  Initialize();

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

  triggerhelper myTriggerHelper;
  _triggerbit = myTriggerHelper.FindTriggerBit(event,foundPaths,indexOfPath);
  if(_triggerbit == 0)return;

  //Get candidate collection
  edm::Handle<edm::View<pat::CompositeCandidate>>candHandle;
  edm::Handle<edm::View<reco::Candidate>>dauHandle;
  edm::Handle<edm::View<pat::Jet>>jetHandle;
  edm::Handle<pat::METCollection> metHandle;
  event.getByLabel(theCandLabel,candHandle);
  event.getByLabel("jets",jetHandle);
  //  event.getByLabel("softLeptons",dauHandle);
  event.getByLabel("slimmedMETs",metHandle);
  const edm::View<pat::CompositeCandidate>* cands = candHandle.product();
  //  const edm::View<reco::Candidate>* daus = dauHandle.product();
  const edm::View<pat::Jet>* jets = jetHandle.product();
  const pat::MET &met = metHandle->front();
  //myNtuple->InitializeVariables();
 
  Nevt_HLT++;
    
  if(jets->size()<2)return;

 _indexevents = event.id().event();
  _runNumber = event.id().run();
  _met = met.sumEt();
  _metphi = met.phi();

  //Do all the stuff here
  //Compute the variables needed for the output and store them in the ntuple
  if(DEBUG)printf("===New Event===\n");

  //Loop over generated b quarks
  if(theisMC)FillbQuarks(event);

  //Loop of softleptons and fill them
  //FillSoftLeptons(daus,theFSR);

  //Loop on Jets
  _numberOfJets = 0;
  if(writeJets)_numberOfJets = FillJet(jets);

  //Loop on pairs
  int candsInd = 0;
  int tmpInd = 0;
  float svMass = -99999;
  
  if(cands->size()>1){
    //findBestCand
    for(edm::View<pat::CompositeCandidate>::const_iterator candi = cands->begin(); candi!=cands->end();++candi){
      const pat::CompositeCandidate& cand = (*candi); 
      if(abs(cand.userFloat("SVfitMass")-125.0)<abs(svMass-125.0)){
	candsInd=tmpInd;
	svMass=cand.userFloat("SVfitMass");
      }
      tmpInd++;
    }
  }
  tmpInd=0;
  for(edm::View<pat::CompositeCandidate>::const_iterator candi = cands->begin(); candi!=cands->end();++candi){
    if(tmpInd!=candsInd){
      tmpInd++;
      continue;
    }
    const pat::CompositeCandidate& cand = (*candi); 
    math::XYZTLorentzVector candp4 = cand.p4();
    //_SVmass.push_back(cand.userFloat("SVfitMass"));
    //_metx.push_back(cand.userFloat("MEt_px"));
    //_mety.push_back(cand.userFloat("MEt_py"));
    
    _SVmass=svMass;//cand.userFloat("SVfitMass");
    _metx=cand.userFloat("MEt_px");
    _mety=cand.userFloat("MEt_py");
    _mothers_PX=candp4.Px();  
    _mothers_PY=candp4.Py();  
    _mothers_PZ=candp4.Pz();  
    _mothers_E=candp4.E();  
    
    const reco::Candidate *dau0 = cand.daughter(0);
    const reco::Candidate *dau1 = cand.daughter(1);
    _dau0_PX= dau0->p4().Px();
    _dau0_PY= dau0->p4().Py();
    _dau0_PZ= dau0->p4().Pz();
    _dau0_E= dau0->p4().E();
    _dau0_decayMode = (Int_t)(userdatahelpers::getUserFloat(dau0,"decayMode"));
    _dau0_Discriminator= userdatahelpers::getUserFloat(dau0,"HPSDiscriminator");
    _dau0_ISO =userdatahelpers::getUserFloat(dau0,"combRelIsoPF");
    _dau0_dxy =userdatahelpers::getUserFloat(dau0,"dxy");
    _dau0_dz =userdatahelpers::getUserFloat(dau0,"dz");
    
    _dau1_PX= dau1->p4().Px();
    _dau1_PY= dau1->p4().Py();
    _dau1_PZ= dau1->p4().Pz();
    _dau1_E= dau1->p4().E();
    _dau1_decayMode = (Int_t)(userdatahelpers::getUserFloat(dau1,"decayMode"));
    _dau1_Discriminator= userdatahelpers::getUserFloat(dau1,"HPSDiscriminator");
    _dau1_ISO =userdatahelpers::getUserFloat(dau1,"combRelIsoPF");
    _dau1_dxy =userdatahelpers::getUserFloat(dau1,"dxy");
    _dau1_dz =userdatahelpers::getUserFloat(dau1,"dz");
  }
  myTree->Fill();
  //return;
}

//Fill jets quantities
int HHTauTauNtuplizer::FillJet(const edm::View<pat::Jet> *jets){
  int nJets=0;
  for(edm::View<pat::Jet>::const_iterator ijet = jets->begin(); ijet!=jets->end();++ijet){
    nJets++;
    _jets.push_back(ijet->p4());
    _jetFlavour.push_back(ijet->partonFlavour());
    _bdiscr.push_back(ijet->bDiscriminator("jetBProbabilityBJetTags"));
    _bdiscr2.push_back(ijet->bDiscriminator("combinedSecondaryVertexBJetTags"));
  }
  return nJets;
}


void HHTauTauNtuplizer::FillbQuarks(const edm::Event& event){
  edm::Handle<edm::View<pat::GenericParticle>>candHandle;
  event.getByLabel("bQuarks",candHandle);
  const edm::View<pat::GenericParticle>* bs = candHandle.product();
  for(edm::View<pat::GenericParticle>::const_iterator ib = bs->begin(); ib!=bs->end();++ib){
    const pat::GenericParticle* cand = &(*ib);
    _bquarks.push_back(cand->p4());
    _bmotmass.push_back(userdatahelpers::getUserFloat(cand,"motHmass"));
  }
}

void HHTauTauNtuplizer::endJob(){
  hCounter->SetBinContent(1,Nevt_Gen);
  hCounter->SetBinContent(2,Nevt_HLT);

  hCounter->GetXaxis()->SetBinLabel(1,"Nevt_Gen");
  hCounter->GetXaxis()->SetBinLabel(2,"Nevt_HLT");
}

void HHTauTauNtuplizer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){
  
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
	  //	  edm::LogInfo("AnalyzeRates")<<"Added path "<<pathName<<" to foundPaths";
    } 
  }
  
}

void HHTauTauNtuplizer::endRun(edm::Run const&, edm::EventSetup const&){
}
void HHTauTauNtuplizer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){
}
void HHTauTauNtuplizer::endLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const& iSetup)
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
// void HHTauTauNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
//   //The following says we do not know what parameters are allowed so do no validation
//   // Please change this to state exactly what you do use, even if it is no parameters
//   edm::ParameterSetDescription desc;
//   desc.setUnknown();
//   descriptions.addDefault(desc);
// }

//define this as a plug-in
DEFINE_FWK_MODULE(HHTauTauNtuplizer);
