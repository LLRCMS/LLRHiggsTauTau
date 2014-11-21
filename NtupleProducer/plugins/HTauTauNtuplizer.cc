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

#include "TLorentzVector.h"

// namespace {
//   bool writePhotons = false;  // Write photons in the tree. 
//   bool writeJets = true;     // Write jets in the tree. 
//   bool writeVetoLep = true;
// }
namespace{
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
  ~HTauTauNtuplizer();  

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
  //virtual void FillJet(const pat::PFJet& jet);

  // ----------member data ---------------------------
  //HTauTauConfigHelper myHelper;
  int theChannel;
  std::string theCandLabel;
  TString theFileName;
  bool theFSR;

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
  std::vector<math::XYZTLorentzVector> _daughters;
  std::vector<const reco::Candidate*> _softLeptons;
  //std::vector<math::XYZTLorentzVector> _daughter2;
  std::vector<Int_t> _indexDau1;
  std::vector<Int_t> _indexDau2;
  std::vector<Int_t> _pdgmot;
  std::vector<Int_t> _pdgdau;
  Int_t _indexevents;
  Int_t _runNumber;
};

// ----Constructor and Desctructor -----
HTauTauNtuplizer::HTauTauNtuplizer(const edm::ParameterSet& pset) {
  theCandLabel = pset.getUntrackedParameter<string>("CandCollection");
  //theChannel = myHelper.channel();
  theFileName = pset.getUntrackedParameter<string>("fileName");
  skipEmptyEvents = pset.getParameter<bool>("skipEmptyEvents");
  theFSR = pset.getParameter<bool>("applyFSR");
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
  _daughters.clear();
  //_daughter2.clear();
  _softLeptons.clear();
  _indexDau1.clear();
  _indexDau2.clear();
  _pdgmot.clear();
  _pdgdau.clear();
  _indexevents=0;
  _runNumber=0;
}

void HTauTauNtuplizer::beginJob(){
  edm::Service<TFileService> fs;
  myTree = fs->make<TTree>("HTauTauTree","HTauTauTree");
  hCounter = fs->make<TH1F>("Counters","Counters",2,0,2);

  //Branches
  myTree->Branch("EventNumber",&_indexevents,"EventNumber/I");
  myTree->Branch("RunNumber",&_runNumber,"RunNumber/I");
  myTree->Branch("mothers",&_mothers);
  myTree->Branch("daughters",&_daughters);
  //myTree->Branch("daughters2",&_daughter2);
  myTree->Branch("PDGIdMothers",&_pdgmot);
  myTree->Branch("PDGIdDaughters",&_pdgdau);
  myTree->Branch("indexDau1",&_indexDau1);
  myTree->Branch("indexDau2",&_indexDau2);
  //myTree->Branch("softLeptons",&_softLeptons);
}

Int_t HTauTauNtuplizer::FindCandIndex(const reco::Candidate& cand,Int_t iCand=0){
  const reco::Candidate *daughter = cand.daughter(iCand);
  for(UInt_t iLeptons=0;iLeptons<_softLeptons.size();iLeptons++){
	//if(daughter==daughterPoint[iLeptons]){
    if(daughter==_softLeptons.at(iLeptons)){
      return iLeptons;
    }
  }
  return -1;
}
// ----Analyzer (main) ----
// ------------ method called for each event  ------------
void HTauTauNtuplizer::analyze(const edm::Event& event, const edm::EventSetup& eSetup)
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


  //Get candidate collection
  //edm::Handle<edm::View<pat::CompositeCandidate>>candHandle;
  //event.getByLabel("barellCand",candHandle);
  //const edm::View<pat::CompositeCandidate>* cands = candHandle.product();
  edm::Handle<edm::View<reco::CompositeCandidate>>candHandle;
  edm::Handle<edm::View<reco::Candidate>>dauHandle;
  event.getByLabel("barellCand",candHandle);
  event.getByLabel("softLeptons",dauHandle);
  const edm::View<reco::CompositeCandidate>* cands = candHandle.product();
  const edm::View<reco::Candidate>* daus = dauHandle.product();

  //myNtuple->InitializeVariables();
  _indexevents = event.id().event();
  _runNumber = event.id().run();
  //Int_t nCands = daus->size()*2;
  //const reco::Candidate *daughterPoint[nCands];

  //Do all the stuff here
  //Compute the variables needed for the output and store them in the ntuple
  int nDaughters=0;
  if(DEBUG)printf("===New Event===\n");
  //Loop of softleptons and fill them
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
    _daughters.push_back(pfour);
    _softLeptons.push_back(cand);
    _pdgdau.push_back(cand->pdgId());
  }

  //Loop on pairs
  for(edm::View<reco::CompositeCandidate>::const_iterator candi = cands->begin(); candi!=cands->end();++candi){
    Npairs++;
    const reco::CompositeCandidate& cand = (*candi); 
    math::XYZTLorentzVector candp4 = cand.p4();
    _pdgmot.push_back(cand.pdgId());

    //if(DEBUG){
      //motherPoint[iMot]=dynamic_cast<const reco::Candidate*>(&*candi);
      //printf("%p %p %p\n",motherPoint[iMot],cand.daughter(0),cand.daughter(1));
    //}

    //We need to find a way to avoid repetitions of daughter in order to save space
    //It would be nice to move this into a separate function (FindCandIndex)
    for(int iCand=0;iCand<2;iCand++){
      //int index=-1;
      int index=FindCandIndex(cand,iCand);
      const reco::Candidate *daughter = cand.daughter(iCand);
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
	      maxPT  = pt;
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

      if(index>=0){//Daughter already in the list!
	if(iCand==0)_indexDau1.push_back(index);
	else _indexDau2.push_back(index);	
      }else {//new lepton
	//daughterPoint[nDaughters]=daughter;
	nDaughters++;
	//_daughters.push_back(pfour);
	//_pdgdau.push_back(daughter->pdgId());
	if(iCand==0)_indexDau1.push_back(_daughters.size()-1);
	else _indexDau2.push_back(_daughters.size()-1);
      }
    }
    _mothers.push_back(candp4);
  
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
