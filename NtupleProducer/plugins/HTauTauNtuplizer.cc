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
#include "LLRHiggsTauTau/NtupleProducer/Utils/OfflineProducerHelper.h"

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
  int FillJet(const edm::View<pat::Jet>* jet);
  void FillSoftLeptons(const edm::View<reco::Candidate> *dauhandler, bool theFSR);
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
  Int_t Npairs;

  //Event Output variables
  Int_t _indexevents;
  Int_t _runNumber;
  Int_t _triggerbit;
  Float_t _met;
  Float_t _metphi;
  
  //Leptons
  std::vector<math::XYZTLorentzVector> _mothers;
  std::vector<math::XYZTLorentzVector> _daughters;
  std::vector<const reco::Candidate*> _softLeptons;
  std::vector<math::XYZTLorentzVector> _bquarks;
  //std::vector<math::XYZTLorentzVector> _daughter2;

  //Mothers output variables
  std::vector<Int_t> _indexDau1;
  std::vector<Int_t> _indexDau2;
  std::vector<Int_t> _genDaughters;
  std::vector<Float_t> _SVmass;
  std::vector<Float_t> _metx;
  std::vector<Float_t> _mety;
  std::vector<Float_t> _bmotmass;
   
  //Leptons variables
  std::vector<Int_t> _pdgdau;
  std::vector<Int_t> _particleType;//0=muon, 1=e, 2=tau
  std::vector<Float_t> _combreliso;
  std::vector<Float_t> _discriminator;//BDT for ele, discriminator for tau
  std::vector<Float_t> _dxy;
  std::vector<Float_t> _dz;
  std::vector<Int_t> _decayType;//for taus only

  //Jets variables
  Int_t _numberOfJets;
  std::vector<math::XYZTLorentzVector> _jets;
  std::vector<Int_t> _jetFlavour;
  std::vector<Float_t> _bdiscr;
  std::vector<Float_t> _bdiscr2;
};

// ----Constructor and Destructor -----
HTauTauNtuplizer::HTauTauNtuplizer(const edm::ParameterSet& pset) {
  theCandLabel = pset.getUntrackedParameter<string>("CandCollection");
  //theChannel = myHelper.channel();
  theFileName = pset.getUntrackedParameter<string>("fileName");
  skipEmptyEvents = pset.getParameter<bool>("skipEmptyEvents");
  theFSR = pset.getParameter<bool>("applyFSR");
  theisMC = pset.getParameter<bool>("IsMC");
  //writeBestCandOnly = pset.getParameter<bool>("onlyBestCandidate");
  //sampleName = pset.getParameter<string>("sampleName");
  Nevt_Gen=0;
  Npairs=0;
  //triggerResultsLabel = InputTag("TriggerResults","","HLT");
  processName= pset.getParameter<edm::InputTag>("triggerResultsLabel");
  Initialize();
}

HTauTauNtuplizer::~HTauTauNtuplizer(){}
//

void HTauTauNtuplizer::Initialize(){
  _mothers.clear();
  _daughters.clear();
  //_daughter2.clear();
  _softLeptons.clear();
  _genDaughters.clear();
  _bquarks.clear();
  _bmotmass.clear();
  _indexDau1.clear();
  _indexDau2.clear();
  _pdgdau.clear();
  _SVmass.clear();
  _metx.clear();
  _mety.clear();
  _particleType.clear();
  _discriminator.clear();
  _dxy.clear();
  _dz.clear();
  _decayType.clear();
  _combreliso.clear();
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

void HTauTauNtuplizer::beginJob(){
  edm::Service<TFileService> fs;
  myTree = fs->make<TTree>("HTauTauTree","HTauTauTree");
  hCounter = fs->make<TH1F>("Counters","Counters",2,0,2);

  //Branches
  myTree->Branch("EventNumber",&_indexevents,"EventNumber/I");
  myTree->Branch("RunNumber",&_runNumber,"RunNumber/I");
  myTree->Branch("triggerbit",&_triggerbit,"triggerbit/I");
  myTree->Branch("met",&_met,"met/F");
  myTree->Branch("metphi",&_metphi,"metphi/F");  
  myTree->Branch("mothers",&_mothers);
  myTree->Branch("daughters",&_daughters);
  if(writeSoftLep)myTree->Branch("softLeptons",&_softLeptons);
  if(theisMC){
    myTree->Branch("genDaughters",&_genDaughters);
    myTree->Branch("bquarks",&_bquarks);
    myTree->Branch("bmotmass",&_bmotmass);
  }
  //myTree->Branch("daughters2",&_daughter2);
  myTree->Branch("SVfitMass",&_SVmass);
  myTree->Branch("METx",&_metx);
  myTree->Branch("METy",&_mety);
  myTree->Branch("PDGIdDaughters",&_pdgdau);
  myTree->Branch("indexDau1",&_indexDau1);
  myTree->Branch("indexDau2",&_indexDau2);
  myTree->Branch("particleType",&_particleType);
  myTree->Branch("discriminator",&_discriminator);
  myTree->Branch("dxy",&_dxy);
  myTree->Branch("dz",&_dz);
  myTree->Branch("decayMode",&_decayType);
  myTree->Branch("combreliso",& _combreliso);
  myTree->Branch("JetsNumber",&_numberOfJets,"JetsNumber/I");
  myTree->Branch("jets",&_jets);
  myTree->Branch("jetFlavour",&_jetFlavour);
  myTree->Branch("bDiscriminator",&_bdiscr);
  myTree->Branch("bCSVscore",&_bdiscr2);
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
  edm::Handle<edm::View<pat::CompositeCandidate>>candHandle;
  edm::Handle<edm::View<reco::Candidate>>dauHandle;
  edm::Handle<edm::View<pat::Jet>>jetHandle;
  edm::Handle<pat::METCollection> metHandle;
  event.getByLabel(theCandLabel,candHandle);
  event.getByLabel("jets",jetHandle);
  event.getByLabel("softLeptons",dauHandle);
  event.getByLabel("slimmedMETs",metHandle);
  const edm::View<pat::CompositeCandidate>* cands = candHandle.product();
  const edm::View<reco::Candidate>* daus = dauHandle.product();
  const edm::View<pat::Jet>* jets = jetHandle.product();
  const pat::MET &met = metHandle->front();
  //myNtuple->InitializeVariables();
  _indexevents = event.id().event();
  _runNumber = event.id().run();
  triggerhelper myTriggerHelper;
  _triggerbit = myTriggerHelper.FindTriggerBit(event,foundPaths,indexOfPath);
  _met = met.sumEt();
  _metphi = met.phi();
    
  //Do all the stuff here
  //Compute the variables needed for the output and store them in the ntuple
  if(DEBUG)printf("===New Event===\n");

  //Loop over generated b quarks
  if(theisMC)FillbQuarks(event);

  //Loop of softleptons and fill them
  FillSoftLeptons(daus,theFSR);

  //Loop on Jets
  _numberOfJets = 0;
  if(writeJets)_numberOfJets = FillJet(jets);

  //Loop on pairs
  for(edm::View<pat::CompositeCandidate>::const_iterator candi = cands->begin(); candi!=cands->end();++candi){
    Npairs++;
    const pat::CompositeCandidate& cand = (*candi); 
    math::XYZTLorentzVector candp4 = cand.p4();
    _SVmass.push_back(cand.userFloat("SVfitMass"));
    _metx.push_back(cand.userFloat("MEt_px"));
    _mety.push_back(cand.userFloat("MEt_py"));
    
    //if(DEBUG){
      //motherPoint[iMot]=dynamic_cast<const reco::Candidate*>(&*candi);
      //printf("%p %p %p\n",motherPoint[iMot],cand.daughter(0),cand.daughter(1));
    //}
    
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

	if(iCand==0)_indexDau1.push_back(index);
	else _indexDau2.push_back(index);	

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

//Fill jets quantities
int HTauTauNtuplizer::FillJet(const edm::View<pat::Jet> *jets){
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

//Fill all leptons (we keep them all for veto purposes
void HTauTauNtuplizer::FillSoftLeptons(const edm::View<reco::Candidate> *daus, bool theFSR){
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
    //math::XYZTLorentzVector pfour(userdatahelpers::getUserFloat(cand,"genPx"),userdatahelpers::getUserFloat(cand,"genPy"),userdatahelpers::getUserFloat(cand,"genPz"),userdatahelpers::getUserFloat(cand,"genE"));
    if(theisMC)_genDaughters.push_back(userdatahelpers::getUserFloat(cand,"fromH"));

    _softLeptons.push_back(cand);//This is needed also for FindCandIndex
    _pdgdau.push_back(cand->pdgId());
    _combreliso.push_back(userdatahelpers::getUserFloat(cand,"combRelIsoPF"));
    _dxy.push_back(userdatahelpers::getUserFloat(cand,"dxy"));
    _dz.push_back(userdatahelpers::getUserFloat(cand,"dz"));
    int type =OfflineProducerHelper::TAU;
    if(cand->isMuon()) type = OfflineProducerHelper::MUON;
    else if(cand->isElectron()) type = OfflineProducerHelper::ELECTRON;
    _particleType.push_back(type);
    float discr=-1;
    int decay=-1;
    if(type==1)discr=userdatahelpers::getUserFloat(cand,"BDT");
    else if(type==2){
      discr=userdatahelpers::getUserFloat(cand,"HPSDiscriminator");
      decay = userdatahelpers::getUserFloat(cand,"decayMode");
    }
    _discriminator.push_back(discr);
    _decayType.push_back(decay);
  }
}

void HTauTauNtuplizer::FillbQuarks(const edm::Event& event){
  edm::Handle<edm::View<pat::GenericParticle>>candHandle;
  event.getByLabel("bQuarks",candHandle);
  const edm::View<pat::GenericParticle>* bs = candHandle.product();
  for(edm::View<pat::GenericParticle>::const_iterator ib = bs->begin(); ib!=bs->end();++ib){
    const pat::GenericParticle* cand = &(*ib);
    _bquarks.push_back(cand->p4());
    _bmotmass.push_back(userdatahelpers::getUserFloat(cand,"motHmass"));
  }
}

void HTauTauNtuplizer::endJob(){
  hCounter->SetBinContent(1,Nevt_Gen);
  hCounter->SetBinContent(2,Npairs);
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
