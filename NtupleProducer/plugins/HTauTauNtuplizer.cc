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
#include "LLRHiggsTauTau/NtupleProducer/Utils/OfflineProducerHelper.h"

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
  
  //Leptons
  //std::vector<TLorentzVector> _mothers;
  std::vector<Float_t> _mothers_px;
  std::vector<Float_t> _mothers_py;
  std::vector<Float_t> _mothers_pz;
  std::vector<Float_t> _mothers_e;
  
  //std::vector<TLorentzVector> _daughters;
  std::vector<Float_t> _daughters_px;
  std::vector<Float_t> _daughters_py;
  std::vector<Float_t> _daughters_pz;
  std::vector<Float_t> _daughters_e;

  std::vector<const reco::Candidate*> _softLeptons;
  
  //std::vector<TLorentzVector> _bquarks;
  std::vector<Float_t> _bquarks_px;
  std::vector<Float_t> _bquarks_py;
  std::vector<Float_t> _bquarks_pz;
  std::vector<Float_t> _bquarks_e;
  std::vector<Int_t> _bquarks_pdg;
  
  //std::vector<math::XYZTLorentzVector> _daughter2;

  //Mothers output variables
  std::vector<Int_t> _indexDau1;
  std::vector<Int_t> _indexDau2;
  std::vector<Int_t> _genDaughters;
  std::vector<Float_t> _SVmass;
  std::vector<Float_t> _metx;
  std::vector<Float_t> _mety;
  std::vector<Float_t> _metCov00;
  std::vector<Float_t> _metCov01;
  std::vector<Float_t> _metCov10;
  std::vector<Float_t> _metCov11;
  std::vector<Float_t> _bmotmass;
   
  //Leptons variables
  std::vector<Int_t> _pdgdau;
  std::vector<Int_t> _particleType;//0=muon, 1=e, 2=tau
  std::vector<Float_t> _combreliso;
  std::vector<Float_t> _discriminator;//BDT for ele, discriminator for tau, muonID for muons (bit 0 loose, 1 soft 2 tight)
  std::vector<Float_t> _dxy;
  std::vector<Float_t> _dz;
  std::vector<Int_t> _decayType;//for taus only
  std::vector<Float_t> _daughters_IetaIeta;
  std::vector<Float_t> _daughters_deltaPhiSuperClusterTrackAtVtx;
  std::vector<Float_t> _daughters_depositR03_tracker;
  std::vector<Float_t> _daughters_depositR03_ecal;
  std::vector<Float_t> _daughters_depositR03_hcal;

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
  
  //genH
  std::vector<Float_t> _genH_px;
  std::vector<Float_t> _genH_py;
  std::vector<Float_t> _genH_pz;
  std::vector<Float_t> _genH_e;
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
  Nevt_PassTrigger = 0;
  Npairs=0;
  //triggerResultsLabel = InputTag("TriggerResults","","HLT");
  processName= pset.getParameter<edm::InputTag>("triggerResultsLabel");
  Initialize();
}

HTauTauNtuplizer::~HTauTauNtuplizer(){}
//

void HTauTauNtuplizer::Initialize(){
  //_mothers.clear();
  _mothers_px.clear();
  _mothers_px.clear();
  _mothers_pz.clear();
  _mothers_e.clear();
  
  //_daughters.clear();
  _daughters_px.clear();
  _daughters_py.clear();
  _daughters_pz.clear();
  _daughters_e.clear();
  _daughters_IetaIeta.clear();
  _daughters_deltaPhiSuperClusterTrackAtVtx.clear();
  _daughters_depositR03_tracker.clear();  
  _daughters_depositR03_ecal.clear();  
  _daughters_depositR03_hcal.clear();  

  //_daughter2.clear();
  _softLeptons.clear();
  _genDaughters.clear();
  
  //_bquarks.clear();
  _bquarks_px.clear();
  _bquarks_py.clear();
  _bquarks_pz.clear();
  _bquarks_e.clear();
  _bquarks_pdg.clear();
  _bmotmass.clear();
  
  _indexDau1.clear();
  _indexDau2.clear();
  _pdgdau.clear();
  _SVmass.clear();
  _metx.clear();
  _mety.clear();
  _metCov00.clear();
  _metCov01.clear();
  _metCov10.clear();
  _metCov11.clear();
  _particleType.clear();
  _discriminator.clear();
  _dxy.clear();
  _dz.clear();
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
  
  _genH_px.clear();
  _genH_py.clear();
  _genH_pz.clear();
  _genH_e.clear();
}

void HTauTauNtuplizer::beginJob(){
  edm::Service<TFileService> fs;
  myTree = fs->make<TTree>("HTauTauTree","HTauTauTree");
  hCounter = fs->make<TH1F>("Counters","Counters",3,0,3);

  //Branches
  myTree->Branch("EventNumber",&_indexevents,"EventNumber/I");
  myTree->Branch("RunNumber",&_runNumber,"RunNumber/I");
  myTree->Branch("lumi",&_lumi,"lumi/I");
  myTree->Branch("triggerbit",&_triggerbit,"triggerbit/I");
  myTree->Branch("metfilterbit",&_metfilterbit,"metfilterbit/I");
  myTree->Branch("met",&_met,"met/F");
  myTree->Branch("metphi",&_metphi,"metphi/F");  
  
  myTree->Branch("mothers_px",&_mothers_px);
  myTree->Branch("mothers_py",&_mothers_py);
  myTree->Branch("mothers_pz",&_mothers_pz);
  myTree->Branch("mothers_e",&_mothers_e);

  myTree->Branch("daughters_px",&_daughters_px);
  myTree->Branch("daughters_py",&_daughters_py);
  myTree->Branch("daughters_pz",&_daughters_pz);
  myTree->Branch("daughters_e",&_daughters_e);

  if(writeSoftLep)myTree->Branch("softLeptons",&_softLeptons);
  if(theisMC){
    myTree->Branch("genDaughters",&_genDaughters);
    myTree->Branch("quarks_px",&_bquarks_px);
    myTree->Branch("quarks_py",&_bquarks_py);
    myTree->Branch("quarks_pz",&_bquarks_pz);
    myTree->Branch("quarks_e",&_bquarks_e);
    myTree->Branch("quarks_pdg",&_bquarks_e);
    myTree->Branch("motmass",&_bmotmass);
    myTree->Branch("MC_weight",&_MC_weight,"MC_weight/F");
    myTree->Branch("genH_px",&_genH_px);
    myTree->Branch("genH_py",&_genH_py);
    myTree->Branch("genH_pz",&_genH_pz);
    myTree->Branch("genH_e",&_genH_e);
  }
  //myTree->Branch("daughters2",&_daughter2);
  myTree->Branch("SVfitMass",&_SVmass);
  myTree->Branch("METx",&_metx);
  myTree->Branch("METy",&_mety);
  myTree->Branch("MET_cov00",&_metCov00);
  myTree->Branch("MET_cov01",&_metCov01);
  myTree->Branch("MET_cov10",&_metCov10);
  myTree->Branch("MET_cov11",&_metCov11);  
  myTree->Branch("PDGIdDaughters",&_pdgdau);
  myTree->Branch("indexDau1",&_indexDau1);
  myTree->Branch("indexDau2",&_indexDau2);
  myTree->Branch("particleType",&_particleType);
  myTree->Branch("discriminator",&_discriminator);
  myTree->Branch("dxy",&_dxy);
  myTree->Branch("dz",&_dz);
  myTree->Branch("decayMode",&_decayType);
  myTree->Branch("combreliso",& _combreliso);
  myTree->Branch("daughters_IetaIeta",&_daughters_IetaIeta);
  myTree->Branch("daughters_deltaPhiSuperClusterTrackAtVtx",&_daughters_deltaPhiSuperClusterTrackAtVtx);
  myTree->Branch("daughters_depositR03_tracker",&_daughters_depositR03_tracker);
  myTree->Branch("daughters_depositR03_ecal",&_daughters_depositR03_ecal);
  myTree->Branch("daughters_depositR03_hcal",&_daughters_depositR03_hcal);
  myTree->Branch("JetsNumber",&_numberOfJets,"JetsNumber/I");
  myTree->Branch("jets_px",&_jets_px);
  myTree->Branch("jets_py",&_jets_py);
  myTree->Branch("jets_pz",&_jets_pz);
  myTree->Branch("jets_e",&_jets_e);
  myTree->Branch("jets_Flavour",&_jets_Flavour);
  myTree->Branch("jets_PUJetID",&_jets_PUJetID);
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
  _metfilterbit = myTriggerHelper.FindMETBit(event);

  //Get candidate collection
  edm::Handle<edm::View<pat::CompositeCandidate>>candHandle;
  edm::Handle<edm::View<reco::Candidate>>dauHandle;
  edm::Handle<edm::View<pat::Jet>>jetHandle;
  edm::Handle<pat::METCollection> metHandle;
  edm::Handle<GenFilterInfo> embeddingWeightHandle;
  edm::Handle<edm::TriggerResults> triggerResults;
  
  event.getByLabel(theCandLabel,candHandle);
  event.getByLabel("jets",jetHandle);
  event.getByLabel("softLeptons",dauHandle);
  event.getByLabel("slimmedMETs",metHandle);
  if(theisMC){
    event.getByLabel(edm::InputTag("generator","minVisPtFilter",""),embeddingWeightHandle);
    _MC_weight = (embeddingWeightHandle.isValid() ? embeddingWeightHandle->filterEfficiency():1.0);
  }

  const edm::View<pat::CompositeCandidate>* cands = candHandle.product();
  const edm::View<reco::Candidate>* daus = dauHandle.product();
  const edm::View<pat::Jet>* jets = jetHandle.product();
  const pat::MET &met = metHandle->front();
  //myNtuple->InitializeVariables();
    
  _indexevents = event.id().event();
  _runNumber = event.id().run();
  _lumi=event.luminosityBlock();
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
    _metCov00.push_back(cand.userFloat("MEt_cov00"));
    _metCov01.push_back(cand.userFloat("MEt_cov01"));
    _metCov10.push_back(cand.userFloat("MEt_cov10"));
    _metCov11.push_back(cand.userFloat("MEt_cov11"));
    
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

        if(iCand==0)_indexDau1.push_back(index);
        else _indexDau2.push_back(index);	

    }
    
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
    _bdiscr2.push_back(ijet->bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags"));
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
    
    _daughters_px.push_back( (float) pfour.X());
    _daughters_py.push_back( (float) pfour.Y());
    _daughters_pz.push_back( (float) pfour.Z());
    _daughters_e.push_back( (float) pfour.T());

    //math::XYZTLorentzVector pfour(userdatahelpers::getUserFloat(cand,"genPx"),userdatahelpers::getUserFloat(cand,"genPy"),userdatahelpers::getUserFloat(cand,"genPz"),userdatahelpers::getUserFloat(cand,"genE"));
    if(theisMC)_genDaughters.push_back(userdatahelpers::getUserFloat(cand,"fromH"));
    _softLeptons.push_back(cand);//This is needed also for FindCandIndex
    _pdgdau.push_back(cand->pdgId());
    _combreliso.push_back(userdatahelpers::getUserFloat(cand,"combRelIsoPF"));
    _dxy.push_back(userdatahelpers::getUserFloat(cand,"dxy"));
    _dz.push_back(userdatahelpers::getUserFloat(cand,"dz"));
    int type = OfflineProducerHelper::TAU;
    if(cand->isMuon()) type = OfflineProducerHelper::MUON;
    else if(cand->isElectron()) type = OfflineProducerHelper::ELECTRON;
    _particleType.push_back(type);
    float discr=-1;
    int decay=-1;
    float ieta=-1,superatvtx=-1,depositTracker=-1,depositEcal=-1,depositHcal=-1;
    if(type==0){
      discr=userdatahelpers::getUserFloat(cand,"muonID");
      depositTracker=userdatahelpers::getUserFloat(cand,"DepositR03TrackerOfficial");
      depositEcal=userdatahelpers::getUserFloat(cand,"DepositR03Ecal");
      depositHcal=userdatahelpers::getUserFloat(cand,"DepositR03Hcal");
    }else if(type==1){
      discr=userdatahelpers::getUserFloat(cand,"BDT");
      ieta=userdatahelpers::getUserFloat(cand,"sigmaIetaIeta");
      superatvtx=userdatahelpers::getUserFloat(cand,"deltaPhiSuperClusterTrackAtVtx");
    }else if(type==2){
      discr=userdatahelpers::getUserFloat(cand,"HPSDiscriminator");
      decay = userdatahelpers::getUserFloat(cand,"decayMode");
    }
    _discriminator.push_back(discr);
    _decayType.push_back(decay);
    _daughters_IetaIeta.push_back(ieta);
    _daughters_deltaPhiSuperClusterTrackAtVtx.push_back(superatvtx);
    _daughters_depositR03_tracker.push_back(depositTracker);
    _daughters_depositR03_ecal.push_back(depositEcal);
    _daughters_depositR03_hcal.push_back(depositHcal);
  }
}

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
	    _genH_px.push_back(packed->px());          
	    _genH_py.push_back(packed->py());          
	    _genH_pz.push_back(packed->pz());          
	    _genH_e.push_back(packed->energy());          
	  }        
	} 

}

void HTauTauNtuplizer::endJob(){
  hCounter->SetBinContent(1,Nevt_Gen);
  hCounter->SetBinContent(2,Nevt_PassTrigger);
  hCounter->SetBinContent(3,Npairs);

  hCounter->GetXaxis()->SetBinLabel(1,"Nevt_Gen");
  hCounter->GetXaxis()->SetBinLabel(2,"Nevt_PassTrigger");
  hCounter->GetXaxis()->SetBinLabel(3,"Npairs");
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
  // Total number of events is the sum of the events in each of these luminosity blocks
  edm::Handle<edm::MergeableCounter> nEventsTotalCounter;
  iLumi.getByLabel("nEventsTotal", nEventsTotalCounter);
  Nevt_Gen += nEventsTotalCounter->value;

  edm::Handle<edm::MergeableCounter> nEventsPassTrigCounter;
  iLumi.getByLabel("nEventsPassTrigger", nEventsPassTrigCounter);
  Nevt_PassTrigger += nEventsPassTrigCounter->value;
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
