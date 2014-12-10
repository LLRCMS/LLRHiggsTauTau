//#include <DataFormats/Candidate/interface/Candidate.h>
#include "TMath.h"
#include <memory>
#include <cmath>
#include <vector>
#include <string>
#include <TNtuple.h>
#include "TLorentzVector.h"
//#include <DataFormats/Math/interface/LorentzVectorFwd.h>

//using namespace reco;
using namespace std;

void ReduceTree(){
  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();
  gSystem->Load("libDataFormatsFWLite.so");
  gSystem->Load("libDataFormatsPatCandidates.so");

  TFile *f = TFile::Open("../test/HTauTauTree.root");
  TTree *tree = ((TDirectory*)f->Get("HTauTauTree"))->Get("HTauTauTree");
 
TFile *fnew = new TFile("HTauTauTree_noReco.root","RECREATE");
 TDirectoryFile* dir = fnew->mkdir("HTauTauTree");
 TTree *myTree = new TTree("HTauTauTree","SlimmedTree");
  /*
  TBranch* b= tree->GetBranch("softLeptons");
  tree->GetListOfBranches()->Remove(b);
  TLeaf* l= tree->GetLeaf("softLeptons");
  tree->GetListOfLeaves()->Remove(l);
  */
  
  Int_t _indexevents;
  Int_t _runNumber;
  Int_t _triggerbit;
  Float_t _met;
  Float_t _metphi;
  
  //Leptons
  vector<math::XYZTLorentzVector> *_mothers;//fsr corrected
  vector<math::XYZTLorentzVector> *_daughters;//fsr corrected
  vector<const reco::Candidate*> *_softLeptons;//should we keep?
  //vector<math::XYZTLorentzVector> _daughter2;

  //Mothers output variables
  vector<Int_t> *_indexDau1;
  vector<Int_t> *_indexDau2;
  vector<Float_t> *_SVmass;
  vector<Float_t> *_metx;
  vector<Float_t> *_mety;
   
  //Leptons variables
  vector<Int_t> *_pdgdau;
  vector<Int_t> *_particleType;//0=muon, 1=e, 2=tau
  vector<Float_t> *_combreliso;
  vector<Float_t> *_discriminator;//BDT for ele, discriminator for tau
  vector<Float_t> *_dxy;
  vector<Float_t> *_dz;
  vector<Int_t> *_decayType;//for taus only

  //Jets variables
  Int_t _numberOfJets;
  vector<math::XYZTLorentzVector> *_jets;
  vector<Float_t> *_bdiscr;
  vector<Float_t> *_bdiscr2;

 
  Int_t myindexevents;
  Int_t myrunNumber;
  Int_t mytriggerbit;
  Float_t mymet;
  Float_t mymetphi;
  
  //Leptons
  vector<math::XYZTLorentzVector> mymothers;//fsr corrected
  vector<math::XYZTLorentzVector> mydaughters;//fsr corrected
  //Mothers output variables
  vector<Int_t> myindexDau1;
  vector<Int_t> myindexDau2;
  vector<Float_t> mySVmass;
  vector<Float_t> mymetx;
  vector<Float_t> mymety;
   
  //Leptons variables
  vector<Int_t> mypdgdau;
  vector<Int_t> myparticleType;//0=muon, 1=e, 2=tau
  vector<Float_t> mycombreliso;
  vector<Float_t> mydiscriminator;//BDT for ele, discriminator for tau
  vector<Float_t> mydxy;
  vector<Float_t> mydz;
  vector<Int_t> mydecayType;//for taus only

  //Jets variables
  Int_t mynumberOfJets;
  vector<math::XYZTLorentzVector> myjets;
  vector<Float_t> mybdiscr;
  vector<Float_t> mybdiscr2;


  TBranch * b_EventNumber;
 TBranch * b_runNumber;
 TBranch * b_triggerbit;
 TBranch * b_met;
 TBranch * b_metphi;
  
  //Leptons
   TBranch * b_mothers;//fsr corrected
   TBranch * b_daughters;//fsr corrected
   TBranch * b_softLeptons;//should we keep?
  //vector<math::XYZTLorentzVector> _daughter2;

  //Mothers output variables
   TBranch * b_indexDau1;
   TBranch * b_indexDau2;
   TBranch * b_SVmass;
   TBranch * b_metx;
   TBranch * b_mety;
  
  
   TBranch * b_pdgdau;
   TBranch * b_particleType;//0=muon, 1=e, 2=tau
   TBranch * b_combreliso;
   TBranch * b_discriminator;//BDT for ele, discriminator for tau
   TBranch * b_dxy;
   TBranch * b_dz;
   TBranch * b_decayType;//for taus only

  //Jets variables
   TBranch * b_numberOfJets;
   TBranch * b_jets;
   TBranch * b_bdiscr;
   TBranch * b_bdiscr2;


   tree->SetBranchAddress("EventNumber",&_indexevents,&b_EventNumber);  
   tree->SetBranchAddress("RunNumber",&_runNumber,&b_runNumber);  
   tree->SetBranchAddress("triggerbit",&_triggerbit,&b_triggerbit);  
   tree->SetBranchAddress("met",&_met,&b_met);  
   tree->SetBranchAddress("metphi",&_metphi,&b_metphi);    
   tree->SetBranchAddress("mothers",&_mothers,&b_mothers);  
   tree->SetBranchAddress("daughters",&_daughters,&b_daughters);  
   tree->SetBranchAddress("softLeptons",&_softLeptons,&b_softLeptons);  
   tree->SetBranchAddress("SVfitMass",&_SVmass,&b_SVmass);  
   tree->SetBranchAddress("METx",&_metx,&b_metx);  
   tree->SetBranchAddress("METy",&_mety,&b_mety);  
   tree->SetBranchAddress("PDGIdDaughters",&_pdgdau,&b_pdgdau);  
   tree->SetBranchAddress("indexDau1",&_indexDau1,&b_indexDau1);  
   tree->SetBranchAddress("indexDau2",&_indexDau2,&b_indexDau2);  
   tree->SetBranchAddress("particleType",&_particleType,&b_particleType);  
   tree->SetBranchAddress("discriminator",&_discriminator,&b_discriminator);  
   tree->SetBranchAddress("dxy",&_dxy,&b_dxy);  
   tree->SetBranchAddress("dz",&_dz,&b_dz);  
   tree->SetBranchAddress("decayMode",&_decayType,&b_decayType);  
   tree->SetBranchAddress("combreliso",& _combreliso,&b_combreliso);  
   tree->SetBranchAddress("JetsNumber",&_numberOfJets,&b_numberOfJets);  
   tree->SetBranchAddress("jets",&_jets,&b_jets);  
   tree->SetBranchAddress("bDiscriminator",&_bdiscr,&b_bdiscr);  
   tree->SetBranchAddress("bCSVscore",&_bdiscr2,&b_bdiscr2);


   myTree->Branch("EventNumber",&myindexevents,"EventNumber/I");
   myTree->Branch("RunNumber",&myrunNumber,"RunNumber/I");
   myTree->Branch("triggerbit",&mytriggerbit,"triggerbit/I");
   myTree->Branch("met",&mymet,"met/F");
   myTree->Branch("metphi",&mymetphi,"metphi/F");  
   myTree->Branch("mothers",&mymothers);
   myTree->Branch("daughters",&mydaughters);
   
   myTree->Branch("SVfitMass",&mySVmass);
   myTree->Branch("METx",&mymetx);
   myTree->Branch("METy",&mymety);
   myTree->Branch("PDGIdDaughters",&mypdgdau);
   myTree->Branch("indexDau1",&myindexDau1);
   myTree->Branch("indexDau2",&myindexDau2);
   myTree->Branch("particleType",&myparticleType);
   myTree->Branch("discriminator",&mydiscriminator);
   myTree->Branch("dxy",&mydxy);
   myTree->Branch("dz",&mydz);
   myTree->Branch("decayMode",&mydecayType);
   myTree->Branch("combreliso",& mycombreliso);
   myTree->Branch("JetsNumber",&mynumberOfJets,"JetsNumber/I");
   myTree->Branch("jets",&myjets);
   myTree->Branch("bDiscriminator",&mybdiscr);
   myTree->Branch("bCSVscore",&mybdiscr2);
   
   int nentries = tree->GetEntriesFast();
   for(int i=0;i<nentries;i++){
     int entry = tree->GetEntry(i);


     myindexevents=_indexevents;
     myrunNumber = _runNumber;
     mytriggerbit=_triggerbit;
     mymet=_met;
     mymetphi=_metphi;
     mymothers=_mothers;
     mydaughters=_daughters;
     mySVmass=_SVmass;
     mymetx=_metx;
     mymety=_mety;
     mypdgdau=_pdgdau;
     myindexDau1=_indexDau1;
     myindexDau2=_indexDau2;
     myparticleType=_particleType;
     mydiscriminator=_discriminator;
     mydxy=_dxy;
     mydz=_dz;
     mydecayType=_decayType;
     mycombreliso=_combreliso;
     mynumberOfJets=_numberOfJets;
     myjets=_jets;
     mybdiscr=_bdiscr;
     mybdiscr2=_bdiscr2;

     myTree->Fill();
   }
 
 dir->cd();  

 myTree->Write();

}
