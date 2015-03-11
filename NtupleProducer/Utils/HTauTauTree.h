// Author: L. Cadamuro (LLR)
// Date:   09 March 2015
//
// Wrapper for Higgs Tau Tau tree output tree
// Can be either included in a interpreted macro or compiled in c++
// (use `root-config --glibs --cflags`)
//
// Create the HTauTauTree object from the pointer to the tree, then access the stored objects from it
// Common TTree functions GetEntry (entry), GetEntries() are implemented

#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
// #include <Rtypes.h>
//#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

using namespace std;

class HTauTauTree {
public :
   TTree          *_tree;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   Int_t           EventNumber;
   Int_t           RunNumber;
   Int_t           lumi;
   Float_t         MC_weight;
   Int_t           triggerbit;
   Float_t         met;
   Float_t         metphi;
   vector<float>   *mothers_px;
   vector<float>   *mothers_py;
   vector<float>   *mothers_pz;
   vector<float>   *mothers_e;
   vector<float>   *daughters_px;
   vector<float>   *daughters_py;
   vector<float>   *daughters_pz;
   vector<float>   *daughters_e;
   vector<int>     *genDaughters;
   vector<float>   *bquarks_px;
   vector<float>   *bquarks_py;
   vector<float>   *bquarks_pz;
   vector<float>   *bquarks_e;
   vector<float>   *bmotmass;
   vector<float>   *SVfitMass;
   vector<float>   *METx;
   vector<float>   *METy;
   vector<int>     *PDGIdDaughters;
   vector<int>     *indexDau1;
   vector<int>     *indexDau2;
   vector<int>     *particleType;
   vector<float>   *discriminator;
   vector<float>   *dxy;
   vector<float>   *dz;
   vector<int>     *decayMode;
   vector<float>   *combreliso;
   Int_t           JetsNumber;
   vector<float>   *jets_px;
   vector<float>   *jets_py;
   vector<float>   *jets_pz;
   vector<float>   *jets_e;
   vector<int>     *jets_Flavour;
   vector<float>   *bDiscriminator;
   vector<float>   *bCSVscore;

   // List of branches
   TBranch        *b_EventNumber;   //!
   TBranch        *b_RunNumber;   //!
   TBranch        *b_lumi           //!
   TBranch        *b_MC_weight           //!
   TBranch        *b_triggerbit;   //!
   TBranch        *b_met;   //!
   TBranch        *b_metphi;   //!
   TBranch        *b_mothers_px;   //!
   TBranch        *b_mothers_py;   //!
   TBranch        *b_mothers_pz;   //!
   TBranch        *b_mothers_e;   //!
   TBranch        *b_daughters_px;   //!
   TBranch        *b_daughters_py;   //!
   TBranch        *b_daughters_pz;   //!
   TBranch        *b_daughters_e;   //!
   TBranch        *b_genDaughters;   //!
   TBranch        *b_bquarks_px;   //!
   TBranch        *b_bquarks_py;   //!
   TBranch        *b_bquarks_pz;   //!
   TBranch        *b_bquarks_e;   //!
   TBranch        *b_bmotmass;   //!
   TBranch        *b_SVfitMass;   //!
   TBranch        *b_METx;   //!
   TBranch        *b_METy;   //!
   TBranch        *b_PDGIdDaughters;   //!
   TBranch        *b_indexDau1;   //!
   TBranch        *b_indexDau2;   //!
   TBranch        *b_particleType;   //!
   TBranch        *b_discriminator;   //!
   TBranch        *b_dxy;   //!
   TBranch        *b_dz;   //!
   TBranch        *b_decayMode;   //!
   TBranch        *b_combreliso;   //!
   TBranch        *b_JetsNumber;   //!
   TBranch        *b_jets_px;   //!
   TBranch        *b_jets_py;   //!
   TBranch        *b_jets_pz;   //!
   TBranch        *b_jets_e;   //!
   TBranch        *b_jets_Flavour;   //!
   TBranch        *b_bDiscriminator;   //!
   TBranch        *b_bCSVscore;   //!
   
   // methods
   HTauTauTree (TTree* tree); //ctor
   ~HTauTauTree();
   void Init(TTree* tree);
   Int_t GetEntry(int entry);
   Long64_t GetEntries();
   TTree* GetTree();
};

HTauTauTree::HTauTauTree (TTree* tree)
{
    Init(tree);
}

HTauTauTree::~HTauTauTree() {}

void HTauTauTree::Init(TTree* tree)
{
   // Set object pointer
   mothers_px = 0;
   mothers_py = 0;
   mothers_pz = 0;
   mothers_e = 0;
   daughters_px = 0;
   daughters_py = 0;
   daughters_pz = 0;
   daughters_e = 0;
   genDaughters = 0;
   bquarks_px = 0;
   bquarks_py = 0;
   bquarks_pz = 0;
   bquarks_e = 0;
   bmotmass = 0;
   SVfitMass = 0;
   METx = 0;
   METy = 0;
   PDGIdDaughters = 0;
   indexDau1 = 0;
   indexDau2 = 0;
   particleType = 0;
   discriminator = 0;
   dxy = 0;
   dz = 0;
   decayMode = 0;
   combreliso = 0;
   jets_px = 0;
   jets_py = 0;
   jets_pz = 0;
   jets_e = 0;
   jets_Flavour = 0;
   bDiscriminator = 0;
   bCSVscore = 0;
   
   // Set branch addresses and branch pointers
   if (!tree) return;
   _tree = tree;
   
   _tree->SetMakeClass(1); // needed especially when compiling

   _tree->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   _tree->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
   _tree->SetBranchAddress("lumi", &lumi, &b_lumi);
   _tree->SetBranchAddress("triggerbit", &triggerbit, &b_triggerbit);
   _tree->SetBranchAddress("met", &met, &b_met);
   _tree->SetBranchAddress("metphi", &metphi, &b_metphi);
   _tree->SetBranchAddress("mothers_px", &mothers_px, &b_mothers_px);
   _tree->SetBranchAddress("mothers_py", &mothers_py, &b_mothers_py);
   _tree->SetBranchAddress("mothers_pz", &mothers_pz, &b_mothers_pz);
   _tree->SetBranchAddress("mothers_e", &mothers_e, &b_mothers_e);
   _tree->SetBranchAddress("daughters_px", &daughters_px, &b_daughters_px);
   _tree->SetBranchAddress("daughters_py", &daughters_py, &b_daughters_py);
   _tree->SetBranchAddress("daughters_pz", &daughters_pz, &b_daughters_pz);
   _tree->SetBranchAddress("daughters_e", &daughters_e, &b_daughters_e);
   _tree->SetBranchAddress("SVfitMass", &SVfitMass, &b_SVfitMass);
   _tree->SetBranchAddress("METx", &METx, &b_METx);
   _tree->SetBranchAddress("METy", &METy, &b_METy);
   _tree->SetBranchAddress("PDGIdDaughters", &PDGIdDaughters, &b_PDGIdDaughters);
   _tree->SetBranchAddress("indexDau1", &indexDau1, &b_indexDau1);
   _tree->SetBranchAddress("indexDau2", &indexDau2, &b_indexDau2);
   _tree->SetBranchAddress("particleType", &particleType, &b_particleType);
   _tree->SetBranchAddress("discriminator", &discriminator, &b_discriminator);
   _tree->SetBranchAddress("dxy", &dxy, &b_dxy);
   _tree->SetBranchAddress("dz", &dz, &b_dz);
   _tree->SetBranchAddress("decayMode", &decayMode, &b_decayMode);
   _tree->SetBranchAddress("combreliso", &combreliso, &b_combreliso);
   _tree->SetBranchAddress("JetsNumber", &JetsNumber, &b_JetsNumber);
   _tree->SetBranchAddress("jets_px", &jets_px, &b_jets_px);
   _tree->SetBranchAddress("jets_py", &jets_py, &b_jets_py);
   _tree->SetBranchAddress("jets_pz", &jets_pz, &b_jets_pz);
   _tree->SetBranchAddress("jets_e", &jets_e, &b_jets_e);
   _tree->SetBranchAddress("jets_Flavour", &jets_Flavour, &b_jets_Flavour);
   _tree->SetBranchAddress("bDiscriminator", &bDiscriminator, &b_bDiscriminator);
   _tree->SetBranchAddress("bCSVscore", &bCSVscore, &b_bCSVscore);

   // MC only
   if(_tree->GetListOfBranches()->FindObject("genDaughters"))
   {
        _tree->SetBranchAddress("genDaughters", &genDaughters, &b_genDaughters);
        _tree->SetBranchAddress("bquarks_px", &bquarks_px, &b_bquarks_px);
        _tree->SetBranchAddress("bquarks_py", &bquarks_py, &b_bquarks_py);
        _tree->SetBranchAddress("bquarks_pz", &bquarks_pz, &b_bquarks_pz);
        _tree->SetBranchAddress("bquarks_e", &bquarks_e, &b_bquarks_e);
        _tree->SetBranchAddress("bmotmass", &bmotmass, &b_bmotmass);
        _tree->SetBranchAddress("MC_weight".&MC_weight,&b_MC_weight);
   }
}

Int_t HTauTauTree::GetEntry(int entry)
{
    return _tree->GetEntry(entry);
} 

Long64_t HTauTauTree::GetEntries()
{
    return _tree->GetEntries();
}

TTree* HTauTauTree::GetTree()
{
    return _tree;
}
