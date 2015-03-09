#include "HTauTauTree.h"
#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"

using namespace std;

void PlotTreeVars()
{
    TFile* fIn = new TFile ("/home/llr/cms/cadamuro/HiggsTauTauFramework/CMSSW_7_2_0/src/LLRHiggsTauTau/NtupleProducer/test/HTauTauAnalysis.root"); 
    TTree* treePtr = (TTree*) fIn -> Get ("HTauTauTree/HTauTauTree");
    
    HTauTauTree* tree = new HTauTauTree (treePtr);
    
    // histograms
    TH1D* bDiscr = new TH1D ("h_bDiscr", "bDiscriminator", 20, 0, 10);
    TH1D* h_mothers_px = new TH1D ("h_mothers_px", "mothers_px", 100, -100, 100);
    
    /*
    // Can activate only some branches like this
    tree->GetTree()->SetBranchStatus("*", 0);
    tree->GetTree()->SetBranchStatus("bDiscriminator", 1);
    tree->GetTree()->SetBranchStatus("mothers_px", 1);
    */

    cout << "TreeEntries: " << tree->GetEntries() << endl;
    for (int iEntry = 0; iEntry < tree->GetEntries() ; iEntry++)
    {
        tree->GetEntry(iEntry);
        
        for (int iJet = 0; iJet < tree->bDiscriminator->size(); iJet++)
        {
            bDiscr->Fill (tree->bDiscriminator->at(iJet));
        }
        
        for (int iMoth = 0; iMoth < tree->mothers_px->size(); iMoth++)
        {
            h_mothers_px->Fill(tree->mothers_px->at(iMoth));
        }
    }
    
    TCanvas * c1 = new TCanvas;
    bDiscr->Draw();
    
    TCanvas * c2 = new TCanvas;
    h_mothers_px->Draw();
    
    //c2->Print ("mothers_px.png", "png");
    
}

int main()
{PlotTreeVars();}
