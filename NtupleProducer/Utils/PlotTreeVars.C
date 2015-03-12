#include "HTauTauTree.h"
#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "OfflineProducerHelper.h"

//For trigger selction tools, please refer to 
//https://github.com/LLRCMS/LLRHiggsTauTau/blob/master/NtupleProducer/Utils/OfflineProducerHelper.h

using namespace std;

void PlotTreeVars()
{
    TFile* fIn = new TFile ("../test/HTauTauAnalysis.root"); 
    
    TTree* treePtr = (TTree*) fIn->Get("HTauTauTree/HTauTauTree");
    TH1F *evCounter = (TH1F*) fIn->Get("HTauTauTree/Counters");
    
    cout<<"Number of Events analysed "<<evCounter->GetBinContent(1)<<endl;
    cout<<"Number of Events Passing at least one trigger "<<evCounter->GetBinContent(2)<<endl;
    
    //HTauTauTree class defines the branches from the Ntuplizer output
    HTauTauTree* tree = new HTauTauTree (treePtr);
    
    //trigger helper 
    OfflineProducerHelper helper;
    /*
    // Can activate only some branches like this
    tree->GetTree()->SetBranchStatus("*", 0);
    tree->GetTree()->SetBranchStatus("bDiscriminator", 1);
    tree->GetTree()->SetBranchStatus("mothers_px", 1);
    */

    // histograms
    TH1F* bDiscr = new TH1F ("h_bDiscr", "bDiscriminator", 20, 0, 10);
    TH1F* h_mothers_px = new TH1F ("h_mothers_px", "mothers_px", 100, -100, 100);
    TH1F* h_mothers_SVmass = new TH1F ("h_mothers_SVfitMass", "mothers_SVfitMass", 100, -100, 800);
    TH2F *h_dauPtvsSVFit = new TH2F("dauPtvsSVMass","dauPtvsSVMass",100,-100,800,100,0,1000);
    cout << "TreeEntries: " << tree->GetEntries() << endl;
    for (int iEntry = 0; iEntry < tree->GetEntries() ; iEntry++)
    {
        tree->GetEntry(iEntry);

        //Check if a specific trigger is fired
        if(helper.IsTriggerFired(tree->triggerbit,"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v1")){
            cout<<"trigger HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v1 fired in event "<<iEntry<<endl;
        }
        //print all fired triggers
        //cout<<"NEW EVENT"<<endl<<endl;
        //helper.printFiredPaths(tree->triggerbit);
        
        //Distribution of the bDiscriminant
        for (int iJet = 0; iJet < tree->bDiscriminator->size(); iJet++)
        {
            bDiscr->Fill (tree->bDiscriminator->at(iJet));
        }
        
        //Distributions for the pairs and for the daughters
        for (int iMoth = 0; iMoth < tree->mothers_px->size(); iMoth++)
        {
            h_mothers_px->Fill(tree->mothers_px->at(iMoth));
            float svmass=tree->SVfitMass->at(iMoth);
            h_mothers_SVmass->Fill(svmass);
            //Distribution of daughters pt as a function of the pair invariant mass
            //indexes of the daughters (in the daughters vector)
            int dau1index = tree->indexDau1->at(iMoth);
            int dau2index = tree->indexDau2->at(iMoth);
            //daughters four momenta
            float dau1px = tree->daughters_px->at(dau1index);
            float dau1py = tree->daughters_py->at(dau1index);
            float dau1pz = tree->daughters_pz->at(dau1index);
            float dau1e = tree->daughters_e->at(dau1index);
            float dau2px = tree->daughters_px->at(dau2index);
            float dau2py = tree->daughters_py->at(dau2index);
            float dau2pz = tree->daughters_pz->at(dau2index);
            float dau2e = tree->daughters_e->at(dau2index);
            TLorentzVector dau1(dau1px, dau1py, dau1pz, dau1e);
            TLorentzVector dau2(dau2px, dau2py, dau2pz, dau2e);
            h_dauPtvsSVFit->Fill(svmass, dau1.Pt());
            h_dauPtvsSVFit->Fill(svmass, dau2.Pt());
        }
        
    }
    
    TCanvas * c1 = new TCanvas;
    bDiscr->Draw();
    
    TCanvas * c2 = new TCanvas;
    h_mothers_px->Draw();
    
    TCanvas * c3 = new TCanvas();
    h_dauPtvsSVFit->Draw("COLZ");
    
    TCanvas *c4 = new TCanvas();
    h_mothers_SVmass->Draw();
    
    //c2->Print ("mothers_px.png", "png");
    
}

// if compiling: c++ -lm -o PlotTreeVars PlotTreeVars.C `root-config --glibs --cflags`
int main()
{PlotTreeVars();}
