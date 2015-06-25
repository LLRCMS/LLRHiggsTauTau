// make use of MC truth to control the efficiency in pair selection
// and the overlap between different categories

#ifndef __CINT__
#include <iostream>
#include "TFile.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#endif

#include "../NtupleProducer/Utils/OfflineProducerHelper.h"
#include "../NtupleProducer/interface/GenFlags.h"
#include "../NtupleProducer/Utils/HTauTauTree.h"
#include "HistoManager.h"

// merge array of selection strings
void Concatenate (TString* sel, int nsel, TString* selConc)
{
    for (int i = 0; i < nsel; i++)
    {
        TString thisName = "";
        for (int j = 0; j <=i; j++)
        {
            thisName += sel[j];
            thisName += ";" ;
        }
        selConc [i] = thisName;
    }
    return;
}


int main()
{
    // lambda 20
    TFile* fIn = new TFile ("/data_CMS/cms/cadamuro/test_submit_to_tier3/HiggsTauTauOutput_HH_Lambda20_NoSvFit_genFix_74X_300000Events_0Skipped_1435155734.19/HH_Lambda20_74X_NoSvFit.root");

    // lambda -4
    //TFile* fIn = new TFile ("/data_CMS/cms/cadamuro/test_submit_to_tier3/HiggsTauTauOutput_HH_Lambdam4_NoSvFit_genFix_74X_300000Events_0Skipped_1435155744.32/HH_Lambdam4_74X_NoSvFit.root");
 
    // lambda 2.46
    //TFile* fIn = new TFile ("/data_CMS/cms/cadamuro/test_submit_to_tier3/HiggsTauTauOutput_HH_Lambda2dot46_NoSvFit_genFix_74X_300000Events_0Skipped_1435182646.12/HH_Lambda2dot46_74X_NoSvFit.root");
    
    const char * decModeNames[6] = {"MuTau", "ETau", "TauTau", "MuMu", "EE", "EMu"}; // follows pairType enum
 
    TTree* treePtr = (TTree*) fIn->Get("HTauTauTree/HTauTauTree");
    TH1F *evCounter = (TH1F*) fIn->Get("HTauTauTree/Counters");

    HTauTauTree* tree = new HTauTauTree (treePtr);
    OfflineProducerHelper helper;

    TString sel[] = {"Vertex", "pTMin", "etaMax", "LepID", "againstEle", "againstMu", "Iso", "OSCharge"};
    int nsel = sizeof(sel)/sizeof(TString);
    TString selConc[nsel];
    Concatenate (sel, nsel, selConc);

    int nEvents = tree->GetEntries();
    //int nEvents = 100000;
    
    int dec[6] = {0,0,0,0,0,0}; // counts MC decays
    int cat[6][7][nsel]; // [MC true decay][associated category+no category][selection step]
    bool thisEvPass [6];
    int nGoodPairs[6]; // count good pairs for each reco pair type
    
    int categorization [6][][nsel];
    
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 7; j++)
            for (int z = 0; z < nsel; z++)
                cat[i][j][z] = 0;

    
    // events are divided according to their MC truth decay
    for (int iEv = 0; iEv < nEvents; iEv++)
    {
        tree->GetEntry(iEv);
        if (iEv % 10000 == 0) cout << iEv << " / " << nEvents << endl;
        // loop on pairs
       
        int decay = helper.MCHiggsTauTauDecayMode (tree);
        if (decay == -1)
        {
            cout << "** WARNING! Couldn't find gen decay mode, skipping event... **" << endl;
            continue;
        }
        dec[decay] ++;
        
        // to apply subsequent selections on this event
        for (int iSel = 0; iSel < nsel; iSel++)
        {
            for (int i = 0; i < 6; i++){
                thisEvPass[i] = false;
                nGoodPairs[i] = 0;
            }
            // loop on all pairs in the event
            for (int iMoth = 0; iMoth < tree->mothers_px->size(); iMoth++)
            {
                int dau1index = tree->indexDau1->at(iMoth);
                int dau2index = tree->indexDau2->at(iMoth);
                int pairType = helper.getPairType (tree->particleType->at(dau1index), tree->particleType->at(dau2index));

                bool PairPassBaseline = helper.pairPassBaseline (tree, iMoth, selConc[iSel]);
                if (PairPassBaseline) nGoodPairs[pairType]++;                
                thisEvPass [pairType] = thisEvPass [pairType] || PairPassBaseline;
            }
            
            // now I know the MC decay and if, for this selection, an event enters a certain category
            
        }

    }
    
    for (int i = 0; i < 6; i++) cout << i << " --> " << dec[i] << endl;

}
