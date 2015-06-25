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


//void MakePlots()
int main()
{
    cout << "Start" << endl;
 
    // lambda 20
    TFile* fIn = new TFile ("/data_CMS/cms/cadamuro/test_submit_to_tier3/HiggsTauTauOutput_HH_Lambda20_NoSvFit_genFix_74X_300000Events_0Skipped_1435155734.19/HH_Lambda20_74X_NoSvFit.root");

    // lambda -4
    //TFile* fIn = new TFile ("/data_CMS/cms/cadamuro/test_submit_to_tier3/HiggsTauTauOutput_HH_Lambdam4_NoSvFit_genFix_74X_300000Events_0Skipped_1435155744.32/HH_Lambdam4_74X_NoSvFit.root");
 
    // lambda 2.46
    //TFile* fIn = new TFile ("/data_CMS/cms/cadamuro/test_submit_to_tier3/HiggsTauTauOutput_HH_Lambda2dot46_NoSvFit_genFix_74X_300000Events_0Skipped_1435182646.12/HH_Lambda2dot46_74X_NoSvFit.root");
 
    // TTbar
    //TFile* fIn = new TFile ("/data_CMS/cms/cadamuro/test_submit_to_tier3/TTBar_74X_NoSvFit_17Giu2015/TTbar_NoSVFit_74X.root");

    const char * decModeNames[6] = {"MuTau", "ETau", "TauTau", "MuMu", "EE", "EMu"}; // follows pairType enum

    
    // histos
    HistoManager* HM [6];
    for (int i = 0; i < 6; i++) HM[i] = new HistoManager (decModeNames[i]); 
    // create list of histos
    for (int i = 0; i < 6; i++) 
    {
        HM[i] -> AddNewHisto ("pTLeg1", "pTLeg1; pT; a.u.", 250, 0, 500);   
        HM[i] -> AddNewHisto ("pTLeg2", "pTLeg2; pT; a.u.", 250, 0, 500);   
        HM[i] -> AddNewHisto ("etaLeg1", "etaLeg1; #eta; a.u.", 100, -3., 3.);   
        HM[i] -> AddNewHisto ("etaLeg2", "etaLeg2; #eta; a.u.", 100, -3., 3.);   
        HM[i] -> AddNewHisto ("isoLeg1", "isoLeg1; iso; a.u.", 100, -10., 90.);
        HM[i] -> AddNewHisto ("isoLeg2", "isoLeg2; iso; a.u.", 100, -10., 90.);
    }
     

    TTree* treePtr = (TTree*) fIn->Get("HTauTauTree/HTauTauTree");
    TH1F *evCounter = (TH1F*) fIn->Get("HTauTauTree/Counters");

    HTauTauTree* tree = new HTauTauTree (treePtr);
    OfflineProducerHelper helper;

    //int nEvents = tree->GetEntries();
    int nEvents = 100000;

    // whatApply: use "OSCharge" (appplies on pairs only)
    // whatApply: use "All", "Iso", "pTMin", "etaMax", "againstEle", "againstMu", "Vertex"; separate various arguments with a semicolon

    // array of subsequent selections; write the selections to apply ("All" applies everything)
    TString sel[] = {"Vertex", "pTMin", "etaMax", "LepID", "againstEle", "againstMu", "OSCharge"};
    int nsel = sizeof(sel)/sizeof(TString);
    TString selConc[nsel];
    Concatenate (sel, nsel, selConc);

    //for (int i = 0; i < nsel; i++)
    //    cout << i << " : " << selConc[i] << endl;

    // counters
    int passReq [6][nsel]; // first index is decay channel, second index is step applied in selection
    int chOverlap [5][nsel]; // for each event, at each selection step, count overlap between categories: 1 - 2 - 3 - 4 - >4 different decay categories accept this event  

    for (int i = 0; i < 6; i++)
        for (int j = 0; j < nsel; j++)
            passReq[i][j] = 0;

    for (int i = 0; i < 5; i++)
        for (int j = 0; j < nsel; j++)
            chOverlap[i][j] = 0;


    bool thisEvPass [6];
   

    for (int iEv = 0; iEv < nEvents; iEv++)
    {
        tree->GetEntry(iEv);
        if (iEv % 10000 == 0) cout << iEv << " / " << nEvents << endl;
        // loop on pairs
        
        // to apply subsequent selections on this event
        for (int iSel = 0; iSel < nsel; iSel++)
        {

            for (int i = 0; i < 6; i++) thisEvPass[i] = false;

            // MC truth
            //int MCDecayTruth = tree->

            // loop on all pairs in the event
            for (int iMoth = 0; iMoth < tree->mothers_px->size(); iMoth++)
            {
                int dau1index = tree->indexDau1->at(iMoth);
                int dau2index = tree->indexDau2->at(iMoth);
                int pairType = helper.getPairType (tree->particleType->at(dau1index), tree->particleType->at(dau2index));

                //cout << tree->particleType->at(dau1index) << " " << tree->particleType->at(dau2index) << endl;
                TLorentzVector v1 (helper.buildDauP4(tree, dau1index));
                TLorentzVector v2 (helper.buildDauP4(tree, dau2index));

                bool PairPassBaseline = helper.pairPassBaseline (tree, iMoth, selConc[iSel]);

                thisEvPass[pairType] = thisEvPass[pairType] || PairPassBaseline;

                if (PairPassBaseline && iSel == nsel-1) // last step    
                {
                    HM[pairType] -> GetHisto ("isoLeg1") -> Fill (tree->daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits->at(dau1index));
                    HM[pairType] -> GetHisto ("isoLeg2") -> Fill (tree->daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits->at(dau2index));
                }

                /*
                HM[pairType] -> GetHisto ("pTLeg1") -> Fill (v1.Pt());
                HM[pairType] -> GetHisto ("pTLeg2") -> Fill (v2.Pt());
                HM[pairType] -> GetHisto ("etaLeg1") -> Fill (v1.Eta());
                HM[pairType] -> GetHisto ("etaLeg2") -> Fill (v2.Eta());
                */
                
                /*
                if (tree->particleType->at(dau1index) == tree->particleType->at(dau2index))
                {
                    TLorentzVector v1 (helper.buildDauP4(tree, dau1index));
                    TLorentzVector v2 (helper.buildDauP4(tree, dau2index));
                    if (tree->particleType->at(dau1index) == OfflineProducerHelper::MUON)
                    {
                        HM.GetHisto ("pT")->Fill (v1.Pt());
                        HM.GetHisto ("eta")->Fill (v1.Eta());
                    }
                }
                */                
            }

            // set counters 
            for (int i = 0; i < 6; i++)
                if (thisEvPass[i]) passReq [i][iSel] ++;
        }
    }
    

    // print all results
    for (int iCh = 0; iCh < 6; iCh++)
    {
        cout << "***  " << decModeNames [iCh] << "  ***" << endl; 
        for (int iSel = 0; iSel < nsel; iSel++)
        {
            cout << sel[iSel] << "    ==> " << passReq [iCh] [iSel] << endl;
        }
        cout << endl << "  ==============================  " << endl;
    }

    
    cout << "Saving histos..." << endl;
    TFile* fOut = new TFile ("Histos_perMode.root", "recreate");
    for (int i = 0; i < 6; i++) HM[i] -> SaveAllToFile (fOut);
    
}
