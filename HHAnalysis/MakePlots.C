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

//void MakePlots()
int main()
{
    cout << "Start" << endl;
    // lambda 1
    TFile* fIn = new TFile ("/home/llr/cms/cadamuro/HiggsTauTauFramework/CMSSW_7_2_3_patch1/src/LLRHiggsTauTau/NtupleProducer/test/HTauTauAnalysis_ntuple_HH_Lambda1_3Giu2015.root");
    
    // lambda 20
    //TFile* fIn = new TFile ("/home/llr/cms/cadamuro/HiggsTauTauFramework/CMSSW_7_2_3_patch1/src/LLRHiggsTauTau/NtupleProducer/test/HTauTauAnalysis_ntuple_HH_Lambda20_3Giu2015.root");

    // histos
    const char * decModeNames[6] = {"MuTau", "ETau", "TauTau", "MuMu", "EE", "EMu"}; // follows pairType enum
    HistoManager* HM [6];
    for (int i = 0; i < 6; i++) HM[i] = new HistoManager (decModeNames[i]); 
    // create list of histos
    for (int i = 0; i < 6; i++) 
    {
        HM[i] -> AddNewHisto ("pTLeg1", "pTLeg1; pT; a.u.", 100, 0, 300);   
        HM[i] -> AddNewHisto ("pTLeg2", "pTLeg2; pT; a.u.", 100, 0, 300);   
        HM[i] -> AddNewHisto ("etaLeg1", "etaLeg1; #eta; a.u.", 100, -3., 3.);   
        HM[i] -> AddNewHisto ("etaLeg2", "etaLeg2; #eta; a.u.", 100, -3., 3.);   
    }

    TTree* treePtr = (TTree*) fIn->Get("Ntuplizer/HTauTauTree");
    TH1F *evCounter = (TH1F*) fIn->Get("Ntuplizer/Counters");
    
    HTauTauTree* tree = new HTauTauTree (treePtr);
    OfflineProducerHelper helper;

    int nEvents = tree->GetEntries();
    for (int iEv = 0; iEv < nEvents; iEv++)
    {
        tree->GetEntry(iEv);
        if (iEv % 10000 == 0) cout << iEv << " / " << nEvents << endl;
        // loop on pairs
        for (int iMoth = 0; iMoth < tree->mothers_px->size(); iMoth++)
        {
            int dau1index = tree->indexDau1->at(iMoth);
            int dau2index = tree->indexDau2->at(iMoth);
            int pairType = helper.getPairType (tree->particleType->at(dau1index), tree->particleType->at(dau2index));
            //cout << tree->particleType->at(dau1index) << " " << tree->particleType->at(dau2index) << endl;
            TLorentzVector v1 (helper.buildDauP4(tree, dau1index));
            TLorentzVector v2 (helper.buildDauP4(tree, dau2index));

            HM[pairType] -> GetHisto ("pTLeg1") -> Fill (v1.Pt());
            HM[pairType] -> GetHisto ("pTLeg2") -> Fill (v2.Pt());
            HM[pairType] -> GetHisto ("etaLeg1") -> Fill (v1.Eta());
            HM[pairType] -> GetHisto ("etaLeg2") -> Fill (v2.Eta());
            
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
    }
    
    cout << "Saving histos..." << endl;
    TFile* fOut = new TFile ("Histos_perMode.root", "recreate");
    for (int i = 0; i < 6; i++) HM[i] -> SaveAllToFile (fOut);
}
