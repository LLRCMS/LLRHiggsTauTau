#include <iostream>
#include "TFile.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TString.h"
#include <utility>
#include "../NtupleProducer/Utils/OfflineProducerHelper.h"
#include "../NtupleProducer/interface/GenFlags.h"
#include "../NtupleProducer/Utils/HTauTauTree.h"
#include "HHKinFit/interface/HHKinFitMaster.h"
//#include "HistoManager.h"

// ================================
// INSTRUCTIONS:
// cmsenv
// ./compileHHKinFitProducer.sh
// export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(pwd)  [just do once if re-compiling]
// launch with ./ComputeHHKinFit
// ================================

using namespace std;

int main()
{
    cout << "Start HH computation..." << endl;

    // lambda 1
    // waiting for the file...
        
    // lambda 20
    TFile* fIn = new TFile ("/data_CMS/cms/cadamuro/test_submit_to_tier3/HiggsTauTauOutput_HH_Lambda20_NoSvfit_Prod5Lug2015_300000Events_0Skipped_1436090974.69/HH_Lambda20_NoSVFit_5Lug2015.root"); TString outName = "HHLambda20_5Lug_";

    // lambda 2.46
    //TFile* fIn = new TFile ("/data_CMS/cms/cadamuro/test_submit_to_tier3/HiggsTauTauOutput_HH_Lambda2dot46_NoSvFit_Prod5Lug2015_300000Events_0Skipped_1436090991.56/HH_Lambda2dot46_NoSVFit_5Lug2015.root"); TString outName = "HHLambda2dot46_";
    
    // lambda -4
    //TFile* fIn = new TFile ("/data_CMS/cms/cadamuro/test_submit_to_tier3/HiggsTauTauOutput_HH_Lambdam4_NoSvFit_Prod5Lug2015_300000Events_0Skipped_1436091004.38/HH_Lambdam4_NoSVFit_5Lug2015.root"); TString outName = "HHLambdam4_";
    
    // ttbar
    //TFile* fIn = new TFile ("/home/llr/cms/cadamuro/TTJets_Govoni.root"); TString outName = "TTBar_";
    
    TTree* treePtr = (TTree*) fIn->Get("HTauTauTree/HTauTauTree");
    TH1F *evCounter = (TH1F*) fIn->Get("HTauTauTree/Counters");

    HTauTauTree* tree = new HTauTauTree (treePtr);
    OfflineProducerHelper helper;

    int nEvents = tree->GetEntries();
    //int nEvents = 10000;

    TH1D* HHMass = new TH1D ("HHMass", "HHMass", 200, 200, 800);

    //define the testd hypotheses
    std::vector<Int_t> hypo_mh1;
    hypo_mh1.push_back(125);
    std::vector<Int_t> hypo_mh2;
    hypo_mh2.push_back(125);

    for (int iEv = 0; iEv < nEvents; iEv++)
    {
        if (iEv % 5000 == 0) cout << iEv << " / " << nEvents << endl;
    
        tree->GetEntry(iEv);
        //define input vectors
        int jet1 = -1;
        int jet2 = -1;
        
        bool jetFound = helper.getBestJets (tree, jet1, jet2, 0); // 0 --> ordered by CVS score
        if (!jetFound) continue;  
        
        // ask for 2 MEDIUM b tagged jets
        if (tree->bCSVscore->at(jet2) <= 0.89)
            continue;
        
        // pritn selected b jets among jet list
        /*
        for (int i = 0; i < tree->bCSVscore->size(); i++)
        {
            cout << tree->bCSVscore->at(i) << " ";
            if (i == jet1 || i == jet2) cout << " <== " << (i == jet2);
            cout << endl;
        }
        cout << endl;        
        */
        
        TLorentzVector b1      = TLorentzVector(tree->jets_px->at(jet1), tree->jets_py->at(jet1), tree->jets_pz->at(jet1), tree->jets_e->at(jet1));
        TLorentzVector b2      = TLorentzVector(tree->jets_px->at(jet2), tree->jets_py->at(jet2), tree->jets_pz->at(jet2), tree->jets_e->at(jet2));
               
        // do just HH kin fit for pair 0, to test
        if (tree->indexDau1->size() == 0) continue; // no pairs here
        
        int iMoth = -1;
        for (int iPair = 0; iPair < tree->mothers_px->size(); iPair++)
        {
            if (helper.pairPassBaseline (tree, iPair) )
            {
                iMoth = iPair; 
                break; // one is enough for this test
            }
        }
        if (iMoth == -1) continue;
        
        int dau1index = tree->indexDau1->at(iMoth); // take first
        int dau2index = tree->indexDau2->at(iMoth);
        TLorentzVector tau1vis = TLorentzVector(tree->daughters_px->at(dau1index), tree->daughters_py->at(dau1index), tree->daughters_pz->at(dau1index), tree->daughters_e->at(dau1index));
        TLorentzVector tau2vis = TLorentzVector(tree->daughters_px->at(dau2index), tree->daughters_py->at(dau2index), tree->daughters_pz->at(dau2index), tree->daughters_e->at(dau2index));
 
        float METx = tree->METx->at(iMoth);
        float METy = tree->METy->at(iMoth);
        float METpt = TMath::Sqrt(METx*METx + METy*METy);
      
        TLorentzVector ptmiss  = TLorentzVector(METx, METy,0, METpt);
        TMatrixD metcov(2,2);
        metcov(0,0)=tree->MET_cov00->at(iMoth);
        metcov(1,0)=tree->MET_cov10->at(iMoth);
        metcov(0,1)=tree->MET_cov01->at(iMoth);    
        metcov(1,1)=tree->MET_cov11->at(iMoth);
    
        //intance of fitter master class
        HHKinFitMaster kinFits = HHKinFitMaster(&b1,&b2,&tau1vis,&tau2vis);
        kinFits.setAdvancedBalance(&ptmiss,metcov);
        //kinFits.setSimpleBalance(ptmiss.Pt(),10); //alternative which uses only the absolute value of ptmiss in the fit
        kinFits.addMh1Hypothesis(hypo_mh1);
        kinFits.addMh2Hypothesis(hypo_mh2);
        kinFits.doFullFit();
    
        //obtain results from different hypotheses
        Double_t chi2_best = kinFits.getBestChi2FullFit();
        Double_t mh_best = kinFits.getBestMHFullFit();
        std::pair<Int_t, Int_t> bestHypo = kinFits.getBestHypoFullFit();
        std::map< std::pair<Int_t, Int_t>, Double_t> fit_results_chi2 = kinFits.getChi2FullFit();
        std::map< std::pair<Int_t, Int_t>, Double_t> fit_results_fitprob = kinFits.getFitProbFullFit();
        std::map< std::pair<Int_t, Int_t>, Double_t> fit_results_mH = kinFits.getMHFullFit();
        std::map< std::pair<Int_t, Int_t>, Double_t> fit_results_pull_b1 = kinFits.getPullB1FullFit();
        std::map< std::pair<Int_t, Int_t>, Double_t> fit_results_pull_b2 = kinFits.getPullB2FullFit();
        std::map< std::pair<Int_t, Int_t>, Double_t> fit_results_pull_balance = kinFits.getPullBalanceFullFit();
        std::map< std::pair<Int_t, Int_t>, Int_t> fit_convergence = kinFits.getConvergenceFullFit();

        //cout << mh_best << " with chi2 best of " << chi2_best << endl;
        HHMass -> Fill (mh_best);
    }
    
    TFile* fOut = new TFile (outName + "HHMass.root", "recreate");
    HHMass->Write();
}
