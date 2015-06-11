#ifndef OfflineHelper_h
#define OfflineHelper_h
/** \class OfflineProducerHelper
 *
 *  Collection of functions to analyze HTauTau primary trees, without CMSSW dependencies.
 *
 *  $Date: 2012/06/06 00:27:43 $
 *  $Revision: 1.3 $
 *  \author G. Ortona - LLR
 */
 
#include <TString.h>
#include "TMath.h"
#include "HTauTauTree.h"
#include "TLorentzVector.h"
#include <iostream>

using namespace std;
 
class OfflineProducerHelper {
 public:
  enum particleType {
  MUON = 0,
  ELECTRON = 1,
  TAU =2
  };

  enum pairType {
    MuHad  = 0,
    EHad   = 1,
    HadHad = 2,
    MuMu   = 3,
    EE     = 4,
    EMu    = 5,
    EEPrompt = 6, // prompt Z->ee/mumu decays
    MuMuPrompt = 7,
    Other  = 8 // for e.g. h->bb
  };

  OfflineProducerHelper();
  int FindTriggerNumber(TString triggername);
  bool IsTriggerFired(int triggerbit, int triggerNumber);
  bool IsTriggerFired(int triggerbit, TString triggerName){return IsTriggerFired(triggerbit, FindTriggerNumber(triggerName));}
  int printFiredPaths(int);
  bool isMuon(int type){if(type == MUON)return true; else return false;}
  bool isElectron(int type){if(type == ELECTRON)return true; else return false;}
  bool isTau(int type){if(type == TAU)return true; else return false;}
  int getPairType (int type1, int type2); // return pair type
  bool checkBit (int word, int bitpos); // check bit "bitpos" in a word
  bool pairPassBaseline (HTauTauTree* tree, int iPair);
  bool eleBaseline (HTauTauTree* tree, int iDau, float ptMin, float relIso); // return true if leptons passes the baseline selections
  bool muBaseline (HTauTauTree* tree, int iDau, float ptMin, float etaMax, float relIso);
  bool tauBaseline (HTauTauTree* tree, int iDau, float ptMin, float etaMax, int againstEleWP, int againstMuWP, float isoRaw3Hits);
  bool tightEleMVAID (float BDT, float fSCeta); // compute tight ele MVA id WP, but isBDT in ntuples has been fixed --> this will be soon deprecated
                                                // approx!! I call it using lepton eta and not superCluster eta
  TLorentzVector buildDauP4 (HTauTauTree* tree, int iDau); // build daughter 4 vector
  ~OfflineProducerHelper(){}

 private:
  static const int nTriggers =19;
  TString triggerlist[nTriggers];
  


};

OfflineProducerHelper::OfflineProducerHelper(){
  TString tmptrigger[nTriggers]={
    "IsoMu17_eta2p1_LooseIsoPFTau20",
    "IsoMu17_eta2p1",
    "IsoMu17_eta2p1_MediumIsoPFTau40_Trk1_eta2p1_Reg",
    "IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1",
    "IsoMu24_eta2p1_IterTrk01",
    "IsoMu24_eta2p1_IterTrk02",
    "IsoMu24_eta2p1_IterTrk02_LooseIsoPFTau20",
    "Ele22_eta2p1_WP85_Gsf_LooseIsoPFTau20",
    "Ele32_eta2p1_WP85_Gsf",
    "Ele32_eta2p1_WP85_Gsf_LooseIsoPFTau20",
    "LooseIsoPFTau50_Trk30_eta2p1_MET120",
    "IsoMu16_eta2p1_CaloMET30_LooseIsoPFTau50_Trk30_eta2p1",
    "IsoMu16_eta2p1_CaloMET30",
    "Mu16_eta2p1_CaloMET30",
    "LooseIsoPFTau50_Trk30_eta2p1",
    "DoubleIsoMu17_eta2p1",
    "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",
    "Ele27_eta2p1_WP85_Gsf_LooseIsoPFTau20",
    "Ele27_eta2p1_WP85_Gsf"
  };
  for(int i=0;i<nTriggers;i++){
    triggerlist[i]=tmptrigger[i];
    triggerlist[i].Prepend("HLT_");
    triggerlist[i].Append("_v1");
  }
}

int OfflineProducerHelper::FindTriggerNumber(TString triggername){
  for(int it=0;it<nTriggers;it++){ 	
  	if(triggerlist[it].CompareTo(triggername.Data())==0)return it;
  	else {
  	    TString newName=triggername.Data();
  	    newName.Prepend("HLT_");
  	    newName.Append("_v1");
  	    if(triggerlist[it].CompareTo(newName.Data())==0)return it;
  	}
  }
  return -1;
}

bool OfflineProducerHelper::IsTriggerFired(int triggerbit, int triggernumber){ 
  if(triggernumber>=0 && triggernumber<nTriggers) return triggerbit & (1 << triggernumber);
  return false;
}

int OfflineProducerHelper::printFiredPaths(int triggerbit){
  int nFired = 0;
  for(int it=0;it<nTriggers;it++){ 	
  	if(IsTriggerFired(triggerbit,it)) {
  	  printf("%s\n",triggerlist[it].Data());
  	  nFired++;
  	  }
  }
  return nFired;
}

bool OfflineProducerHelper::checkBit (int word, int bitpos)
{
    bool res = word & (1 << bitpos);
    return res;
}

int OfflineProducerHelper::getPairType (int type1, int type2)
{
    int nmu = 0;
    int nele = 0;
    int ntau = 0;
    
    if (isMuon (type1) )     nmu++;
    if (isElectron (type1) ) nele++;
    if (isTau (type1) )      ntau++;

    if (isMuon (type2) )     nmu++;
    if (isElectron (type2) ) nele++;
    if (isTau (type2) )      ntau++;

    if (nmu == 1 && nele == 0 && ntau == 1) return (int) pairType::MuHad;
    if (nmu == 0 && nele == 1 && ntau == 1) return (int) pairType::EHad;
    if (nmu == 0 && nele == 0 && ntau == 2) return (int) pairType::HadHad;
    if (nmu == 2 && nele == 0 && ntau == 0) return (int) pairType::MuMu;
    if (nmu == 0 && nele == 2 && ntau == 0) return (int) pairType::EE;
    if (nmu == 1 && nele == 1 && ntau == 0) return (int) pairType::EMu;
    
    return -1;

}

bool OfflineProducerHelper::pairPassBaseline (HTauTauTree* tree, int iPair)
{
    int dau1index = tree->indexDau1->at(iPair);
    int dau2index = tree->indexDau2->at(iPair);
    int type1 = tree->particleType->at(dau1index);
    int type2 = tree->particleType->at(dau2index);
    int pairType = getPairType (type1, type2);

    // pairs are always ordered as: e mu | e tau | mu tau  (e < mu < tau)
    // if same type of particle, highest pt one is the first
    if (pairType == MuHad)
    {
        bool leg1 = muBaseline (tree, dau1index, 18., 2.1, 0.1);
        bool leg2 = tauBaseline (tree, dau2index, 20., 2.3, 0, 1, 1.5);
        return (leg1 && leg2);
    }

    if (pairType == EHad)
    {
        bool leg1 = eleBaseline (tree, dau1index, 23., 0.1);
        bool leg2 = tauBaseline (tree, dau2index, 20., 2.3, 3, 0, 1.5);
        return (leg1 && leg2);
    }

    // ordered by pT and not by most isolated, but baseline asked in sync is the same...
    if (pairType == HadHad)
    {
        bool leg1 = tauBaseline (tree, dau1index, 45., 2.1, 0, 0, 1.0);
        bool leg2 = tauBaseline (tree, dau2index, 45., 2.1, 0, 0, 1.0);
        return (leg1 && leg2);
    }

    if (pairType == EMu)
    {
        bool leg1 = eleBaseline (tree, dau1index, 13., 0.15);
        bool leg2 = muBaseline (tree, dau2index, 9., 2.4, 0.15);
        return (leg1 && leg2);
    }
    
    // e e, mu mu missing for the moment...
    if (pairType == EE) return false;
    if (pairType == MuMu) return false;
    
    return false;
        
}


bool OfflineProducerHelper::eleBaseline (HTauTauTree* tree, int iDau, float ptMin, float relIso)
{
    float px = tree->daughters_px->at(iDau);
    float py = tree->daughters_py->at(iDau);
    float pz = tree->daughters_pz->at(iDau);
    float e =  tree->daughters_e->at(iDau);

    TLorentzVector p4 (px, py, pz, e);
    
    bool vertexS = (tree->dxy->at(iDau) < 0.045 && tree->dz->at(iDau) < 0.2);
    //bool idS = checkBit (tree->daughters_iseleCUT->at(iDau), 3); // 3 is TIGHT ele id CUT BASED
    //bool idS = tree->daughters_iseleBDT->at(iDau); // use it in ntuples produced after 11 June 2015, contains tight WP bool  
    bool idS = tightEleMVAID (tree->discriminator->at(iDau), TMath::Abs(p4.Eta())); // APPROX! Using lepton eta and not super cluster eta, discriminator contains ele BDT  
    bool isoS = (tree->combreliso->at(iDau) < relIso);
    bool ptS = (p4.Pt() > ptMin);
    bool etaS = (fabs(p4.Eta()) < 2.5);
    
    bool totalS = (vertexS && idS && isoS && ptS && etaS);
    return totalS;
    
}

bool OfflineProducerHelper::muBaseline (HTauTauTree* tree, int iDau, float ptMin, float etaMax, float relIso)
{
    float px = tree->daughters_px->at(iDau);
    float py = tree->daughters_py->at(iDau);
    float pz = tree->daughters_pz->at(iDau);
    float e =  tree->daughters_e->at(iDau);

    TLorentzVector p4 (px, py, pz, e);
    int discr = tree->daughters_muonID->at(iDau);
        
    bool vertexS = (tree->dxy->at(iDau) < 0.045 && tree->dz->at(iDau) < 0.2);
    bool idS = checkBit (discr, 2); // bit 2 is MEDIUM mu id
    bool isoS = (tree->combreliso->at(iDau) < relIso);
    bool ptS = (p4.Pt() > ptMin);
    bool etaS = (fabs(p4.Eta()) < etaMax);
    
    bool totalS = (vertexS && idS && isoS && ptS && etaS);
    return totalS;
}

// againstEleWP: 0 = VLoose, 1 = Loose, 2 = Medium, 3 = Tight, 4 = VTight [all are MVA discr]
// againstMuWP: 0 = Loose, 1 = Tight
bool OfflineProducerHelper::tauBaseline (HTauTauTree* tree, int iDau, float ptMin, float etaMax, int againstEleWP, int againstMuWP, float isoRaw3Hits)
{
    float px = tree->daughters_px->at(iDau);
    float py = tree->daughters_py->at(iDau);
    float pz = tree->daughters_pz->at(iDau);
    float e =  tree->daughters_e->at(iDau);

    TLorentzVector p4 (px, py, pz, e);

    if (againstEleWP < 0 || againstEleWP > 4) {
        cout << " ** OfflineProducerHelper::tauBaseline: againstEleWP must be between 0 and 4 --> using 0" << endl;
        againstEleWP = 0;
    } 

    if (againstMuWP < 0 || againstEleWP > 1) {
        cout << " ** OfflineProducerHelper::tauBaseline: againstMuWP must be between 0 and 1 --> using 0" << endl;
        againstMuWP = 0;
    }
    
    int agEleVal = 0;
    int agMuVal = 0;
    
    // ag ele:
    if (againstEleWP == 0)      agEleVal = tree->daughters_againstElectronVLooseMVA5->at(iDau);
    else if (againstEleWP == 1) agEleVal = tree->daughters_againstElectronLooseMVA5->at(iDau);
    else if (againstEleWP == 2) agEleVal = tree->daughters_againstElectronMediumMVA5->at(iDau);
    else if (againstEleWP == 3) agEleVal = tree->daughters_againstElectronTightMVA5->at(iDau);
    else if (againstEleWP == 4) agEleVal = tree->daughters_againstElectronVTightMVA5->at(iDau);   

    // ag mu:
    if (againstMuWP == 0)      agMuVal = tree->daughters_againstMuonLoose3->at(iDau);
    else if (againstMuWP == 1) agMuVal = tree->daughters_againstMuonTight3->at(iDau);

    bool dmfS = (tree->daughters_decayModeFindingOldDMs->at(iDau) == 1 || tree->daughters_decayModeFindingNewDMs->at(iDau) == 1);
    bool vertexS = (tree->dxy->at(iDau) < 0.045 && tree->dz->at(iDau) < 0.2);
    bool agEleS = (agEleVal == 1); 
    bool agMuS  = (agMuVal == 1); 
    bool isoS = (tree->daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits->at(iDau) < isoRaw3Hits);
    bool ptS = (p4.Pt() > ptMin);
    bool etaS = (fabs(p4.Eta()) < etaMax);

    bool totalS = (dmfS && vertexS && agEleS && agMuS && isoS && ptS && etaS);
    return totalS;    
}

// branch isBDT has been updated, this will soon be unnecessary (and eventually obsolete)
bool OfflineProducerHelper::tightEleMVAID (float BDT, float fSCeta)
{
    // POG_MVA_ID_Run2_NonTrig_Tight PHYS14
    if (fSCeta < 0) fSCeta *= -1;
    bool isBDT = false;
    if (fSCeta <0.8) isBDT = (BDT>0.73);
    else if (fSCeta < 1.479) isBDT = (BDT>0.57);
    else isBDT = (BDT >0.05);
    
    return isBDT;
}

TLorentzVector OfflineProducerHelper::buildDauP4 (HTauTauTree* tree, int iDau)
{
    float px = tree->daughters_px->at(iDau);
    float py = tree->daughters_py->at(iDau);
    float pz = tree->daughters_pz->at(iDau);
    float e =  tree->daughters_e->at(iDau);

    TLorentzVector p4 (px, py, pz, e);
    return p4;
}

#endif
