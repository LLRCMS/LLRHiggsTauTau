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
 
class OfflineProducerHelper {
 public:
  enum particleType {
  MUON = 0,
  ELECTRON = 1,
  TAU =2
  };

  OfflineProducerHelper();
  int FindTriggerNumber(TString triggername);
  bool IsTriggerFired(int triggerbit, int triggerNumber);
  bool IsTriggerFired(int triggerbit, TString triggerName){return IsTriggerFired(triggerbit, FindTriggerNumber(triggerName));}
  int printFiredPaths(int);
  bool isMuon(int type){if(type == MUON)return true; else return false;}
  bool isElectron(int type){if(type == ELECTRON)return true; else return false;}
  bool isTau(int type){if(type == TAU)return true; else return false;}

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
  	if(triggerlist[it].Data()==triggername.Data())return it;
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


#endif
