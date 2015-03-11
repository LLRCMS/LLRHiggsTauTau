/** \class triggerhelper
 *
 *  No description available.
 *
 *  $Date: 2012/06/06 00:27:43 $
 *  $Revision: 1.3 $
 *  \author G. Ortona - LLR
 */
 
#include <DataFormats/Common/interface/TriggerResults.h>
#include <FWCore/Common/interface/TriggerNames.h>
#include "LLRHiggsTauTau/NtupleProducer/interface/triggerhelper.h"
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/TriggerNamesService.h>

using namespace std;
using namespace edm;

triggerhelper::triggerhelper(){
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
  TString tmpMETfilters[nMETs]={
                   "Flag_CSCTightHaloFilter", 
 "Flag_EcalDeadCellTriggerPrimitiveFilter",
                    "Flag_HBHENoiseFilter",  
                         "Flag_METFilters",  
                "Flag_ecalLaserCorrFilter",  
                      "Flag_eeBadScFilter",  
                       "Flag_goodVertices",  
               "Flag_hcalLaserEventFilter",  
              "Flag_trackingFailureFilter",  
                      "Flag_trkPOGFilters",  
     "Flag_trkPOG_logErrorTooManyClusters",  
            "Flag_trkPOG_manystripclus53X",  
         "Flag_trkPOG_toomanystripclus53X",  
  };
  for(int i=0;i<nMETs;i++)metlist[i]=tmpMETfilters[i];
}

int triggerhelper::FindTriggerBit(const edm::Event& event, const vector<string> foundPaths, const vector<int> indexOfPaths){
  
  int bit =0;
  
  edm::Handle<edm::TriggerResults> triggerResults;
  event.getByLabel(InputTag("TriggerResults","","HLT"), triggerResults);
      
  // for(int it=0;it<nTriggers;it++){
  //   for(int j=0;j<(int)triggerNames->size();j++)cout<<triggerNames->triggerName(j)<<endl;
  //   unsigned i =  triggerNames->triggerIndex(triggerlist[it].Data());
  //   if (i< triggerNames->size()){
  //     bit |= 1 << it;   
  //   }
  // }
  // return bit;
  
  for(int it=0;it<nTriggers;it++){
    //for(int j=0;j<(int)foundPaths.size();j++){
    for(int j=0;j<(int)foundPaths.size();j++){
      if(triggerlist[it].CompareTo(foundPaths.at(j))==0){
	      if(triggerResults->accept(indexOfPaths[j]))bit |= 1 <<it;
	      break;
      }
    }
  }
  //printf("bit: %d\n",bit);
  return bit;
}

int triggerhelper::FindMETBit(const edm::Event& event){

  int bit =0;
  edm::Handle<edm::TriggerResults> metFilterBits;
  event.getByLabel(InputTag("TriggerResults","","PAT"), metFilterBits);
  //edm::EDGetTokenT<edm::TriggerResults> metFilterBitsToken_;//(consumes<edm::TriggerResults>(<edm::InputTag>("metFilterBits")));
  //event.getByToken(metFilterBitsToken_, metFilterBits);
  const edm::TriggerNames &metNames = event.triggerNames(*metFilterBits);
  for(int im=0;im<nMETs;im++){
    for(unsigned int i = 0; i < metFilterBits->size(); ++i){
      if(metlist[im].CompareTo(metNames.triggerName(i))==0){
	      if(metFilterBits->accept(i))bit |= 1 <<im;
	      break;
      }
    }
  }
  return bit;
}

int triggerhelper::FindTriggerNumber(TString triggername, bool isTrigger){ 
  int nLoop = nTriggers;
  TString *list = triggerlist;
  if(!isTrigger){
    list = metlist;
    nLoop = nMETs;
  }
  for(int it=0;it<nLoop;it++){ 	
  	if(list[it].CompareTo(triggername.Data())==0)return it;
  }
  return -1;
}

bool triggerhelper::IsTriggerFired(int triggerbit, int triggernumber, bool isTrigger){ 
  int nLoop = nTriggers;
  if(!isTrigger)nLoop = nMETs;
  if(triggernumber>=0 && triggernumber<nLoop) return triggerbit & (1 << triggernumber);
  return false;
}

int triggerhelper::printFiredPaths(int triggerbit, bool isTrigger){
  int nFired = 0;
  int nLoop = nTriggers;
  TString *list = triggerlist;
  if(!isTrigger){
    list = metlist;
    nLoop = nMETs;
  }  
  for(int it=0;it<nLoop;it++){ 	
  	if(IsTriggerFired(triggerbit,it,isTrigger)) {
  	  printf("%s\n",list[it].Data());
  	  nFired++;
  	  }
  }
  return nFired;
}
