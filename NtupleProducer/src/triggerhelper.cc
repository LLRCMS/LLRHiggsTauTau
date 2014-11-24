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
  for(int i=0;i<19;i++)triggerlist[i]=tmptrigger[i];
}

int triggerhelper::FindTriggerBit(const edm::Event& event){
  int bit =0;

  edm::Handle<edm::TriggerResults> triggerResults;
  const edm::TriggerNames* triggerNames;

  event.getByLabel(InputTag("TriggerResults"), triggerResults);
  triggerNames = &(event.triggerNames(*triggerResults));
  //from this I'm getting the paths, not the triggers :(


  //  edm::Handle<edm::TriggerResults> myTriggerResults;
  //  const edm::Triggernames* myTriggerNames=0;

  for(int it=0;it<nTriggers;it++){
    for(int j=0;j<(int)triggerNames->size();j++)cout<<triggerNames->triggerName(j)<<endl;
    unsigned i =  triggerNames->triggerIndex(triggerlist[it].Data());
    if (i< triggerNames->size()){
      bit |= 1 << it;   
    }
  }
  return bit;
}
