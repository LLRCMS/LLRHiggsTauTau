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
#include <iostream>
#include <vector>

using namespace std;
using namespace edm;

triggerhelper::triggerhelper(std::vector<std::string> HLTPaths) : nTriggers(HLTPaths.size())
{

  //cout << "nTriggers: " << nTriggers << endl;
  bool found = false;
  triggerlist = new TString [nTriggers];
  //triggerMap = vector<triggerMapper>;
  
  for (int iHLT = 0; iHLT < nTriggers; ++iHLT)
  {
    triggerlist[iHLT] = HLTPaths.at(iHLT);
    TString dummyList[1] = {"dummy"};
    if(triggerlist[iHLT]=="HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1"){
    	triggerMap.push_back(triggerMapper(triggerlist[iHLT],
                      "hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter",
					"hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23",triggerMapper::kemu));
	found=true;
    }
    if(triggerlist[iHLT]=="HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v1"){
      triggerMap.push_back(triggerMapper(triggerlist[iHLT]
				      ,"hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter",
				      "hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8",triggerMapper::kemu));
      found=true;
    }
    if(triggerlist[iHLT]=="HLT_IsoMu24_eta2p1_v1"){
      triggerMap.push_back(triggerMapper(triggerlist[iHLT]
				      ,"",
				      "hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09",triggerMapper::kemu));
      found=true;
    }
    if(triggerlist[iHLT]=="HLT_IsoMu27_v1"){
      triggerMap.push_back(triggerMapper(triggerlist[iHLT]
				      ,"",
				      "hltL3crIsoL1sMu25L1f0L2f10QL3f27QL3trkIsoFiltered0p09",triggerMapper::kemu));
      found=true;
    }
    if(triggerlist[iHLT]=="HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v1"){
      TString listfilt2[2]={"hltPFTau20TrackLooseIsoAgainstMuon ","hltOverlapFilterIsoMu17LooseIsoPFTau20"};
      TString listfilt1[2]={"hltL3crIsoL1sMu16erTauJet20erL1f0L2f10QL3f17QL3trkIsoFiltered0p09","hltOverlapFilterIsoMu17LooseIsoPFTau20"};
      triggerMap.push_back(triggerMapper(triggerlist[iHLT]
				      ,listfilt1
				      ,listfilt2,2,2,triggerMapper::kmutau));
      found=true;
    }
    if(triggerlist[iHLT]=="HLT_IsoMu24_eta2p1_v1"){
      triggerMap.push_back(triggerMapper(triggerlist[iHLT]
				      ,"hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09"
				      ,"",triggerMapper::kmutau));
      found=true;
    }
    if(triggerlist[iHLT]=="HLT_IsoMu27_v1"){
      triggerMap.push_back(triggerMapper(triggerlist[iHLT]
				      ,"hltL3crIsoL1sMu25L1f0L2f10QL3f27QL3trkIsoFiltered0p09"
				      ,"",triggerMapper::kmutau));
      found=true;
    }
    if(triggerlist[iHLT]=="HLT_Ele22_eta2p1_WP75_Gsf_LooseIsoPFTau20_v1"){
      TString listfilt2[2]={"hltPFTau20TrackLooseIso ","hltOverlapFilterIsoEle22WP75GsfLooseIsoPFTau20"};
      TString listfilt1[2]={"hltEle22WP75L1IsoEG20erTau20erGsfTrackIsoFilter ","hltOverlapFilterIsoEle22WP75GsfLooseIsoPFTau20"};
      triggerMap.push_back(triggerMapper(triggerlist[iHLT]
				      ,listfilt1
				      ,listfilt2,2,2,triggerMapper::ketau));
      found=true;
    }
    if(triggerlist[iHLT]=="HLT_Ele32_eta2p1_WP75_Gsf_v1"){
      triggerMap.push_back(triggerMapper(triggerlist[iHLT]
				      ,"hltEle32WP75GsfTrackIsoFilter"
				      ,"",triggerMapper::ketau));
      found=true;
    }
    if(triggerlist[iHLT]=="HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v1"){
      triggerMap.push_back(triggerMapper(triggerlist[iHLT]
				      ,"hltDoublePFTau40TrackPt1MediumIsolationDz02Reg"
				      ,"hltDoublePFTau40TrackPt1MediumIsolationDz02Reg",triggerMapper::ktautau));
      found=true;
    }
    if(!found)triggerMap.push_back(triggerMapper (triggerlist[iHLT] , dummyList, dummyList, 1, 1, triggerMapper::kemu));
  }


  /*
  triggerMap[0]=triggerMapper("HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v1"
                      ,"hltMu23Ele12GsfTrackIsoLegEle12GsfCaloIdTrackIdIsoMediumWPFilter",
                      "hltL1Mu12EG7L3IsoMuFiltered23",triggerMapper::kemu);
                      
  triggerMap[1]=triggerMapper("HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v1"
                      ,"hltMu8Ele23GsfTrackIsoLegEle23GsfCaloIdTrackIdIsoMediumWPFilter",
                      "hltL1sL1Mu5EG20ORL1Mu5IsoEG18L3IsoFiltered8",triggerMapper::kemu);
  
  TString listfilt2[2]={"hltL1sMu16erTauJet20er","hltOverlapFilterIsoMu17LooseIsoPFTau20"};
  TString listfilt1[1]={"hltOverlapFilterIsoMu17LooseIsoPFTau20"};
  triggerMap[2]=triggerMapper("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v1"
                      ,listfilt1
                      ,listfilt2,1,2,triggerMapper::kmutau);
  
  triggerMap[3]=triggerMapper("HLT_IsoMu24_eta2p1_IterTrk02_v1"
                      ,"hltL3crIsoL1sMu20Eta2p1L1f0L2f20QL3f24QL3crIsoRhoFiltered0p15IterTrk02",
                      "",triggerMapper::kmutau);
  
  listfilt2[0]="hltL1sL1IsoEG20erTauJet20er ";
  listfilt2[1]="hltOverlapFilterIsoEle22WP85GsfLooseIsoPFTau20";
  listfilt1[0]="hltOverlapFilterIsoEle22WP85GsfLooseIsoPFTau20";
  triggerMap[4]=triggerMapper("HLT_Ele22_eta2p1_WP85_Gsf_LooseIsoPFTau20_v1"
                      ,listfilt1,
                      listfilt2,1,2,triggerMapper::ketau);
  
  triggerMap[5]=triggerMapper("HLT_Ele27_eta2p1_WP85_Gsf_v1"
                      ,"hltEle27WP85GsfTrackIsoFilter",
                      "",triggerMapper::ketau);
  
  TString taulist[3]={"hltL1sDoubleTauJet36erORDoubleTauJet68er", "hltDoubleL2IsoTau35eta2p1", "hltDoublePFTau40TrackPt1MediumIsolationDz02Reg"};
  triggerMap[6]=triggerMapper("HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v1"
                      ,taulist,
                      taulist,3,3,triggerMapper::ktautau);
  */

  /*
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
  };*/
  
  /*
  for(int i=0;i<nTriggers;i++){
    triggerlist[i]=triggerMap[i].GetHLTPath();
    //pathKind[i]=-1;
    //triggerlist[i].Prepend("HLT_");
    //triggerlist[i].Append("_v1");
  }
  */

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
/*    
  pathKind[0]=kemuPath;
  pathKind[1]=kemuPath;
  pathKind[2]=kmutauPath;
  pathKind[3]=kmutauPath;
  pathKind[4]=ketauPath;
  pathKind[5]=ketauPath;
  pathKind[6]=ktautauPath;

  //Careful: e<mu<tau
  nFiltersforPath_leg1={1,1,1,1,1,1,3};
  nFiltersforPath_leg2={1,1,2,0,2,0,0};
  filterPaths_leg1={
    "hltMu23Ele12GsfTrackIsoLegEle12GsfCaloIdTrackIdIsoMediumWPFilter",
    "hltMu8Ele23GsfTrackIsoLegEle23GsfCaloIdTrackIdIsoMediumWPFilter",
    "hltOverlapFilterIsoMu17LooseIsoPFTau20",
    "hltL3crIsoL1sMu20Eta2p1L1f0L2f20QL3f24QL3crIsoRhoFiltered0p15IterTrk02",
    "hltOverlapFilterIsoEle22WP85GsfLooseIsoPFTau20",
    "hltEle27WP85GsfTrackIsoFilter",
    "hltL1sDoubleTauJet36erORDoubleTauJet68er ",
    "hltDoubleL2IsoTau35eta2p1",
    "hltDoublePFTau40TrackPt1MediumIsolationDz02Reg"
  };
  filterPaths_leg2={
    "hltL1Mu12EG7L3IsoMuFiltered23",
    "hltL1sL1Mu5EG20ORL1Mu5IsoEG18L3IsoFiltered8",
    "hltL1sMu16erTauJet20er ",
    "hltOverlapFilterIsoMu17LooseIsoPFTau20",
    "hltL1sL1IsoEG20erTauJet20er ",
    "hltOverlapFilterIsoEle22WP85GsfLooseIsoPFTau20",
  };*/
}

triggerhelper::~triggerhelper(){
  delete triggerlist;
  //delete triggerMap;
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
  event.getByLabel(InputTag("TriggerResults","","HLT"), metFilterBits);
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

triggerMapper triggerhelper::GetTriggerMap(TString path){
  for(int i=0;i<nTriggers;i++){
    if(triggerMap[i].GetHLTPath().CompareTo(path.Data())==0)return triggerMap[i];
  }
  return triggerMapper();
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

TString triggerhelper::printTriggerName(int ntrigger){
  if(ntrigger<nTriggers) return triggerlist[ntrigger];
  return "";
}

