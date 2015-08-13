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

triggerhelper::triggerhelper(vector<string> HLTPaths) //: nTriggers(HLTPaths.size())
{

  //cout << "nTriggers: " << nTriggers << endl;
  triggerlist=HLTPaths;
  string tmpMETfilters[nMETs]={
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
         "Flag_trkPOG_toomanystripclus53X"  
  };
  for(int i=0;i<nMETs;i++)metlist[i]=tmpMETfilters[i];

  
}

triggerhelper::triggerhelper(TH1F* hCounter){

  for(int itr=1;itr<=hCounter->GetNbinsX();itr++){
    TString binlabel = hCounter->GetXaxis()->GetBinLabel(itr);
    if(binlabel.BeginsWith("HLT"))triggerlist.push_back(hCounter->GetXaxis()->GetBinLabel(itr));
  }

  string tmpMETfilters[nMETs]={
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
         "Flag_trkPOG_toomanystripclus53X"
  };
  for(int i=0;i<nMETs;i++)metlist[i]=tmpMETfilters[i];


}

triggerhelper::triggerhelper()//:nTriggers(0)
{
  string tmpMETfilters[nMETs]={
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
         "Flag_trkPOG_toomanystripclus53X"  
  };
  for(int i=0;i<nMETs;i++)metlist[i]=tmpMETfilters[i];

}

triggerhelper::~triggerhelper(){
  //delete triggerlist;
  //delete triggerMap;
}

void triggerhelper::addTriggerMap(string hlt, vector<string> path1, vector<string> path2, int channel){
  triggerlist.push_back(hlt);
  const int n1 = path1.size();
  const int n2 = path2.size();
  triggerMapper map(hlt,path1,path2,n1,n2,channel);
  triggerMap.push_back(map);
  //nTriggers++;
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
 
  for(int it=0;it<(int)triggerlist.size();it++){
    //for(int j=0;j<(int)foundPaths.size();j++){
    for(int j=0;j<(int)foundPaths.size();j++){
      string myString (triggerlist.at(it));
      if(myString.compare(foundPaths.at(j))==0){
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
      if(metlist[im].compare(metNames.triggerName(i))==0){
	      if(metFilterBits->accept(i))bit |= 1 <<im;
	      break;
      }
    }
  }
  return bit;
}

triggerMapper triggerhelper::GetTriggerMap(string path){
  for(int i=0;i<(int)triggerMap.size();i++){
    if(triggerMap.at(i).GetHLTPath().compare(path)==0)return triggerMap.at(i);
  }
  return triggerMapper();
}


int triggerhelper::FindTriggerNumber(string triggername, bool isTrigger){ 
  if (isTrigger) return FindTriggerNumberTrig (triggername);
  else return FindTriggerNumberMET (triggername);
}

int triggerhelper::FindTriggerNumberMET (string triggername)
{
  for(int it=0; it < nMETs; it++){   
    if(metlist[it].compare(triggername)==0) return it;
  }
  return -1;
}

int triggerhelper::FindTriggerNumberTrig (string triggername)
{
  for(unsigned int it=0; it < triggerlist.size(); it++){   
    if(triggerlist.at(it).compare(triggername)==0) return it;
  }
  return -1;
}


bool triggerhelper::IsTriggerFired(int triggerbit, int triggernumber, bool isTrigger){ 
  int nLoop = triggerlist.size();
  if(!isTrigger)nLoop = nMETs;
  if(triggernumber>=0 && triggernumber<nLoop) return triggerbit & (1 << triggernumber);
  return false;
}


int triggerhelper::printFiredPaths(int triggerbit, bool isTrigger){
  if (isTrigger) return printFiredPathsTrig (triggerbit);
  else return printFiredPathsMET (triggerbit);
}

int triggerhelper::printFiredPathsMET(int triggerbit)
{
  int nFired = 0;
  for (unsigned int it = 0; it < nMETs; it++)
  {
      printf("%s\n",metlist[it].c_str());
      nFired++;    
  }
  return nFired;
}

int triggerhelper::printFiredPathsTrig(int triggerbit)
{
  int nFired = 0;
  for (unsigned int it = 0; it < triggerlist.size(); it++)
  {
      printf("%s\n",triggerlist.at(it).c_str());
      nFired++;    
  }  
  return nFired;
}


string triggerhelper::printTriggerName(int ntrigger){
  if(ntrigger<(int)triggerlist.size()) return triggerlist.at(ntrigger);
  string dummy ("");
  return dummy;
}

