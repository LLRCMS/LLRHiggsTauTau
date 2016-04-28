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
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "Flag_CSCTightHalo2015Filter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_goodVertices",
    "Flag_eeBadScFilter"
  };
  for(int i=0;i<nMETs;i++)metlist[i]=tmpMETfilters[i];

  
}

triggerhelper::triggerhelper(TH1F* hCounter){

  for(int itr=1;itr<=hCounter->GetNbinsX();itr++){
    TString binlabel = hCounter->GetXaxis()->GetBinLabel(itr);
    if(binlabel.BeginsWith("HLT"))triggerlist.push_back(hCounter->GetXaxis()->GetBinLabel(itr));
  }

  string tmpMETfilters[nMETs]={
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "Flag_CSCTightHalo2015Filter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_goodVertices",
    "Flag_eeBadScFilter"
  };
  for(int i=0;i<nMETs;i++)metlist[i]=tmpMETfilters[i];


}

triggerhelper::triggerhelper()//:nTriggers(0)
{
  string tmpMETfilters[nMETs]={
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "Flag_CSCTightHalo2015Filter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_goodVertices",
    "Flag_eeBadScFilter"
  };
  for(int i=0;i<nMETs;i++)metlist[i]=tmpMETfilters[i];

}

triggerhelper::~triggerhelper(){
  //delete triggerlist;
  //delete triggerMap;
}

void triggerhelper::addTriggerMap(string hlt, vector<string> path1, vector<string> path2, int channel){
  
  if (find(triggerlist.begin(), triggerlist.end(), hlt) != triggerlist.end() )
  {
    cout << "** triggerHelper :: Warning: path " << hlt << " already added, skipping" << endl;
    return;
  }

  triggerlist.push_back(hlt);
  //const int n1 = path1.size();
  //const int n2 = path2.size();
  //triggerMapper map(hlt,path1,path2,n1,n2,channel);
  triggerMapper map(hlt,path1,path2,channel);
  triggerMap.push_back(map);
  //nTriggers++;
}

void triggerhelper::addTriggerMap(string hlt,vector<string> path1, vector<string> path2, int leg1ID, int leg2ID)
{
  if (find(triggerlist.begin(), triggerlist.end(), hlt) != triggerlist.end() )
  {
    cout << "** triggerHelper :: Warning: path " << hlt << " already added, skipping" << endl;
    return;
  }

  triggerlist.push_back(hlt);
  triggerMapper map(hlt,path1,path2,leg1ID, leg2ID);
  triggerMap.push_back(map);
}

Long64_t triggerhelper::FindTriggerBit(const edm::Event& event, const vector<string> foundPaths, const vector<int> indexOfPaths){
  
  Long64_t bit =0;
  
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
    for(int j=0;j<(int)foundPaths.size();j++){
      //string myString (triggerlist.at(it));
      //if(myString.compare(foundPaths.at(j))==0) // full trigger name match required
      
      string toCheckTrigger  = triggerlist.at(it) ;
      string elemAllTriggers = foundPaths.at(j) ;
      if (elemAllTriggers.find(toCheckTrigger) != std::string::npos) // equivalent to wildcard at the end or beginning of triggername 
      {
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
  bool error = false;
  edm::Handle<edm::TriggerResults> metFilterBits;
  event.getByLabel(InputTag("TriggerResults","","HLT"), metFilterBits);
  //edm::EDGetTokenT<edm::TriggerResults> metFilterBitsToken_;//(consumes<edm::TriggerResults>(<edm::InputTag>("metFilterBits")));
  //event.getByToken(metFilterBitsToken_, metFilterBits);
  const edm::TriggerNames &metNames = event.triggerNames(*metFilterBits);
  for(int im=0;im<nMETs;im++){
    bool foundFilter = false;
    for(unsigned int i = 0; i < metFilterBits->size(); ++i){      
      if(metlist[im].compare(metNames.triggerName(i))==0){
	foundFilter = true;
	if ( metFilterBits->accept(i) ) bit |= 1 <<im;
	break;
      }
    }
    if ( !foundFilter ) {
      cout << "** triggerHelper :: Failed to find MET filter " << metlist[im] << endl;
      error = true;
    }
  }
  if ( error ) bit *= -1;
  return bit;
}

triggerMapper triggerhelper::GetTriggerMap(string path){
  for(int i=0;i<(int)triggerMap.size();i++){
    //if(triggerMap.at(i).GetHLTPath().compare(path)==0)return triggerMap.at(i); // full name comparison
    if (path.find(triggerMap.at(i).GetHLTPath()) != std::string::npos) return triggerMap.at(i); // use av equivalent of "contains" for versioning wildcard
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
    //if(triggerlist.at(it).compare(triggername)==0) return it; // full nale matching
    if (triggername.find(triggerlist.at(it)) != std::string::npos) return it; // just check that the input name contains the wanted name
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

