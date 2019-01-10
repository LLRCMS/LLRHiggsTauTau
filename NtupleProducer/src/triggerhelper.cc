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
    /* //FRA: 2016 data 80X
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "Flag_CSCTightHalo2015Filter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_goodVertices",
    "Flag_eeBadScFilter"*/
    
    //FRA: Fall17 94X
    "Flag_goodVertices",
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_globalSuperTightHalo2016Filter",
    "Flag_BadPFMuonFilter",
    "Flag_BadChargedCandidateFilter",
    "Flag_eeBadScFilter",
    "Flag_ecalBadCalibFilter"
  };
  for(int i=0;i<nMETs;i++)metlist[i]=tmpMETfilters[i];

  
}

triggerhelper::triggerhelper(TH1F* hCounter){

  for(int itr=1;itr<=hCounter->GetNbinsX();itr++){
    TString binlabel = hCounter->GetXaxis()->GetBinLabel(itr);
    if(binlabel.BeginsWith("HLT"))triggerlist.push_back(hCounter->GetXaxis()->GetBinLabel(itr));
  }

  string tmpMETfilters[nMETs]={
    /* //FRA: 2016 data 80X
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "Flag_CSCTightHalo2015Filter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_goodVertices",
    "Flag_eeBadScFilter"*/
    
    //FRA: Fall17 94X
    "Flag_goodVertices",
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_globalTightHalo2016Filter",
    "Flag_BadPFMuonFilter",
    "Flag_BadChargedCandidateFilter",
    "Flag_eeBadScFilter",
    "Flag_ecalBadCalibFilter"
  };
  for(int i=0;i<nMETs;i++)metlist[i]=tmpMETfilters[i];


}

triggerhelper::triggerhelper()//:nTriggers(0)
{
  string tmpMETfilters[nMETs]={
    /* //FRA: 2016 data 80X
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "Flag_CSCTightHalo2015Filter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_goodVertices",
    "Flag_eeBadScFilter"*/
    
    //FRA: Fall17 94X
    "Flag_goodVertices",
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_globalTightHalo2016Filter",
    "Flag_BadPFMuonFilter",
    "Flag_BadChargedCandidateFilter",
    "Flag_eeBadScFilter",
    "Flag_ecalBadCalibFilter"
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

//FRA : added path3 and path4 for VBF matching with jets
void triggerhelper::addTriggerMap(string hlt,vector<string> path1, vector<string> path2, vector<string> path3, vector<string> path4, int leg1ID, int leg2ID)
{

  if (find(triggerlist.begin(), triggerlist.end(), hlt) != triggerlist.end() )
  {
    cout << "** triggerHelper :: Warning: path " << hlt << " already added, skipping" << endl;
    return;
  }

  triggerlist.push_back(hlt);
  triggerMapper map(hlt, path1, path2, path3, path4, leg1ID, leg2ID);
  triggerMap.push_back(map);
}

//FRA : added pt cuts for leg1 and leg2
void triggerhelper::addTriggerMap(string hlt,vector<string> path1, vector<string> path2, vector<string> path3, vector<string> path4, int leg1ID, int leg2ID, double pt1, double pt2)
{

  if (find(triggerlist.begin(), triggerlist.end(), hlt) != triggerlist.end() )
  {
    cout << "** triggerHelper :: Warning: path " << hlt << " already added, skipping" << endl;
    return;
  }
  
  triggerlist.push_back(hlt);
  triggerMapper map(hlt, path1, path2, path3, path4, leg1ID, leg2ID, pt1, pt2);
  triggerMap.push_back(map);
}

Long64_t triggerhelper::FindTriggerBit(const edm::Event& event, const vector<string> foundPaths, const vector<int> indexOfPaths, const edm::Handle<edm::TriggerResults>& triggerResults){
  
  Long64_t bit =0;
  
  //edm::Handle<edm::TriggerResults> triggerResults;
  //event.getByLabel(InputTag("TriggerResults","","HLT"), triggerResults);
      
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
	if(triggerResults->accept(indexOfPaths[j]))bit |= long(1) <<it;
	break;
      }
    }
  }
  //printf("bit: %d\n",bit);
  return bit;
}

int triggerhelper::FindMETBit(const edm::Event& event, edm::EDGetTokenT<edm::TriggerResults> metFilterBitsToken){

  int bit =0;
  edm::Handle<edm::TriggerResults> metFilterBits;
  event.getByToken(metFilterBitsToken, metFilterBits); // PAT, RECO, HLT, SIM ?
  //edm::EDGetTokenT<edm::TriggerResults> metFilterBitsToken_;//(consumes<edm::TriggerResults>(<edm::InputTag>("metFilterBits")));
  //event.getByToken(metFilterBitsToken_, metFilterBits);
  const edm::TriggerNames &metNames = event.triggerNames(*metFilterBits);
  for(int im=0;im<nMETs;im++){
    // cout << "DEB looking for " << metlist[im] << endl;
    bool foundFilter = false;
    for(unsigned int i = 0; i < metFilterBits->size(); ++i){      
      // cout << " --> and testing " << metNames.triggerName(i) << endl;
      if(metlist[im].compare(metNames.triggerName(i))==0){
        foundFilter = true;
        if ( metFilterBits->accept(i) ) bit |= 1 <<im;
        break;
      }
    }
    if ( !foundFilter ) {
      cout << "** triggerHelper :: Failed to find MET filter " << metlist[im] << endl;
    }
  }
  return bit;
}

triggerMapper triggerhelper::GetTriggerMap(string path){

  for(int i=0;i<(int)triggerMap.size();i++){
    //if(triggerMap.at(i).GetHLTPath().compare(path)==0)return triggerMap.at(i); // full name comparison
    if (path.find(triggerMap.at(i).GetHLTPath()) != std::string::npos)
    {
      return triggerMap.at(i); // use av equivalent of "contains" for versioning wildcard
    }
  }
  cout << "** trigger mapper : error : path " << path << " not found in triggerMap stored , return empty map" << endl;
  return triggerMapper();
}

triggerMapper triggerhelper::GetTriggerMap(int idx){
  if (idx < (int)triggerMap.size()) return triggerMap.at(idx);
  cout << "** trigger mapper : error : index idx exceeds size of triggerMap stored = " << triggerMap.size() << " , return empty map" << endl;
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


bool triggerhelper::IsTriggerFired(Long64_t triggerbit, int triggernumber, bool isTrigger){ 
  int nLoop = triggerlist.size();
  Long64_t bitdigit = 1;
  if(!isTrigger)nLoop = nMETs;
  if(triggernumber>=0 && triggernumber<nLoop) return triggerbit & (bitdigit << triggernumber);
  return false;
}


int triggerhelper::printFiredPaths(Long64_t triggerbit, bool isTrigger){
  if (isTrigger) return printFiredPathsTrig (triggerbit);
  else return printFiredPathsMET (triggerbit);
}

int triggerhelper::printFiredPathsMET(Long64_t triggerbit)
{
  int nFired = 0;
  for (unsigned int it = 0; it < nMETs; it++)
  {
      printf("%s\n",metlist[it].c_str());
      nFired++;    
  }
  return nFired;
}

int triggerhelper::printFiredPathsTrig(Long64_t triggerbit)
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

