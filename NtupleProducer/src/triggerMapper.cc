/** \class triggerMapper
 *
 *  Maps tau trigger paths and filters to match to it
 *
 *  $Date: 2012/06/06 00:27:43 $
 *  $Revision: 1.3 $
 *  \author G. Ortona - LLR
 */
 
#include <DataFormats/Common/interface/TriggerResults.h>
#include <FWCore/Common/interface/TriggerNames.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/TriggerNamesService.h>
#include "LLRHiggsTauTau/NtupleProducer/interface/triggerMapper.h"
#include <iostream>

using namespace std;
using namespace edm;

triggerMapper::triggerMapper(){
  HLT="";
  filter_leg1.push_back("");
  filter_leg2.push_back("");
}

triggerMapper::triggerMapper(TString HLTtrigger, TString *filters1, TString *filters2,int n1, int n2, int c){
  HLT=HLTtrigger.Data();
  for(int i=0;i<n1;i++)filter_leg1.push_back(filters1[i].Data());
  for(int i=0;i<n1;i++)filter_leg2.push_back(filters2[i].Data());
  channel = c;
}

triggerMapper::triggerMapper(TString HLTtrigger, TString filter1, TString filter2, int c){
  HLT=HLTtrigger.Data();
  filter_leg1.push_back(filter1.Data());
  filter_leg2.push_back(filter2.Data());
  channel=c;
}

triggerMapper::~triggerMapper(){
  //HLT="";//delete HLT;
  filter_leg1.clear();
  filter_leg2.clear();
}

TString triggerMapper::Getfilter(bool isleg1,int iFilter){
  TString result="";
  if(isleg1){
    if(iFilter>(int)(filter_leg1.size()))cout<<"TRIGGER MAPPER ERROR NFILTER "<<iFilter<<endl;
    else result=filter_leg1.at(iFilter).Data();
  }else {
    if(iFilter>(int)(filter_leg2.size()))cout<<"TRIGGER MAPPER ERROR NFILTER "<<iFilter<<endl;
    else result=filter_leg2.at(iFilter).Data();
  }
  return result;
}