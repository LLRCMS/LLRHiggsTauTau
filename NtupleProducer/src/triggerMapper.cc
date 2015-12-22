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
#include <vector>
#include <string>

using namespace std;
using namespace edm;

triggerMapper::triggerMapper(){
  HLT="";
  filter_leg1.push_back("");
  filter_leg2.push_back("");
}

triggerMapper::triggerMapper(string HLTtrigger, std::vector<string> filters1, std::vector<string> filters2, int c){
  HLT=HLTtrigger;
  int n1 = filter_leg1.size();
  int n2 = filter_leg2.size();
  for(int i=0;i<n1;i++)filter_leg1.push_back(filters1.at(i));
  for(int i=0;i<n2;i++)filter_leg2.push_back(filters2.at(i));
  
  channel = c;  
  if (channel == (int) kemu)
  {
    leg1ID = (int) kele;
    leg2ID = (int) kmu;
  }
  else if (channel == (int) ketau)
  {
    leg1ID = (int) kele;
    leg2ID = (int) ktau;
  }
  else if (channel == (int) kmutau)
  {
    leg1ID = (int) kmu;
    leg2ID = (int) ktau;
  }
  else if (channel == (int) ktautau)
  {
    leg1ID = (int) ktau;
    leg2ID = (int) ktau;
  }
  else
  {
    leg1ID = (int) kOther;
    leg2ID = (int) kOther;
  }
}

triggerMapper::triggerMapper(string HLTtrigger, std::vector<std::string> filters1, std::vector<std::string> filters2, int theleg1ID, int theleg2ID)
{
  HLT=HLTtrigger;
  int n1 = filter_leg1.size();
  int n2 = filter_leg2.size();
  for(int i=0;i<n1;i++)filter_leg1.push_back(filters1.at(i));
  for(int i=0;i<n2;i++)filter_leg2.push_back(filters2.at(i));

  leg1ID = theleg1ID;
  leg2ID = theleg2ID;

  if (leg1ID == (int) kele && leg2ID == (int) kmu)  channel = (int) kemu;
  else if (leg1ID == (int) kele && leg2ID == (int) ktau) channel = (int) ketau;
  else if (leg1ID == (int) kmu  && leg2ID == (int) ktau) channel = (int) kmutau;
  else if (leg1ID == (int) ktau && leg2ID == (int) ktau) channel = (int) ktautau;
  else channel = 999;

}


triggerMapper::triggerMapper(string HLTtrigger, string filter1, string filter2, int c){
  HLT=HLTtrigger;
  filter_leg1.push_back(filter1);
  filter_leg2.push_back(filter2);
  channel=c;
  if (channel == (int) kemu)
  {
    leg1ID = (int) kele;
    leg2ID = (int) kmu;
  }
  else if (channel == (int) ketau)
  {
    leg1ID = (int) kele;
    leg2ID = (int) ktau;
  }
  else if (channel == (int) kmutau)
  {
    leg1ID = (int) kmu;
    leg2ID = (int) ktau;
  }
  else if (channel == (int) ktautau)
  {
    leg1ID = (int) ktau;
    leg2ID = (int) ktau;
  }
  else
  {
    leg1ID = (int) kOther;
    leg2ID = (int) kOther;
  }
}

triggerMapper::~triggerMapper(){
  //HLT="";//delete HLT;
  filter_leg1.clear();
  filter_leg2.clear();
}

string triggerMapper::Getfilter(bool isleg1,int iFilter){
  string result("");
  if(isleg1){
    if(iFilter>(int)(filter_leg1.size()))cout<<"TRIGGER MAPPER ERROR NFILTER "<<iFilter<<endl;
    else result=filter_leg1.at(iFilter);
  }else {
    if(iFilter>(int)(filter_leg2.size()))cout<<"TRIGGER MAPPER ERROR NFILTER "<<iFilter<<endl;
    else result=filter_leg2.at(iFilter);
  }
  return result;
}

int triggerMapper::GetLegFromID(int ID)
{
  if (ID == leg1ID) return 1;
  if (ID == leg2ID) return 2;
  return 0;
}
