#ifndef TriggerMapper_h
#define TriggerMapper_h

/** \class triggerMapper
 *
 *  No description available.
 *
 *  $Date: 2012/06/06 00:27:43 $
 *  $Revision: 1.3 $
 *  \author G. Ortona - LLR
 */


#include "FWCore/Framework/interface/Event.h"

#include <TString.h>
#include <string>
#include <vector>
#include <sstream>

using namespace std;

class triggerMapper {
 public:
  enum triggerChannel {kemu=0, ketau=1,kmutau=2,ktautau=3};
  triggerMapper();
  triggerMapper(TString, TString*,TString*,int,int,int);
  //triggerMapper(TString, TString*,TString*,int,int);
  triggerMapper(TString, TString,TString,int);
  
  TString GetHLTPath(){return HLT.Data();}
  int GetNfiltersleg1(){return filter_leg1.size();}
  int GetNfiltersleg2(){return filter_leg2.size();}
  TString Getfilter(bool isleg1,int iFilter);
  int GetTriggerChannel(){return channel;}
  
  ~triggerMapper();

 private:
  TString HLT;
  std::vector<TString> filter_leg1;
  std::vector<TString> filter_leg2;
  int channel;
};
#endif
