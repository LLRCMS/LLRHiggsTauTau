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
  triggerMapper(string, std::vector<std::string>,std::vector<std::string>, int, int, int);
  //triggerMapper(TString, TString*,TString*,int,int);
  triggerMapper(string, string, string, int);
  
  string GetHLTPath(){return HLT;}
  int GetNfiltersleg1(){return filter_leg1.size();}
  int GetNfiltersleg2(){return filter_leg2.size();}
  string Getfilter(bool isleg1, int iFilter);
  int GetTriggerChannel(){return channel;}
  
  ~triggerMapper();

 private:
  string HLT;
  std::vector<string> filter_leg1;
  std::vector<string> filter_leg2;
  int channel;
};
#endif
