#ifndef TriggerHelper_h
#define TriggerHelper_h

/** \class triggerhelper
 *
 *  No description available.
 *
 *  $Date: 2012/06/06 00:27:43 $
 *  $Revision: 1.3 $
 *  \author G. Ortona - LLR
 */


#include "FWCore/Framework/interface/Event.h"
#include "LLRHiggsTauTau/NtupleProducer/interface/triggerMapper.h"

#include <TString.h>
#include <string>
#include <vector>
#include <sstream>

using namespace std;

class triggerhelper {
 public:
  triggerhelper();
  int FindTriggerBit(const edm::Event&, const vector<string>, const vector<int>);
  int FindMETBit(const edm::Event&);
  int FindTriggerNumber(TString triggername,bool istrigger=true);
  bool IsTriggerFired(int triggerbit, int triggerNumber,bool istrigger=true);
  bool IsTriggerFired(int triggerbit, TString triggerName,bool istrigger=true){return IsTriggerFired(triggerbit, FindTriggerNumber(triggerName));}
  int printFiredPaths(int triggerbit,bool istrigger=true);
  triggerMapper GetTriggerMap(TString trigger);
  int GetNTriggers(){return nTriggers;}
  TString printTriggerName(int ntrigger);

  ~triggerhelper(){}

 private:
  static const int nTriggers =7;
  TString triggerlist[nTriggers];
  triggerMapper triggerMap[nTriggers];
  static const int nMETs =13;
  TString metlist[nMETs];

};
#endif
