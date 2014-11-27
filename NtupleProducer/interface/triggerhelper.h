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

#include <TString.h>
#include <string>
#include <vector>
#include <sstream>

using namespace std;

class triggerhelper {
 public:
  triggerhelper();
  int FindTriggerBit(const edm::Event&, const vector<string>, const vector<int>);
  bool IsTriggerFired(int triggerbit){return 0;}
  int printFiredPaths(){return 0;}
  bool checkIfPathIsFired(TString path){return 0;}
  bool checkIfPathIsFired(int path){return 0;}

  ~triggerhelper(){}

 private:
  static const int nTriggers =19;
  TString triggerlist[nTriggers];

};
#endif
