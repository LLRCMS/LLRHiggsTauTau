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


//#include "FWCore/Framework/interface/Event.h"

#include <TString.h>
#include <string>
#include <vector>
#include <sstream>
#include <utility>
#include <iostream> //FRA

using namespace std;

class triggerMapper {
 public:
  enum triggerChannel {kemu=0, ketau=1,kmutau=2,ktautau=3, kmumu=4, kee=5};
  enum legIDtype {kele=11, kmu=13, ktau=15, kOther=999}; 
  triggerMapper();
  triggerMapper(const triggerMapper& trigMap);
  triggerMapper(string, std::vector<std::string>,std::vector<std::string>, int);
  triggerMapper(string, std::vector<std::string>,std::vector<std::string>, int theleg1ID, int theleg2ID);
  triggerMapper(string, std::vector<std::string>,std::vector<std::string>, std::vector<std::string>, std::vector<std::string>, int theleg1ID, int theleg2ID); //FRA
  //triggerMapper(string, std::vector<std::string>,std::vector<std::string>, std::vector<std::string>, std::vector<std::string>, int theleg1ID, int theleg2ID, double, double); //FRA
  //triggerMapper(TString, TString*,TString*,int,int);
  triggerMapper(string, string, string, int);
  
  string GetHLTPath(){return HLT;}
  int GetNfiltersleg1(){return filter_leg1.size();}
  int GetNfiltersleg2(){return filter_leg2.size();}
  int GetNfiltersleg3(){return filter_leg3.size();} //FRA
  int GetNfiltersleg4(){return filter_leg4.size();} //FRA
  string Getfilter(bool isleg1, int iFilter);
  string GetfilterVBF(bool isleg3, int iFilter); //FRA
  //double GetPtCut1() {return pt1;} //FRA
  //double GetPtCut2() {return pt2;} //FRA
  std::vector<std::string> Getfilters(bool isleg1){return (isleg1 ? filter_leg1 : filter_leg2); };
  void Setfilters(std::vector<std::string>);
  int GetTriggerChannel(){return channel;}
  pair <int, int> GetTriggerLegsID(){return make_pair(leg1ID, leg2ID);}
  int GetLegFromID(int ID); //0: not fourd; 1 = leg1; 2 = leg2

  ~triggerMapper();

 private:
  string HLT;
  std::vector<string> filter_leg1;
  std::vector<string> filter_leg2;
  std::vector<string> filter_leg3; //FRA
  std::vector<string> filter_leg4; //FRA
  int channel; // final decay channel: mu tau, ele tau, ...
  int leg1ID;  // abs(pdgID) of leg 1
  int leg2ID;  // abs(pdgID) of leg 2
  //double pt1;  // pT cut for leg1 //FRA
  //double pt2;  // pT cut for leg2 //FRA
};
#endif
