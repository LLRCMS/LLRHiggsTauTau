#ifndef LeptonIsoHelper_h
#define LeptonIsoHelper_h

/** \class LeptonIsoHelper
 *
 *  Helper for computing lepon isolation
 *
 *  $Date: 2012/06/10 17:40:50 $
 *  $Revision: 1.3 $
 *  \author N. Amapane
 */

#include <FWCore/Utilities/interface/InputTag.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Tau.h>
#include "DataFormats/PatCandidates/interface/Jet.h"
#include <algorithm>



namespace SelfVetoPolicy
{
   enum SelfVetoPolicy
     {
	selfVetoNone=0, selfVetoAll=1, selfVetoFirst=2
     };
}



namespace LeptonIsoHelper {

  /// Set the rho tag and the EA targed based on setup, for muons and electrons
  edm::InputTag getMuRhoTag(int sampleType, int setup);

  edm::InputTag getEleRhoTag(int sampleType, int setup);
  
  /// Compute combRelIsoPF for a mu
  float combRelIsoPF(int sampleType, int setup, double rho, const pat::Muon& mu, bool dr03=false, float fsr=0);
  
  /// Compute combRelIsoPF for a tau
  float combRelIsoPF(const pat::Tau& mu);

  /// Compute combRelIsoPF for an ele
  float combRelIsoPF(int sampleType, int setup, double rho, const pat::Electron& ele, float fsr=0);
  
  /// Generic version, assuming that Lep is a PATObject; calls one of the above
  float combRelIsoPF(int sampleType, int setup, double rho, const reco::Candidate* lep, float fsr=0);


  
  //MiniRel isolation
  void PFIso_particles(const edm::View<pat::PackedCandidate>* pfCands, std::vector<const pat::PackedCandidate *> & PFIso_charged, std::vector<const pat::PackedCandidate *> & PFIso_neutral);
  float isoSumRaw(const reco::Candidate* cand, const std::vector<const pat::PackedCandidate *> pfCands_Iso, float dR, float innerR, float threshold, SelfVetoPolicy::SelfVetoPolicy selfVeto, int pdgId=-1);
  float PfIsoCharged(const reco::Candidate* cand, const std::vector<const pat::PackedCandidate *> pfCands_charged, float miniIsoR);
  float PfIsoNeutral(const reco::Candidate* cand, const std::vector<const pat::PackedCandidate *> pfCands_neutral, float miniIsoR);
  std::pair<float,float> miniRelIso_ChargedNeutral(const reco::Candidate* cand, const std::vector<const pat::PackedCandidate *> pfCands_charged, const std::vector<const pat::PackedCandidate *> pfCands_neutral, float rho, int year);

  //Used for lepton MVA
  int jetNDauChargedMVASel(const reco::Candidate* cand, pat::Jet jet);
  float jetPtRel(const reco::Candidate& cand, const pat::Jet& jet, std::string JECname);
  float jetPtRatio(const reco::Candidate& cand, const pat::Jet& jet, std::string JECname);

}
#endif


