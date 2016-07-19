/* 
**
** Helpers for gen info
** 
** 
** \date:    13 May 2015
** \author:  L. Cadamuro (LLR)
*/

#ifndef GenHelper_h
#define GenHelper_h

#include "TVector3.h"
#include "TLorentzVector.h"

#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include <vector>

namespace genhelper {

    enum HZDecay {
        MuHad  = 0,
        EHad   = 1,
        HadHad = 2,
        MuMu   = 3,
        EE     = 4,
        EMu    = 5,
        EEPrompt = 6, // prompt Z->ee/mumu decays
        MuMuPrompt = 7,
        Other  = 8 // for e.g. h->bb
    };

    enum WDecay {
      Had      = 0, // W->qqbar
      MuPrompt = 1,
      EPrompt  = 2,
      TauMu    = 3, // W->tau->mu
      TauE     = 4, // W->tau->e
      TauHad   = 5, // W->tau->tauh
      other  = 6
    };

    bool IsLastCopy (const reco::GenParticle& part); // return true if particle has no sons with its same pdgId to reject showering clones
    bool IsFirstCopy (const reco::GenParticle& part, const bool checkAbsPdg = false); // return true if particle has no mothers with its same pdgId to handle showering clones
    
    int GetTauDecay (const reco::GenParticle& part); // 0: tau->mu; 1: tau->ele; 2: tau->had
    int GetTauDecay (const reco::Candidate* part); // 0: tau->mu; 1: tau->ele; 2: tau->had

    const reco::Candidate* GetFirstCopy (const reco::Candidate* part); // follow all the replicated particle chain until the first clone
    const reco::Candidate* GetLastCopy (const reco::Candidate* part); // follow all the replicated particle chain until the last clone
    HZDecay GetHZDecay (const reco::Candidate* part); // return final state for H/Z -> see enum for code
    WDecay GetWDecay (const reco::Candidate* part); // return final state for W -> see enum for code
    WDecay GetTopDecay (const reco::Candidate* part); // return final state for top (= final state for W) -> see enum for code    

    reco::GenParticle GetTauHad (const reco::Candidate* part); // build had tau by summing sons without nu
    reco::GenParticle GetTauHadNeutrals (const reco::Candidate* part); // build neutral component of had tau by summing sons without nu
    
    const reco::Candidate* IsFromID (const reco::Candidate* part, int targetPDGId); // find if is son of a certain particle (select by targetPDGId); if not found, return NULL, else return its pointer
    int GetIndexInOutput (const reco::Candidate* part, std::vector<const reco::Candidate *> cands);

    typedef reco::GenParticleCollection::const_iterator IG;
    typedef reco::GenParticleRefVector::const_iterator IGR;
    TVector3 ImpactParameter(const TVector3& pv, const TVector3& sv, const TLorentzVector& p4);//Calculate generator level impact parameter
    void GetTausDaughters(const reco::GenParticle& tau, reco::GenParticleRefVector& products, bool ignoreNus, bool direct);
    void FindDescendents(const reco::GenParticle& base, reco::GenParticleRefVector& descendents, int status, int pdgId=0, bool skipPhotonsPi0AndFSR=false);
    const reco::GenParticleRef GetLeadChParticle(const reco::GenParticleRefVector& products);
    int getDetailedTauDecayMode(const reco::GenParticleRefVector& products);
    
}
#endif
