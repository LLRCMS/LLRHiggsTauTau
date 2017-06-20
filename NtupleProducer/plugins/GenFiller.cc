/* \class GenFiller
**
** Create a collection of filtered gen level objects
** (including the creation of hadronic taus)
** 
** \date:    13 May 2015
** \author:  L. Cadamuro (LLR)
*/

#define DEBUG false

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Utilities/interface/InputTag.h>
#include <DataFormats/PatCandidates/interface/GenericParticle.h>
#include <DataFormats/PatCandidates/interface/PackedGenParticle.h>
#include <DataFormats/HepMCCandidate/interface/GenStatusFlags.h>
#include "LLRHiggsTauTau/NtupleProducer/interface/GenHelper.h"
//#include <LLRHiggsTauTau/NtupleProducer/interface/DaughterDataHelpers.h>

#include <vector>

using namespace edm;
using namespace std;
using namespace reco;

class GenFiller : public edm::EDProducer {
 public:
  /// Constructor
  explicit GenFiller(const edm::ParameterSet&);
  /// Destructor
  ~GenFiller(){};  

 private:
  virtual void beginJob(){};  
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};

  bool IsInteresting (const GenParticle& p);
  int makeFlagVector (const GenParticle* p); // put all gen flags in a single int word
  bool isVBFParton(const GenParticle& p);

  //edm::InputTag src_;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > src_;
  std::vector<const reco::Candidate *> cands_;
  //std::vector<reco::GenParticle> tauHadcands_; // gen H tau build in this class
  std::vector<int> tauHadcandsMothers_; // contains the index in the cands_ vector of the tauh mother
  const bool storeLightFlavAndGlu_;
};

// ------------------------------------------------------------------

GenFiller::GenFiller(const edm::ParameterSet& iConfig):
src_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("src"))),
storeLightFlavAndGlu_(iConfig.getParameter<bool>("storeLightFlavAndGlu"))
{
    //src_ = iConfig.getParameter<InputTag>("src");
    produces<pat::GenericParticleCollection>();
}

void GenFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    Handle <edm::View<reco::GenParticle> > genHandle;
    iEvent.getByToken (src_, genHandle);
    cands_.clear();
    //tauHadcands_.clear();
    tauHadcandsMothers_.clear();
    
    // output collection
    //auto_ptr<pat::GenericParticleCollection> result( new pat::GenericParticleCollection );
	std::unique_ptr<pat::GenericParticleCollection> result( new pat::GenericParticleCollection );
	
    unsigned int Ngen = genHandle->size();

    // fill vector of interesting Candidate* object, also used later for indexes
    // tau decays are analyzed here to build tauH candidates and put them in the list
    for (unsigned int iGen = 0; iGen < Ngen; iGen++)
    {
        const GenParticle& genP = (*genHandle)[iGen];
        if (IsInteresting (genP))
        {
            cands_.push_back (&genP); 
        }
    }    
    
    // loop on all previously filtered Gen Particles and establish relations between them + set flags
    unsigned int NGenSel = cands_.size();
    if (DEBUG) cout << "SELECTED PARTICLES: " << NGenSel << endl;
    for (unsigned int iGen = 0; iGen < NGenSel; iGen++)
    {        
        const reco::Candidate* genP = cands_.at(iGen);
        const GenParticle* genPClone = (GenParticle*) cands_.at(iGen);
        int PdgId = genP->pdgId();
        int APdgId = abs(PdgId);

        //only select some particles
        pat::GenericParticle filtGenP (*genP); // to be saved in output
        if (DEBUG) cout << iGen << " | id: " << genP->pdgId () << "   pt: " << genP->pt() << "   eta: " << genP->eta() << " | px: " << genP->px() << " , eta: " << genP->eta() << endl;
        
        // ------------------- general info flag on particles
        filtGenP.addUserInt ("generalGenFlags", makeFlagVector (genPClone));

        // ------------------- tau decay flags
        if (APdgId == 15)
        {
            int decay = genhelper::GetTauDecay(genP);
            filtGenP.addUserInt ("tauGenDecayMode", decay);
            if (decay == 2) tauHadcandsMothers_.push_back (iGen); // for later usage for tauh
            if (DEBUG) cout << "   --> tau decay: " << decay << endl;
        }

        // -------------------- H/Z decay mode: set final state and is is prompt
        if (APdgId == 25 || APdgId == 23 || APdgId == 36)
        {
            genhelper::HZDecay decay = genhelper::GetHZDecay (genP);
            filtGenP.addUserInt ("HZDecayMode", static_cast<int> (decay));
            if (DEBUG) cout << "   --> H/Z decay: " << decay << endl;
        }

        // -------------------- W decay mode: set final state and is is prompt
        if (APdgId == 24)
        {
            genhelper::WDecay decay = genhelper::GetWDecay (genP);
            filtGenP.addUserInt ("WDecayMode", static_cast<int> (decay));
            if (DEBUG) cout << "   --> W decay: " << decay << endl;
        }

        // -------------------- top decay mode: set final state and is is prompt
        if (APdgId == 6)
        {
            genhelper::WDecay decay = genhelper::GetTopDecay (genP);
            filtGenP.addUserInt ("TopDecayMode", static_cast<int> (decay));
            if (DEBUG) cout << "   --> Top decay: " << decay << endl;
        }

        // ------------------- additional flags for b (first & last) ==> TO DO if needed
        /*
        
        if (APdgId == 5)
        {
            bool IsFirst = 
            //filtGenP.addUserInt ("tauGenDecayMode", decay);
            //if (decay == 2) tauHadcandsMothers_.push_back (iGen); // for later usage for tauh
            //cout << "   --> tau decay: " << decay << endl;
        }
        */
        
        // ------------------- mother info flag
        
        const reco::Candidate* MothPtr;

        // H        
        MothPtr = genhelper::IsFromID (genP, 25);
        if (MothPtr != NULL) // save space, only add userfloats when valid
        {
            //filtGenP.addUserInt ("fromH", 1);
            filtGenP.addUserInt ("HMothIndex", genhelper::GetIndexInOutput(MothPtr, cands_));       
            if (DEBUG) cout << "   --> fromH: 1, indexH: " << genhelper::GetIndexInOutput(MothPtr, cands_) << endl;
        }

        // H MSSM
        MothPtr = genhelper::IsFromID (genP, 36);
        if (MothPtr != NULL) // save space, only add userfloats when valid
        {
            //filtGenP.addUserInt ("fromH", 1);
            filtGenP.addUserInt ("MSSMHMothIndex", genhelper::GetIndexInOutput(MothPtr, cands_));       
            if (DEBUG) cout << "   --> fromH(MSSM): 1, indexH: " << genhelper::GetIndexInOutput(MothPtr, cands_) << endl;
        }
	
        // Z        
        MothPtr = genhelper::IsFromID (genP, 23);
        if (MothPtr != NULL) // save space, only add userfloats when valid
        {
            //filtGenP.addUserInt ("fromZ", 1);
            filtGenP.addUserInt ("ZMothIndex", genhelper::GetIndexInOutput(MothPtr, cands_));       
            if (DEBUG) cout << "   --> fromZ: 1, indexZ: " << genhelper::GetIndexInOutput(MothPtr, cands_) << endl;
        }
        
        // top
        MothPtr = genhelper::IsFromID (genP, 6);
        if (MothPtr != NULL)
        {
            //filtGenP.addUserInt ("fromTop", 1);
            filtGenP.addUserInt ("TopMothIndex", genhelper::GetIndexInOutput(MothPtr, cands_));       
            if (DEBUG) cout << "   --> fromTop: 1, indexTop: " << genhelper::GetIndexInOutput(MothPtr, cands_) << endl;
        }

        // W        
        MothPtr = genhelper::IsFromID (genP, 24);
        if (MothPtr != NULL) // save space, only add userfloats when valid
        {
            //filtGenP.addUserInt ("fromW", 1);
            filtGenP.addUserInt ("WMothIndex", genhelper::GetIndexInOutput(MothPtr, cands_));       
            if (DEBUG) cout << "   --> fromW: 1, indexW: " << genhelper::GetIndexInOutput(MothPtr, cands_) << endl;
        }

        // b
        MothPtr = genhelper::IsFromID (genP, 5);
        if (MothPtr != NULL)
        {
            //filtGenP.addUserInt ("fromb", 1);
            filtGenP.addUserInt ("bMothIndex", genhelper::GetIndexInOutput(MothPtr, cands_));       
            if (DEBUG) cout << "   --> fromb: 1, indexb: " << genhelper::GetIndexInOutput(MothPtr, cands_) << endl;
        }

        // tau
        MothPtr = genhelper::IsFromID (genP, 15);
        if (MothPtr != NULL)
        {
            //filtGenP.addUserInt ("fromTau", 1);
            filtGenP.addUserInt ("TauMothIndex", genhelper::GetIndexInOutput(MothPtr, cands_));       
            if (DEBUG) cout << "   --> fromTau: 1, indexTau: " << genhelper::GetIndexInOutput(MothPtr, cands_) << endl;
        }


	if(filtGenP.hasUserInt("tauGenDecayMode") &&
	   (filtGenP.hasUserInt("ZMothIndex") || filtGenP.hasUserInt("HMothIndex") || filtGenP.hasUserInt("MSSMHMothIndex"))){

	  TVector3 aPVGenPoint = TVector3(genPClone->vx(), genPClone->vy(), genPClone->vz());
	  reco::GenParticleRefVector tauDaughters;
	  genhelper::GetTausDaughters(*genPClone,tauDaughters,true,false);
	  int detailedDecayMode = genhelper::getDetailedTauDecayMode(tauDaughters);
	  filtGenP.addUserInt("tauGenDetailedDecayMode", detailedDecayMode);
	  
	  reco::GenParticleRef leadChParticleRef = genhelper::GetLeadChParticle(tauDaughters);

	  TLorentzVector p4LeadingChParticle(leadChParticleRef->px(),
					     leadChParticleRef->py(),
					     leadChParticleRef->pz(),
					     leadChParticleRef->energy());
	  TVector3 tauDecayVertex(leadChParticleRef->vx(), leadChParticleRef->vy(), leadChParticleRef->vz());
	  TVector3 pca = genhelper::ImpactParameter(aPVGenPoint, tauDecayVertex, p4LeadingChParticle);
	  filtGenP.addUserFloat("pca_x",pca.X());
	  filtGenP.addUserFloat("pca_y",pca.Y());
	  filtGenP.addUserFloat("pca_z",pca.Z());
	}

        result->push_back (filtGenP);
    }


    // finally, do the hadronic tau (leave as last step not to spoil internal ordering of mother / daughter vector)    
    if (DEBUG) cout << "BUILDING tauH, size: " << tauHadcandsMothers_.size() << endl;
    for (unsigned int iTauH = 0; iTauH < tauHadcandsMothers_.size(); iTauH++)
    {
        int tauMothInd = tauHadcandsMothers_.at(iTauH);
        pat::GenericParticle tauH (genhelper::GetTauHad(cands_.at(tauMothInd)));
	pat::GenericParticle tauH_neutral (genhelper::GetTauHadNeutrals(cands_.at(tauMothInd)));
        tauH.addUserInt ("TauMothIndex", tauMothInd);
	tauH_neutral.addUserInt ("TauMothIndex", tauMothInd);
        
        // copy all the other flags from original tau
        pat::GenericParticle& tauMothGenP = result->at(tauMothInd);
        if (tauMothGenP.hasUserInt("HMothIndex") ){
	  tauH.addUserInt ("HMothIndex", tauMothGenP.userInt ("HMothIndex"));
	  tauH_neutral.addUserInt ("HMothIndex", tauMothGenP.userInt ("HMothIndex"));
	}
        if (tauMothGenP.hasUserInt("MSSMHMothIndex") ){
	  tauH.addUserInt ("MSSMHMothIndex", tauMothGenP.userInt ("MSSMHMothIndex"));
	  tauH_neutral.addUserInt ("MSSMHMothIndex", tauMothGenP.userInt ("MSSMHMothIndex"));
	}
        if (tauMothGenP.hasUserInt("TopMothIndex") ){
	  tauH.addUserInt ("TopMothIndex", tauMothGenP.userInt ("TopMothIndex"));
	  tauH_neutral.addUserInt ("TopMothIndex", tauMothGenP.userInt ("TopMothIndex"));
	}
        if (tauMothGenP.hasUserInt("bMothIndex") ){
	  tauH.addUserInt ("bMothIndex", tauMothGenP.userInt ("bMothIndex"));
	  tauH_neutral.addUserInt ("bMothIndex", tauMothGenP.userInt ("bMothIndex"));
	}
        if (tauMothGenP.hasUserInt("WMothIndex") ){
	  tauH.addUserInt ("WMothIndex", tauMothGenP.userInt ("WMothIndex"));
	  tauH_neutral.addUserInt ("WMothIndex", tauMothGenP.userInt ("WMothIndex"));
	}
        if (tauMothGenP.hasUserInt("ZMothIndex") ){
	  tauH.addUserInt ("ZMothIndex", tauMothGenP.userInt ("ZMothIndex"));
	  tauH_neutral.addUserInt ("ZMothIndex", tauMothGenP.userInt ("ZMothIndex"));
	}

	
        
        // many flags change of meaning w.r.t. mother tau, put everything to 0 (can be changed in future)
        int tauhFlags = 0;        
        tauH.addUserInt ("generalGenFlags", tauhFlags); // remember! TauH inherits ALL the flags from
	tauH_neutral.addUserInt ("generalGenFlags", tauhFlags); // remember! TauH inherits ALL the flags from 
        
        if (DEBUG){
	  cout << " ++ " << iTauH << " id: " << tauH.pdgId() << " | pt: " << tauH.pt() << " | eta: " << tauH.eta() << endl;
	  cout << " ++ " << iTauH << " id: " << tauH_neutral.pdgId() << " | pt: " << tauH_neutral.pt() << " | eta: " << tauH_neutral.eta() << endl;
	}
        result->push_back (tauH);
	result->push_back (tauH_neutral);
    }        
    //iEvent.put(result);
	iEvent.put(std::move(result));
}

// set of requirement(s) defining which particle must be saved
bool GenFiller::IsInteresting (const GenParticle& p)
{
    int APdgId = abs(p.pdgId());

    bool IsLast = genhelper::IsLastCopy(p);
    bool GoodPdgId = (APdgId == 25 || APdgId == 36 || APdgId == 23 || APdgId == 24 ||// bosons
                      APdgId == 1000022 || APdgId == 1000023 || APdgId == 1000025  ||// SUSY particles
                      APdgId == 6  || // quarks
                      APdgId == 11 || APdgId == 12 || APdgId == 13 || APdgId == 14 || APdgId == 15 || APdgId == 16); // leptons

    if(isVBFParton(p)) return true ;
            
    if (IsLast && GoodPdgId) return true;
    
    // case of b quarks, just save first one (too many showering products)
    bool IsFirst = genhelper::IsFirstCopy(p, true);
    bool GoodFirstPdg = (APdgId == 5 || APdgId == 6 || APdgId == 11 || APdgId == 13 || APdgId == 15 || APdgId==25);
    if (storeLightFlavAndGlu_) // also light flavors and quarks
        GoodFirstPdg = (GoodFirstPdg || APdgId == 1 || APdgId == 2 || APdgId == 3 || APdgId == 4 || APdgId == 21);
    
    if (GoodFirstPdg && IsFirst) return true;

    // if ((APdgId == 1 || APdgId == 2 || APdgId == 3 || APdgId == 4 || APdgId == 5 || APdgId == 6 || APdgId == 11 || APdgId == 13 || APdgId == 15 || APdgId == 21 || APdgId==25) && IsFirst) return true; // for b, save also the first copy in the list
    // if (APdgId == 5 && IsFirst) return true; // for b, save also the first copy in the list
                                             // check abs id to avoid problems if b -> (b bar) b as the bbar woudl result as first

    return false;
}


int GenFiller::makeFlagVector (const GenParticle* p)
{
    int flags = 0;
    const GenStatusFlags& fl = p->statusFlags();
    
    if (fl.isPrompt())                  flags |= (1 << 0);
    if (fl.isDecayedLeptonHadron())     flags |= (1 << 1);
    if (fl.isTauDecayProduct())         flags |= (1 << 2);
    if (fl.isPromptTauDecayProduct())   flags |= (1 << 3);
    if (fl.isDirectTauDecayProduct())   flags |= (1 << 4);
    if (fl.isDirectPromptTauDecayProduct())       flags |= (1 << 5);
    if (fl.isDirectHadronDecayProduct())          flags |= (1 << 6);
    if (fl.isHardProcess())             flags |= (1 << 7);
    if (fl.fromHardProcess())           flags |= (1 << 8);
    if (fl.isHardProcessTauDecayProduct())       flags |= (1 << 9);
    if (fl.isDirectHardProcessTauDecayProduct()) flags |= (1 << 10);
    if (fl.fromHardProcessBeforeFSR())  flags |= (1 << 11);
    if (fl.isFirstCopy())               flags |= (1 << 12);
    if (fl.isLastCopy())                flags |= (1 << 13);
    if (fl.isLastCopyBeforeFSR())       flags |= (1 << 14);
    if (isVBFParton(*p))                 flags |= (1 << 15);   
    return flags;
}

bool GenFiller::isVBFParton(const GenParticle& p)
{
  int APdgId = abs(p.pdgId());
  bool IsVBFPartonPdgId = (APdgId == 1 || APdgId == 2 || APdgId == 3 || APdgId == 4 || APdgId == 5 || APdgId == 21);// quark or gluon
  bool FoundHiggs = false;
  
  if(IsVBFPartonPdgId)
    {
      for(unsigned int iMother = 0 ; iMother < p.numberOfMothers() ; ++iMother)
	{
	  const reco::Candidate* Mother = p.mother(iMother);
	  for(unsigned int iDaughter = 0 ; iDaughter < Mother->numberOfDaughters() ; ++iDaughter)
	    {
	      const reco::Candidate* Daughter = Mother->daughter(iDaughter);
	      if(Daughter->pdgId()==25) FoundHiggs = true ;
	      if(FoundHiggs) break;
	    }
	}
    }
  return FoundHiggs ;
}


#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(GenFiller);
