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

  edm::InputTag src_;
  std::vector<const reco::Candidate *> cands_;
  //std::vector<reco::GenParticle> tauHadcands_; // gen H tau build in this class
  std::vector<int> tauHadcandsMothers_; // contains the index in the cands_ vector of the tauh mother
};

// ------------------------------------------------------------------

GenFiller::GenFiller(const edm::ParameterSet& iConfig)
{
    src_ = iConfig.getParameter<InputTag>("src");
    produces<pat::GenericParticleCollection>();
}

void GenFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    Handle <edm::View<reco::GenParticle> > genHandle;
    iEvent.getByLabel (src_, genHandle);
    cands_.clear();
    //tauHadcands_.clear();
    tauHadcandsMothers_.clear();
    
    // output collection
    auto_ptr<pat::GenericParticleCollection> result( new pat::GenericParticleCollection );
    
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
        int PdgId = genP->pdgId();
        int APdgId = abs(PdgId);

        //only select some particles
        pat::GenericParticle filtGenP (*genP); // to be saved in output
        if (DEBUG) cout << iGen << " | id: " << genP->pdgId () << "   pt: " << genP->pt() << "   eta: " << genP->eta() << " | px: " << genP->px() << " , eta: " << genP->eta() << endl;
        
        // ------------------- tau decay flags
        if (APdgId == 15)
        {
            int decay = genhelper::GetTauDecay(genP);
            filtGenP.addUserInt ("tauGenDecayMode", decay);
            if (decay == 2) tauHadcandsMothers_.push_back (iGen); // for later usage for tauh
            if (DEBUG) cout << "   --> tau decay: " << decay << endl;
        }

        // -------------------- H/Z decay mode: set final state and is is prompt
        if (APdgId == 25 || APdgId == 23)
        {
            genhelper::HZDecay decay = genhelper::GetHZDecay (genP);
            filtGenP.addUserInt ("HZDecayMode", static_cast<int> (decay));
            if (DEBUG) cout << "   --> H/Z decay: " << decay << endl;
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

        // tau
        MothPtr = genhelper::IsFromID (genP, 15);
        if (MothPtr != NULL)
        {
            //filtGenP.addUserInt ("fromTau", 1);
            filtGenP.addUserInt ("TauMothIndex", genhelper::GetIndexInOutput(MothPtr, cands_));       
            if (DEBUG) cout << "   --> fromTau: 1, indexTau: " << genhelper::GetIndexInOutput(MothPtr, cands_) << endl;
        }

        result->push_back (filtGenP);
    }


    // finally, do the hadronic tau (leave as last step not to spoil internal ordering of mother / daughter vector)    
    if (DEBUG) cout << "BUILDING tauH, size: " << tauHadcandsMothers_.size() << endl;
    for (unsigned int iTauH = 0; iTauH < tauHadcandsMothers_.size(); iTauH++)
    {
        int tauMothInd = tauHadcandsMothers_.at(iTauH);
        pat::GenericParticle tauH (genhelper::GetTauHad(cands_.at(tauMothInd)));       
        tauH.addUserInt ("TauMothIndex", tauMothInd);
        if (DEBUG) cout << " ++ " << iTauH << " id: " << tauH.pdgId() << " | pt: " << tauH.pt() << " | eta: " << tauH.eta() << endl;
        result->push_back (tauH);
    }        
    iEvent.put(result);
}

// set of requirement(s) defining which particle must be saved
bool GenFiller::IsInteresting (const GenParticle& p)
{
    int APdgId = abs(p.pdgId());
    
    bool IsLast = genhelper::IsLastCopy(p);
    //bool IsFirst = genhelper::IsFirstCopy(p);
    bool GoodPdgId = (APdgId == 25 || APdgId == 23 || // bosons
                      APdgId == 5 || APdgId == 6 || // quarks
                      APdgId == 11 || APdgId == 12 || APdgId == 13 || APdgId == 14 || APdgId == 15 || APdgId == 16); // leptons
            
    if (IsLast && GoodPdgId) return true;
    
    /*
    if (APdgId == 5 && IsFirst) return true; // for b, save also the first copy in the list
                                             // --> note: might give problems if b -> (b bar) b as the bbar woudl result as first
    */
            
    return false;
}



#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(GenFiller);