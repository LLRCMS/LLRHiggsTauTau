/* \class ZrecoilCorrectionProducer
**
** This class applies the Z-recoil corrections 
** to correct simulated events for data/MC differences in MET response and resolution.
** 
** input type is collection of uncorrected MET objects
**
** output type is collection of corrected MET objects
**  
** \date:    25 March 2016
** \author:  C. Veelken (Tallinn)
*/

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Utilities/interface/InputTag.h>
#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/MET.h>
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include <TLorentzVector.h>

#include <vector>
#include <string>

using namespace edm;
using namespace std;
using namespace reco;

typedef edm::View<reco::Candidate> CandidateView;

// ------------------------------------------------------------------

class ZrecoilCorrectionProducer : public edm::EDProducer {
 public:
  /// Constructor
  explicit ZrecoilCorrectionProducer(const edm::ParameterSet&);
    
  /// Destructor
  ~ZrecoilCorrectionProducer();  

 private:
  virtual void beginJob(){};  
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};

  edm::InputTag srcPairs_;
  edm::InputTag srcMEt_;
  edm::InputTag srcGenParticles_;
  edm::InputTag srcJets_;
  std::string correction_;
  RecoilCorrector* recoilCorrector_;
};

// ------------------------------------------------------------------


ZrecoilCorrectionProducer::ZrecoilCorrectionProducer(const edm::ParameterSet& cfg)
  : recoilCorrector_(0)
{
  srcPairs_ = cfg.getParameter<edm::InputTag>("srcPairs");
  consumes<CompositeCandidateView>(srcPairs_);
  srcMEt_ = cfg.getParameter<edm::InputTag>("srcMEt");
  consumes<pat::METCollection>(srcMEt_);
  srcGenParticles_ = cfg.getParameter<edm::InputTag>("srcGenParticles");
  consumes<reco::GenParticleCollection>(srcGenParticles_);
  srcJets_ = cfg.getParameter<edm::InputTag>("srcJets");
  consumes<CandidateView>(srcJets_);
  correction_ = cfg.getParameter<std::string>("correction");
  recoilCorrector_ = new RecoilCorrector(correction_);

  produces<pat::METCollection>();  
}  

ZrecoilCorrectionProducer::~ZrecoilCorrectionProducer()
{
  delete recoilCorrector_;
}

void ZrecoilCorrectionProducer::produce(edm::Event& evt, const edm::EventSetup& es)
{  
  //std::auto_ptr<pat::METCollection> result(new pat::METCollection());
  std::unique_ptr<pat::METCollection> result(new pat::METCollection());

  // retrieve gen info
  edm::Handle<reco::GenParticleCollection> genParticles;
  evt.getByLabel(srcGenParticles_, genParticles);

  reco::Candidate::LorentzVector genBosonP4;
  reco::Candidate::LorentzVector genBosonVisP4;
  for ( reco::GenParticleCollection::const_iterator genParticle = genParticles->begin();
	genParticle != genParticles->end(); ++genParticle ) {
    bool fromHardProcessFinalState = genParticle->fromHardProcessFinalState();
    int absPdgId = abs(genParticle->pdgId());
    bool isElectron = (absPdgId == 11);
    bool isMuon = (absPdgId == 13);
    bool isNeutrino = (absPdgId == 12 || absPdgId == 14 || absPdgId == 16);
    bool isDirectHardProcessTauDecayProduct = genParticle->isDirectHardProcessTauDecayProductFinalState();
    if ( (fromHardProcessFinalState && (isMuon || isElectron || isNeutrino)) || isDirectHardProcessTauDecayProduct ) {
      genBosonP4 += genParticle->p4();
    }
    if ( (fromHardProcessFinalState && (isMuon || isElectron)) || (isDirectHardProcessTauDecayProduct && !isNeutrino) ) { 
      genBosonVisP4 += genParticle->p4();
    }
  }
  
  bool isW = false;
  bool isZ = false;
  bool isHiggs = false;
  for ( reco::GenParticleCollection::const_iterator genParticle = genParticles->begin();
	genParticle != genParticles->end(); ++genParticle ) {
    int absPdgId = abs(genParticle->pdgId());
    if      ( absPdgId == 24                                     ) isW = true;
    else if ( absPdgId == 22 || absPdgId == 23                   ) isZ = true;
    else if ( absPdgId == 25 || absPdgId == 35 || absPdgId == 36 ) isHiggs = true;
    bool isChargedLepton = (absPdgId == 11 || absPdgId == 13 || absPdgId == 15);
    if ( isW || isZ || isHiggs || isChargedLepton ) break;
  }
  if ( !(isW || isZ || isHiggs) ) isZ = true; // CV: off-shell Z/gamma* -> ll events may have no Z-boson in genParticle list

  // loop on al lepton pairs
  edm::Handle<CompositeCandidateView> pairs;
  evt.getByLabel(srcPairs_, pairs);

  edm::Handle<CandidateView> jets;
  evt.getByLabel(srcJets_, jets);

  edm::Handle<pat::METCollection> uncorrMEt;
  evt.getByLabel(srcMEt_, uncorrMEt);

  unsigned int nPairs = pairs->size();
  for (unsigned int iPair = 0; iPair < nPairs; ++iPair)
  {
    const CompositeCandidate& pair = (*pairs)[iPair];
    const Candidate *l1 = pair.daughter(0);
    const Candidate *l2 = pair.daughter(1);

    int nJets = 0;
    for ( CandidateView::const_iterator jet = jets->begin(); jet != jets->end(); ++jet )
    {
      bool isLepton1 = ( deltaR(*jet, *l1) < 0.5 ? true : false );
      bool isLepton2 = ( deltaR(*jet, *l2) < 0.5 ? true : false );

      if ( !isLepton1 ) ++nJets;
      if ( !isLepton2 ) ++nJets;
    }
    if ( isW ) ++nJets; // CV: add jet that fakes the hadronic tau candidate

    const pat::MET& theUncorrMEt = uncorrMEt->at(iPair);
  
    float corrMEtPx, corrMEtPy;
    recoilCorrector_->CorrectByMeanResolution(
      theUncorrMEt.px(),
      theUncorrMEt.py(),
      genBosonP4.px(),
      genBosonP4.py(),
      genBosonVisP4.px(),
      genBosonVisP4.py(),
      nJets,
      corrMEtPx,
      corrMEtPy
    );
    reco::Candidate::LorentzVector corrMEtP4(corrMEtPx, corrMEtPy, 0., sqrt(corrMEtPx*corrMEtPx + corrMEtPy*corrMEtPy));

    pat::MET corrMEt(theUncorrMEt);
    corrMEt.setP4(corrMEtP4);
    corrMEt.setSignificanceMatrix(theUncorrMEt.getSignificanceMatrix());
    corrMEt.addUserFloat("uncorrPx", theUncorrMEt.px());
    corrMEt.addUserFloat("uncorrPy", theUncorrMEt.py());
    result->push_back(corrMEt);
  }

  //evt.put(result);
  evt.put(std::move(result));
}

#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(ZrecoilCorrectionProducer);
