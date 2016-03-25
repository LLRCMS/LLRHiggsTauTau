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
#include <DataFormats/PatCandidates/interface/MET.h>
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include <TLorentzVector.h>

#include <vector>
#include <string>

using namespace edm;
using namespace std;
using namespace reco;

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

  edm::InputTag srcLeptons_;
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
  srcLeptons_ = cfg.getParameter<edm::InputTag>("srcLeptons");
  srcMEt_ = cfg.getParameter<edm::InputTag>("srcMEt");
  srcGenParticles_ = cfg.getParameter<edm::InputTag>("srcGenParticles");
  srcJets_ = cfg.getParameter<edm::InputTag>("srcJets");
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

  typedef edm::View<reco::Candidate> CandidateView;
  edm::Handle<CandidateView> leptons;
  evt.getByLabel(srcLeptons_, leptons);

  edm::Handle<CandidateView> jets;
  evt.getByLabel(srcJets_, jets);
  int nJets = 0;
  for ( CandidateView::const_iterator jet = jets->begin();
	jet != jets->end(); ++jet ) {
    bool isLepton = false;
    for ( CandidateView::const_iterator lepton = leptons->begin();
	  lepton != leptons->end(); ++lepton ) {
      double dR = deltaR(jet->p4(), lepton->p4());
      if ( dR < 0.5 ) isLepton = true;
    }
    if ( !isLepton ) ++nJets;
  }
  if ( isW ) ++nJets; // CV: add jet that fakes the hadronic tau candidate

  edm::Handle<pat::METCollection> uncorrMEt;
  evt.getByLabel(srcMEt_, uncorrMEt);
  const pat::MET& theUncorrMEt = uncorrMEt->at(0);
  
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

  std::auto_ptr<pat::METCollection> result(new pat::METCollection());
  result->push_back(corrMEt);

  evt.put(result);
}

#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(ZrecoilCorrectionProducer);
