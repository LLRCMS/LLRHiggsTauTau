#ifndef RecoMET_METPUSubtraction_PFMETProducerMVAInputPairs_h
#define RecoMET_METPUSubtraction_PFMETProducerMVAInputPairs_h

/** \class PFMETProducerMVAInputPairs
 *
 * Produce PFMET objects computed by MVA 
 *
 * \authors Phil Harris, CERN
 *          Christian Veelken, LLR
 * \modified by L. Cadamuro, LLR, to handle reco::CompositeCandidate as input
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "RecoMET/METAlgorithms/interface/METAlgo.h"
#include "RecoMET/METAlgorithms/interface/PFSpecificAlgo.h"
#include "RecoMET/METPUSubtraction/interface/PFMETAlgorithmMVA.h"
#include "RecoMET/METPUSubtraction/interface/mvaMEtUtilities.h"

#include "RecoJets/JetProducers/interface/PileupJetIdAlgo.h"

#include <vector>

namespace reco
{
  class PFMETProducerMVAInputPairs : public edm::stream::EDProducer<>
  {
   public:

    PFMETProducerMVAInputPairs(const edm::ParameterSet&); 
    ~PFMETProducerMVAInputPairs();

   private:
  
    void produce(edm::Event&, const edm::EventSetup&);

    // auxiliary functions
    std::vector<mvaMEtUtilities::JetInfo> computeJetInfo(const reco::PFJetCollection&, const reco::PFJetCollection&, const reco::VertexCollection&, const reco::Vertex*, 
							 const JetCorrector &iCorr,edm::Event & iEvent,const edm::EventSetup &iSetup,
							 std::vector<mvaMEtUtilities::leptonInfo> &iLeptons,std::vector<mvaMEtUtilities::pfCandInfo> &iCands);
    
    std::vector<mvaMEtUtilities::pfCandInfo> computePFCandidateInfo(const reco::PFCandidateCollection&, const reco::Vertex*);
    std::vector<reco::Vertex::Point> computeVertexInfo(const reco::VertexCollection&);
    double chargedFrac(const reco::Candidate *iCand,const reco::PFCandidateCollection& pfCandidates,const reco::Vertex* hardScatterVertex);
    
    bool   passPFLooseId(const reco::PFJet *iJet);
    bool   istau        (const reco::Candidate *iCand);
    double chargedFracInCone(const reco::Candidate *iCand,const reco::PFCandidateCollection& pfCandidates,const reco::Vertex* hardScatterVertex,double iDRMax=0.2);

   // configuration parameter
    edm::EDGetTokenT<reco::PFJetCollection> srcCorrJets_;
    edm::EDGetTokenT<reco::PFJetCollection> srcUncorrJets_;
    edm::EDGetTokenT<reco::PFCandidateCollection> srcPFCandidates_;
    edm::EDGetTokenT<edm::View<reco::Candidate> > srcPFCandidatesView_;
    edm::EDGetTokenT<reco::VertexCollection> srcVertices_;
    //typedef std::vector<edm::InputTag> vInputTag;
    edm::EDGetTokenT<pat::CompositeCandidate> srcLeptons_;
    int minNumLeptons_; // CV: option to skip MVA MET computation in case there are less than specified number of leptons in the event
    edm::EDGetTokenT<edm::Handle<double> > srcRho_;

    std::string correctorLabel_;
    bool isOld42_ ;
    bool useType1_;
    
    double globalThreshold_;

    double minCorrJetPt_;

    METAlgo metAlgo_;
    PFSpecificAlgo pfMEtSpecificAlgo_;
    PFMETAlgorithmMVA mvaMEtAlgo_;
    bool mvaMEtAlgo_isInitialized_;
    PileupJetIdAlgo mvaJetIdAlgo_;

    int verbosity_;
  };
}

#endif
