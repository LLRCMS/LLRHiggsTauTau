/** \class MuFiller
 *
 *  No description available.
 *
 *  $Date: 2013/05/14 10:08:19 $
 *  $Revision: 1.18 $
 *  \author N. Amapane (Torino)
 *  \author S. Bolognesi (JHU)
 *  \author C. Botta (CERN)
 *  \author S. Casasso (Torino)
 */

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/ESHandle.h>

#include <DataFormats/PatCandidates/interface/Muon.h>
#include <Muon/MuonAnalysisTools/interface/MuonEffectiveArea.h>
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalTrajectoryExtrapolatorToLine.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"

#include <LLRHiggsTauTau/NtupleProducer/interface/DaughterDataHelpers.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/CutSet.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/LeptonIsoHelper.h>
//#include <ZZAnalysis/AnalysisStep/interface/SIPCalculator.h>


#include <vector>
#include <string>

 using namespace edm;
 using namespace std;
 using namespace reco;

 class MuFiller : public edm::EDProducer {
 public:
  /// Constructor
  explicit MuFiller(const edm::ParameterSet&);

  /// Destructor
  ~MuFiller(){};  

private:
  virtual void beginJob(){};  
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};

  //const edm::InputTag theCandidateTag;
  //edm::EDGetTokenT<edm::View<pat::Muon> > theCandidateTag;
  edm::EDGetTokenT<pat::MuonRefVector> theCandidateTag;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > theGenTag ;
  edm::EDGetTokenT<double> theRhoTag ;
  edm::EDGetTokenT<vector<Vertex>> theVtxTag ;
  int sampleType;
  int setup;
  const StringCutObjectSelector<pat::Muon, true> cut;
  const CutSet<pat::Muon> flags;
  //SIPCalculator *sipCalculator_;

};


MuFiller::MuFiller(const edm::ParameterSet& iConfig) :
//theCandidateTag(iConfig.getParameter<InputTag>("src")),
//theCandidateTag(consumes<edm::View<pat::Muon> >(iConfig.getParameter<InputTag>("src"))),
theCandidateTag(consumes<pat::MuonRefVector>(iConfig.getParameter<InputTag>("src"))),
theGenTag(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genCollection"))),
theRhoTag(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoCollection"))),
theVtxTag(consumes<vector<Vertex>>(iConfig.getParameter<edm::InputTag>("vtxCollection"))),
sampleType(iConfig.getParameter<int>("sampleType")),
setup(iConfig.getParameter<int>("setup")),
cut(iConfig.getParameter<std::string>("cut")),
flags(iConfig.getParameter<edm::ParameterSet>("flags"))
{
  //sipCalculator_ = new SIPCalculator();
  produces<pat::MuonCollection>();
}


void
MuFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //Initialize SIP calculator
  //sipCalculator_->initialize(iSetup);

  //--- Get leptons and rho
  edm::Handle<pat::MuonRefVector> muonHandle;
  //Handle<edm::View<reco::GenParticle> > genHandle;
  //iEvent.getByLabel(theCandidateTag, muonHandle);
  iEvent.getByToken(theCandidateTag, muonHandle);

  //InputTag theRhoTag = LeptonIsoHelper::getMuRhoTag(sampleType, setup);
  edm::Handle<double> rhoHandle;
  //iEvent.getByLabel(theRhoTag, rhoHandle);//RH;
  iEvent.getByToken(theRhoTag, rhoHandle);
  double rho = *rhoHandle;//RH
  //double rho = 1.;

  edm::Handle<vector<Vertex> >  vertexs;
  iEvent.getByToken(theVtxTag,vertexs);
  //iEvent.getByLabel("offlineSlimmedPrimaryVertices",vertexs);

  // Output collection
  //auto_ptr<pat::MuonCollection> result( new pat::MuonCollection() );
  std::unique_ptr<pat::MuonCollection> result( new pat::MuonCollection() );

  for (unsigned int i = 0; i< muonHandle->size(); ++i){
    //---Clone the pat::Muon
    pat::Muon l(*((*muonHandle)[i].get()));

    //--- PF ISO
    float PFChargedHadIso   = l.chargedHadronIso();
    float PFNeutralHadIso   = l.neutralHadronIso();
    float PFPhotonIso       = l.photonIso();
    
    float combRelIsoPF = LeptonIsoHelper::combRelIsoPF(sampleType, setup, rho, l);
    float combRelIsoPF03 = LeptonIsoHelper::combRelIsoPF(sampleType, setup, rho, l, true);
    //--- SIP, dxy, dz
    float IP      = std::abs(l.dB(pat::Muon::PV3D));
    float IPError = l.edB(pat::Muon::PV3D);
    float SIP     = IP/IPError;

    float dxy = 999.;
    float dz  = 999.;
    float dxy_innerTrack = 999.;
    float dz_innerTrack  = 999.;
    const Vertex* vertex = 0;

    if (vertexs->size()>0) {
      vertex = &(vertexs->front());
      dxy = (l.muonBestTrack()->dxy(vertex->position()));
      dz  = (l.muonBestTrack()->dz(vertex->position()));
      dxy_innerTrack = (l.innerTrack()->dxy(vertex->position()));
      dz_innerTrack  = (l.innerTrack()->dz(vertex->position()));
    }

    float rel_error_trackpt = l.muonBestTrack()->ptError()/l.muonBestTrack()->pt();    

    /*
    //--- Trigger matching
    bool HLTMatch = ((!l.triggerObjectMatchesByFilter("hltSingleMu13L3Filtered17").empty())||
		     ((!l.triggerObjectMatchesByFilter("hltSingleMu13L3Filtered13").empty()) && 
		      (l.triggerObjectMatchesByFilter("hltSingleMu13L3Filtered13").at(0).pt()>17)) || 
		     ((!l.triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered5").empty()) && 
		      (l.triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered5").at(0).pt()>17)) || 
		     ((!l.triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered7").empty()) && 
		      (l.triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered7").at(0).pt()>17)));
    //FIXME*/

    //--- Embed user variables
    l.addUserFloat("PFChargedHadIso",PFChargedHadIso);
    l.addUserFloat("PFNeutralHadIso",PFNeutralHadIso);
    l.addUserFloat("PFPhotonIso",PFPhotonIso);
    l.addUserFloat("combRelIsoPF",combRelIsoPF);
    l.addUserFloat("combRelIsoPF03",combRelIsoPF03);
    l.addUserFloat("rho",rho);
    l.addUserFloat("DepositR03TrackerOfficial",l.isolationR03().sumPt);
    l.addUserFloat("DepositR03ECal",l.isolationR03().emEt);
    l.addUserFloat("DepositR03Hcal",l.isolationR03().hadEt);
    //    if(vertexs->size()>0) {
    //      SIP=sipCalculator_->calculate(l,vertexs->front());
    //    }
    l.addUserFloat("SIP",SIP);
    l.addUserFloat("segmentCompatibility",l.segmentCompatibility());
    
    l.addUserFloat("dxy",dxy);
    l.addUserFloat("dz",dz);
    l.addUserFloat("dxy_innerTrack",dxy_innerTrack);
    l.addUserFloat("dz_innerTrack",dz_innerTrack);
    l.addUserFloat("rel_error_trackpt",rel_error_trackpt);


    //l.addUserFloat("HLTMatch", HLTMatch);
    // l.addUserCand("MCMatch",genMatch); // FIXME
    int idbit=0;
    if(l.isLooseMuon())idbit |= 1 << 0;
    bool global =  l.isGlobalMuon();
    if(l.isMediumMuon())idbit |= 1 << 2;
    int normX2= 999;
    if(vertex){
      if(l.globalTrack().isNonnull()) normX2 = l.globalTrack()->normalizedChi2();
      if(l.isSoftMuon(vertexs->front())) idbit |= 1 << 1;
      if(l.isTightMuon(vertexs->front())) idbit |= 1 << 3;
      if(l.isHighPtMuon(vertexs->front())) idbit |= 1 << 4;
    }
    if(l.isLooseMuon() && global && normX2<10 && l.globalTrack().isNonnull() && l.innerTrack().isNonnull()){
      if(l.globalTrack()->hitPattern().numberOfValidMuonHits() > 0 && l.numberOfMatchedStations()>1){
	if(l.innerTrack()->hitPattern().numberOfValidPixelHits()>0 && l.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5)idbit |= 1 << 5;
      }
    }

    bool goodGlob = l.isGlobalMuon() && 
      l.globalTrack()->normalizedChi2() < 3 && 
      l.combinedQuality().chi2LocalPosition < 12 && 
      l.combinedQuality().trkKink < 20; 
    bool isMedium = muon::isLooseMuon(l) && 
      l.innerTrack()->validFraction() > 0.49 && 
      muon::segmentCompatibility(l) > (goodGlob ? 0.303 : 0.451); 
    
    if(isMedium) idbit |= 1 << 6;

    l.addUserInt("muonID",idbit);
    //--- isPFMuon flag - in old samples, l.isPFMuon() is not functional, so this has to be filled
    //    beforehand with the module PATPFMuonEmbedder.
    if(!l.hasUserFloat("isPFMuon")) {
      l.addUserFloat("isPFMuon",l.isPFMuon());
    }
    if(!l.hasUserFloat("isGlobalMuon")) {
      l.addUserFloat("isGlobalMuon",l.isGlobalMuon());
    }
    if(!l.hasUserFloat("isTrackerMuon")) {
      l.addUserFloat("isTrackerMuon",l.isTrackerMuon());
    }
    
    //--- MC parent code 
    const reco::GenParticle* genL= l.genParticleRef().get();
    float px=0,py=0,pz=0,E=0,fromH=0;
    int status=99999, id=99999;
    
    //printf("%p\n",genL);
    if(genL){
      px =genL->p4().Px();
      py =genL->p4().Py();
      pz =genL->p4().Pz();
      E =genL->p4().E();
      status =genL->status();
      id =genL->pdgId();
      
      //cout << "Mu filler: " << i << " [px, id] = " << l.px() << " , " << l.pdgId() << " | (px, py, pz, e) " << px << " " << py << " " << pz << " " << E << " | ID: " << genL->pdgId() << " | status: " << genL->status() << endl;
      
      //search if it comes from H
      Handle<edm::View<reco::GenParticle> > genHandle;
      iEvent.getByToken(theGenTag, genHandle);
      int nmot = genL->numberOfMothers();
      for (int im = 0; im<nmot&&fromH==0; ++im){
	for(unsigned int ipruned = 0; ipruned< genHandle->size(); ++ipruned){
	  if(ipruned==i)continue;
	  if(userdatahelpers::isAncestor(&(*genHandle)[ipruned],genL)){
	    int pdgmot = (&(*genHandle)[ipruned])->pdgId();
	    if(abs(pdgmot)==25){
	      fromH=1;
	      break;
	    }
	  }
	}
      }
    }
    l.addUserFloat("genPx",px);
    l.addUserFloat("genPy",py);
    l.addUserFloat("genPz",pz);
    l.addUserFloat("genE",E);    
    l.addUserInt("status", status);
    l.addUserInt("id", id);
    l.addUserFloat("fromH",fromH);
    
    //     MCHistoryTools mch(iEvent);
    //     if (mch.isMC()) {
    //       int MCParentCode = 0;//FIXME: does not work on cmg mch.getParentCode((l.genParticleRef()).get());
    //       l.addUserFloat("MCParentCode",MCParentCode);
    //     }
    
    //--- Check selection cut. Being done here, flags are not available; but this way we 
    //    avoid wasting time on rejected leptons.
    if (!cut(l)) continue;

    //--- Embed flags (ie flags specified in the "flags" pset)
    for(CutSet<pat::Muon>::const_iterator flag = flags.begin(); flag != flags.end(); ++flag) {
      l.addUserFloat(flag->first,int((*(flag->second))(l)));
    }
    
    result->push_back(l);
  }
  //iEvent.put(result);
  iEvent.put(std::move(result));
}

#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(MuFiller);

