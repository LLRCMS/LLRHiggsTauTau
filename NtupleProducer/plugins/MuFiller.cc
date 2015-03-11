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

  const edm::InputTag theCandidateTag;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > theGenTag ;
  int sampleType;
  int setup;
  const StringCutObjectSelector<pat::Muon, true> cut;
  const CutSet<pat::Muon> flags;
  //SIPCalculator *sipCalculator_;

};


MuFiller::MuFiller(const edm::ParameterSet& iConfig) :
  theCandidateTag(iConfig.getParameter<InputTag>("src")),
  theGenTag(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genCollection"))),
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
  iEvent.getByLabel(theCandidateTag, muonHandle);

  InputTag theRhoTag = LeptonIsoHelper::getMuRhoTag(sampleType, setup);
  edm::Handle<double> rhoHandle;
  iEvent.getByLabel(theRhoTag, rhoHandle);//RH;
  double rho = *rhoHandle;//RH
  //double rho = 1.;

  edm::Handle<vector<Vertex> >  vertexs;
  iEvent.getByLabel("goodPrimaryVertices",vertexs);

  // Output collection
  auto_ptr<pat::MuonCollection> result( new pat::MuonCollection() );

  for (unsigned int i = 0; i< muonHandle->size(); ++i){
    //---Clone the pat::Muon
    pat::Muon l(*((*muonHandle)[i].get()));

    //--- PF ISO
    float PFChargedHadIso   = l.chargedHadronIso();
    float PFNeutralHadIso   = l.neutralHadronIso();
    float PFPhotonIso       = l.photonIso();
    
    float combRelIsoPF = LeptonIsoHelper::combRelIsoPF(sampleType, setup, rho, l);

    //--- SIP, dxy, dz
    float IP      = std::abs(l.dB(pat::Muon::PV3D));
    float IPError = l.edB(pat::Muon::PV3D);
    float SIP     = IP/IPError;

    float dxy = 999.;
    float dz  = 999.;
    const Vertex* vertex = 0;
    if (vertexs->size()>0) {
      vertex = &(vertexs->front());
      dxy = (l.innerTrack()->dxy(vertex->position()));
      dz  = (l.innerTrack()->dz(vertex->position()));
    }
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
    l.addUserFloat("rho",rho);
    l.addUserFloat("DepositR03TrackerOfficial",l.isolationR03().sumPt);
    l.addUserFloat("DepositR03ECal",l.isolationR03().emEt);
    l.addUserFloat("DepositR03Hcal",l.isolationR03().hadEt);
    //    if(vertexs->size()>0) {
    //      SIP=sipCalculator_->calculate(l,vertexs->front());
    //    }
    l.addUserFloat("SIP",SIP);

    l.addUserFloat("dxy",dxy);
    l.addUserFloat("dz",dz);
    //l.addUserFloat("HLTMatch", HLTMatch);
    // l.addUserCand("MCMatch",genMatch); // FIXME
    int idbit=0;
    if(l.isLooseMuon())idbit |= 1 << 0;
    if(vertex){
        if(l.isSoftMuon(vertexs->front())) idbit |= 1 << 1;
        if(l.isTightMuon(vertexs->front())) idbit |= 1 << 2;
    }
    l.addUserFloat("muonID",idbit);
    //--- isPFMuon flag - in old samples, l.isPFMuon() is not functional, so this has to be filled
    //    beforehand with the module PATPFMuonEmbedder.
    if(!l.hasUserFloat("isPFMuon")) {
      l.addUserFloat("isPFMuon",l.isPFMuon());
    }
    
    //--- MC parent code 
    const reco::GenParticle* genL= l.genParticleRef().get();
    float px=0,py=0,pz=0,E=0,fromH=0;
    //printf("%p\n",genL);
    if(genL){
      px =genL->p4().Px();
      py =genL->p4().Py();
      pz =genL->p4().Pz();
      E =genL->p4().E();

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
  iEvent.put(result);
}


#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(MuFiller);

