/** \class EleFiller
 *
 *  No description available.
 *
 *  $Date: 2013/05/24 15:42:42 $
 *  $Revision: 1.28 $
 *  \author N. Amapane (Torino)
 *  \author S. Bolognesi (JHU)
 *  \author C. Botta (CERN)
 *  \author S. Casasso (Torino)
 *  \author G. Ortona (LLR)
 *  \author G. Topo (LLR)
 */
 
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
 #include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
#include "Math/VectorUtil.h"
#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include "DataFormats/Common/interface/ValueMap.h"
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/DaughterDataHelpers.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/CutSet.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/LeptonIsoHelper.h>
//#include "EgammaAnalysis/ElectronTools/interface/EGammaMvaEleEstimatorCSA14.h"

#include <vector>
#include <string>

using namespace edm;
using namespace std;
using namespace reco;


//bool recomputeBDT = false;

class EleFiller : public edm::EDProducer {
 public:
  /// Constructor
  explicit EleFiller(const edm::ParameterSet&);
    
  /// Destructor
  ~EleFiller(){
    //delete bdt;
  };  

 private:
  virtual void beginJob(){};  
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};

  //const edm::InputTag theCandidateTag;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > theGenTag ;
  edm::EDGetTokenT<double> theRhoTag ;
  edm::EDGetTokenT<vector<Vertex>> theVtxTag ;
  int sampleType;
  int setup;
  const StringCutObjectSelector<pat::Electron, true> cut;
  //const CutSet<pat::Electron> flags;
  //EGammaMvaEleEstimatorCSA14* myMVATrig;
  edm::EDGetTokenT<edm::View<pat::Electron> > electronCollectionToken_;

  //BDTId* bdt;
};


EleFiller::EleFiller(const edm::ParameterSet& iConfig) :
  //theCandidateTag(iConfig.getParameter<InputTag>("src")),
  theGenTag(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genCollection"))),
  theRhoTag(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoCollection"))),
  theVtxTag(consumes<vector<Vertex>>(iConfig.getParameter<edm::InputTag>("vtxCollection"))),
  sampleType(iConfig.getParameter<int>("sampleType")),
  setup(iConfig.getParameter<int>("setup")),
  cut(iConfig.getParameter<std::string>("cut")),
  //flags(iConfig.getParameter<ParameterSet>("flags")),//,
  //myMVATrig(0),
  electronCollectionToken_(consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("src")))
 
  //bdt(0)
{
  produces<pat::ElectronCollection>();
}

  void EleFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
  {  
  // Get leptons and rho
 
  edm::Handle<edm::View<pat::Electron> > electrons;
 
  iEvent.getByToken(electronCollectionToken_, electrons);

 
  edm::Handle<double> rhoHandle;
  iEvent.getByToken(theRhoTag, rhoHandle); 
  double rho = *rhoHandle;

  edm::Handle<vector<Vertex> >  vertexs;

  iEvent.getByToken(theVtxTag,vertexs);
    
//**********************


  // Output collection
  std::unique_ptr<pat::ElectronCollection> result( new pat::ElectronCollection() );

  //unsigned int i=0;
  //for (unsigned int i = 0; i< electronHandle->size(); ++i){
  for( View<pat::Electron>::const_iterator el = electrons->begin(); el != electrons->end(); el++){
    const unsigned int i = distance (electrons->begin(), el);
    const auto ele = electrons->ptrAt(i);

    //---Clone the pat::Electron
 
    pat::Electron l(*el);
    
    //--- PF ISO
    float PFChargedHadIso   = l.chargedHadronIso();
    float PFNeutralHadIso   = l.neutralHadronIso();
    float PFPhotonIso       = l.photonIso();

    //float fSCeta = fabs(l.eta()); 
    float fSCeta = l.superCluster()->eta();

    float combRelIsoPF = (l.pfIsolationVariables().sumChargedHadronPt + max(l.pfIsolationVariables().sumNeutralHadronEt +l.pfIsolationVariables().sumPhotonEt - 0.5 * l.pfIsolationVariables().sumPUPt, 0.0)) / l.pt();
//LeptonIsoHelper::combRelIsoPF(sampleType, setup, rho, l);

    //--- SIP, dxy, dz
    float IP      = fabs(l.dB(pat::Electron::PV3D));
    float IPError = l.edB(pat::Electron::PV3D);
    float SIP     = IP/IPError;

    float dxy = 999.;
    float dz  = 999.;
    const Vertex* vertex = 0;
    if (vertexs->size()>0) {
      vertex = &(vertexs->front());
      dxy = (l.gsfTrack()->dxy(vertex->position()));
      dz  = (l.gsfTrack()->dz(vertex->position()));
    } 

    float rel_error_trackpt = l.gsfTrack()->ptError()/l.gsfTrack()->pt();   
    
    float isEleID80         = l.electronID("mvaEleID-Fall17-iso-V2-wp80");
    float isEleNoIsoID80    = l.electronID("mvaEleID-Fall17-noIso-V2-wp80");
    float isEleID90         = l.electronID("mvaEleID-Fall17-iso-V2-wp90");
    float isEleNoIsoID90    = l.electronID("mvaEleID-Fall17-noIso-V2-wp90");
    float isEleIDLoose      = l.electronID("mvaEleID-Fall17-iso-V2-wpLoose");
    float isEleNoIsoIDLoose = l.electronID("mvaEleID-Fall17-noIso-V2-wpLoose");

    float mvaValue_Iso = l.userFloat("ElectronMVAEstimatorRun2Fall17IsoV2Values");
    float mvaValue_HZZ = l.userFloat("ElectronMVAEstimatorRun2Spring16HZZV1Values");

    bool isconversionveto=l.passConversionVeto();
    
    //-- Missing hit  
    int missinghit = l.gsfTrack()->hitPattern().numberOfAllHits(HitPattern::MISSING_INNER_HITS);
    int missinglosthit = l.gsfTrack()->hitPattern().numberOfLostHits(HitPattern::MISSING_INNER_HITS);
	  
    //--- 3 charge assignement
    bool isGsfCtfScPixChargeConsistent=l.isGsfCtfScPixChargeConsistent();
    bool isGsfScPixChargeConsistent=l.isGsfScPixChargeConsistent();

    float IoEmIoP_ttH = 9e9;
    if(l.ecalEnergy()>0)
      IoEmIoP_ttH = (1.0/l.ecalEnergy() - l.eSuperClusterOverP()/l.ecalEnergy());
      
    //--- Scale and smearing corrections and uncertainties - https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2#Energy_Scale_and_Smearing

    float ecalTrkEnergyPostCorr    = l.userFloat("ecalTrkEnergyPostCorr");
    float ecalTrkEnergyErrPostCorr = l.userFloat("ecalTrkEnergyErrPostCorr");
    float energyScaleUp            = l.userFloat("energyScaleUp");
    float energyScaleDown          = l.userFloat("energyScaleDown");
    float energySigmaUp            = l.userFloat("energySigmaUp");
    float energySigmaDown          = l.userFloat("energySigmaDown"); 

    //--- Embed user variables
    l.addUserFloat("PFChargedHadIso",PFChargedHadIso);
    l.addUserFloat("PFNeutralHadIso",PFNeutralHadIso);
    l.addUserFloat("PFPhotonIso",PFPhotonIso);
    l.addUserFloat("combRelIsoPF",combRelIsoPF);
    l.addUserFloat("rho",rho);
    l.addUserFloat("SIP",SIP);
    l.addUserFloat("dxy",dxy);
    l.addUserFloat("dz",dz);
    l.addUserFloat("isEleID80",isEleID80);
    l.addUserFloat("isEleNoIsoID80",isEleNoIsoID80);
    l.addUserFloat("isEleID90",isEleID90);
    l.addUserFloat("isEleNoIsoID90",isEleNoIsoID90);
    l.addUserFloat("isEleIDLoose",isEleIDLoose);
    l.addUserFloat("isEleNoIsoIDLoose",isEleNoIsoIDLoose);
    l.addUserFloat("mvaValue_Iso",mvaValue_Iso);
    l.addUserFloat("mvaValue_HZZ",mvaValue_HZZ);
    l.addUserFloat("rel_error_trackpt",rel_error_trackpt);
    l.addUserInt("isConversionVeto",(isconversionveto ? 1 : 0));
    //l.addUserFloat("HLTMatch", HLTMatch);
    l.addUserInt("missingHit", missinghit);
    l.addUserInt("missingLostHit", missinglosthit);
    l.addUserInt("isGsfCtfScPixChargeConsistent",(isGsfCtfScPixChargeConsistent ? 1 : 0));
    l.addUserInt("isGsfScPixChargeConsistent",(isGsfScPixChargeConsistent ? 1 : 0));
    l.addUserFloat("sigmaIetaIeta",l.sigmaIetaIeta());
    l.addUserFloat("full5x5_sigmaIetaIeta",l.full5x5_sigmaIetaIeta());
    l.addUserFloat("hOverE",l.hcalOverEcal());
    l.addUserFloat("deltaEtaSuperClusterTrackAtVtx",l.deltaEtaSuperClusterTrackAtVtx());
    l.addUserFloat("deltaPhiSuperClusterTrackAtVtx",l.deltaPhiSuperClusterTrackAtVtx());
    l.addUserFloat("IoEmIoP",(1.0/l.ecalEnergy())-(1.0/l.p()));
    l.addUserFloat("IoEmIoP_ttH",IoEmIoP_ttH);
    l.addUserFloat("SCeta", fSCeta);
    l.addUserInt("isEB", int(l.isEB()));
    l.addUserFloat("ecalTrkEnergyPostCorr",ecalTrkEnergyPostCorr);
    l.addUserFloat("ecalTrkEnergyErrPostCorr",ecalTrkEnergyErrPostCorr);
    l.addUserFloat("energyScaleUp",energyScaleUp);
    l.addUserFloat("energyScaleDown",energyScaleDown);
    l.addUserFloat("energySigmaUp",energySigmaUp);
    l.addUserFloat("energySigmaDown",energySigmaDown); 
   
    //--- MC info
    const reco::GenParticle* genL= l.genParticleRef().get();
    float px=0,py=0,pz=0,E=0,fromH=0;
    int status=99999, id=99999;
    if(genL){
      px =genL->p4().Px();
      py =genL->p4().Py();
      pz =genL->p4().Pz();
      E =genL->p4().E();
      status =genL->status();
      id =genL->pdgId();
      
      //cout << "Ele filler: " << el - electrons->begin() << " [px, id] = " << l.px() << " , " << l.pdgId() << " | (px, py, pz, e) " << px << " " << py << " " << pz << " " << E << " | ID: " << genL->pdgId() << " | status: " << genL->status() << endl;

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


    //--- Check selection cut. Being done here, flags are not available; but this way we 
    //    avoid wasting time on rejected leptons.
    if (!cut(l)) continue;

    //--- Embed flags (ie flags specified in the "flags" pset)
    //for(CutSet<pat::Electron>::const_iterator flag = flags.begin(); flag != flags.end(); ++flag) {
    //  l.addUserFloat(flag->first,int((*(flag->second))( l)));
    // }

    result->push_back(l);

  }

  iEvent.put(std::move(result));
}


#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(EleFiller);

