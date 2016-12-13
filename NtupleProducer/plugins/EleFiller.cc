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
  const CutSet<pat::Electron> flags;
  //EGammaMvaEleEstimatorCSA14* myMVATrig;
  edm::EDGetTokenT<edm::View<pat::Electron> > electronCollectionToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronVetoIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronLooseIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronMediumIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronTightIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > mvaValuesMapToken_;
  edm::EDGetTokenT<edm::ValueMap<int> > mvaCategoriesMapToken_;//*******



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
  flags(iConfig.getParameter<ParameterSet>("flags")),//,
  //myMVATrig(0),
  electronCollectionToken_(consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("src"))),
  electronVetoIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronVetoIdMap"))),
  electronLooseIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronLooseIdMap"))),
  electronMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronMediumIdMap"))),
  electronTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronTightIdMap"))),
  eleMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
  eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap"))),
  mvaValuesMapToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaValuesMap"))),
  mvaCategoriesMapToken_(consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("mvaCategoriesMap"))) //****************

 
  //bdt(0)
{/*
  //if (recomputeBDT) bdt = new BDTId;
  //Recompute BDT
  std::vector<std::string> myManualCatWeigths;
  myManualCatWeigths.push_back("EgammaAnalysis/ElectronTools/data/PHYS14FIX/EIDmva_EB1_5_oldscenario2phys14FIX_BDT.weights.xml");
  myManualCatWeigths.push_back("EgammaAnalysis/ElectronTools/data/PHYS14FIX/EIDmva_EB2_5_oldscenario2phys14FIX_BDT.weights.xml");
  myManualCatWeigths.push_back("EgammaAnalysis/ElectronTools/data/PHYS14FIX/EIDmva_EE_5_oldscenario2phys14FIX_BDT.weights.xml");
  myManualCatWeigths.push_back("EgammaAnalysis/ElectronTools/data/PHYS14FIX/EIDmva_EB1_10_oldscenario2phys14FIX_BDT.weights.xml");
  myManualCatWeigths.push_back("EgammaAnalysis/ElectronTools/data/PHYS14FIX/EIDmva_EB2_10_oldscenario2phys14FIX_BDT.weights.xml");
  myManualCatWeigths.push_back("EgammaAnalysis/ElectronTools/data/PHYS14FIX/EIDmva_EE_10_oldscenario2phys14FIX_BDT.weights.xml");

  vector<string> myManualCatWeigthsTrig;
  string the_path;
  for (unsigned i = 0 ; i < myManualCatWeigths.size() ; i++){
    the_path = edm::FileInPath ( myManualCatWeigths[i] ).fullPath();
    myManualCatWeigthsTrig.push_back(the_path);
  }
  myMVATrig = new EGammaMvaEleEstimatorCSA14();
  myMVATrig->initialize("BDT",
	  EGammaMvaEleEstimatorCSA14::kNonTrigPhys14,
	  true,
	  myManualCatWeigthsTrig);
  */
  produces<pat::ElectronCollection>();
}

  void EleFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
  {  
  // Get leptons and rho
  //edm::Handle<pat::ElectronRefVector> electronHandle;
  //iEvent.getByLabel(theCandidateTag, electronHandle);
  edm::Handle<edm::View<pat::Electron> > electrons;
  //edm::Handle<edm::View<reco::GsfElectron> > electrons;
  iEvent.getByToken(electronCollectionToken_, electrons);

  //InputTag theRhoTag = LeptonIsoHelper::getEleRhoTag(sampleType,setup);
  edm::Handle<double> rhoHandle;
  iEvent.getByToken(theRhoTag, rhoHandle); 
  double rho = *rhoHandle;

  edm::Handle<vector<Vertex> >  vertexs;
  //iEvent.getByLabel("offlineSlimmedPrimaryVertices",vertexs);
  iEvent.getByToken(theVtxTag,vertexs);
    
  edm::Handle<edm::ValueMap<bool> > veto_id_decisions;
  edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
  edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
  edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
  iEvent.getByToken(electronVetoIdMapToken_,veto_id_decisions);
  iEvent.getByToken(electronLooseIdMapToken_,loose_id_decisions);
  iEvent.getByToken(electronMediumIdMapToken_,medium_id_decisions);
  iEvent.getByToken(electronTightIdMapToken_,tight_id_decisions);

//*********************
  //Handle<edm::View<reco::GenParticle> > genParticles;
  edm::Handle<edm::ValueMap<bool> > medium_id_decisions2;
  edm::Handle<edm::ValueMap<bool> > tight_id_decisions2; 
  iEvent.getByToken(eleMediumIdMapToken_,medium_id_decisions2);
  iEvent.getByToken(eleTightIdMapToken_,tight_id_decisions2);
  edm::Handle<edm::ValueMap<float> > mvaValues;
  edm::Handle<edm::ValueMap<int> > mvaCategories;
  iEvent.getByToken(mvaValuesMapToken_,mvaValues);
  iEvent.getByToken(mvaCategoriesMapToken_,mvaCategories);

//**********************


  // Output collection
  auto_ptr<pat::ElectronCollection> result( new pat::ElectronCollection() );

  //unsigned int i=0;
  //for (unsigned int i = 0; i< electronHandle->size(); ++i){
  for( View<pat::Electron>::const_iterator el = electrons->begin(); el != electrons->end(); el++){
    const unsigned int i = distance (electrons->begin(), el);
    const auto ele = electrons->ptrAt(i);

    //---Clone the pat::Electron
    //pat::Electron l(*((*electrons)[i].get()));
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

    bool isconversionveto=l.passConversionVeto();
    
    //-- Missing hit  
    int missinghit = l.gsfTrack()->hitPattern().numberOfHits(HitPattern::MISSING_INNER_HITS);

    //--- 3 charge assignement
    bool isGsfCtfScPixChargeConsistent=l.isGsfCtfScPixChargeConsistent();
    bool isGsfScPixChargeConsistent=l.isGsfScPixChargeConsistent();

    float IoEmIoP_ttH = 9e9;
    if(l.ecalEnergy()>0)
      IoEmIoP_ttH = (1.0/l.ecalEnergy() - l.eSuperClusterOverP()/l.ecalEnergy());

    //--- Embed user variables
    l.addUserFloat("PFChargedHadIso",PFChargedHadIso);
    l.addUserFloat("PFNeutralHadIso",PFNeutralHadIso);
    l.addUserFloat("PFPhotonIso",PFPhotonIso);
    l.addUserFloat("combRelIsoPF",combRelIsoPF);
    l.addUserFloat("rho",rho);
    l.addUserFloat("SIP",SIP);
    l.addUserFloat("dxy",dxy);
    l.addUserFloat("dz",dz);
    l.addUserFloat("rel_error_trackpt",rel_error_trackpt);
    l.addUserInt("isConversionVeto",(isconversionveto ? 1 : 0));
    //l.addUserFloat("HLTMatch", HLTMatch);
    l.addUserInt("missingHit", missinghit);
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
    const Ptr<pat::Electron> elPtr(electrons, el - electrons->begin() );
    int eleCUT=0;
    if((*veto_id_decisions)[ elPtr ])eleCUT |= 1 << 0;
    if((*loose_id_decisions)[ elPtr ])eleCUT |= 1 << 1;
    if((*medium_id_decisions)[ elPtr ])eleCUT |= 1 << 2;
    if((*tight_id_decisions)[ elPtr ])eleCUT |= 1 << 3;
    l.addUserInt("isCUT",eleCUT);

    int isEleID80 = (*medium_id_decisions2)[ele];
    int  isEleID90  = (*tight_id_decisions2)[ele];
    //std::cout<<(*medium_id_decisions2)[ele]<<isEleID90<<endl;
    float eleMVAvalue=(*mvaValues)[ele];
    l.addUserInt("isEleID80",isEleID80);
    l.addUserInt("isEleID90",isEleID90);
    l.addUserFloat("eleMVAvalue",eleMVAvalue);
    l.addUserFloat("BDT",eleMVAvalue); //I know, it's duplicated, but I don't want to change to change all the downstream code...
    bool isBDT = false;
    float afSCeta = fabs(fSCeta);
    if(afSCeta < 2.4){
      if (afSCeta <0.8){
        if (l.pt()>10) isBDT = (eleMVAvalue>0.913286);
        else isBDT = (eleMVAvalue>-0.083313);
      } else{
        if (l.pt()>10) isBDT = (eleMVAvalue>0.805013);
        else isBDT = (eleMVAvalue>-0.235222);
      }
    }else{
     if (l.pt()>10) isBDT = (eleMVAvalue>0.358969);
     else isBDT = (eleMVAvalue>-0.67099);   
    }
    l.addUserInt("isBDT",(isBDT ? 1 : 0));

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

//     MCHistoryTools mch(iEvent);
//     if (mch.isMC()) {
//       int MCParentCode = 0;
//       //      int MCParentCode = mch.getParentCode(&l); //FIXME: does not work on cmg
//       l.addUserFloat("MCParentCode",MCParentCode);
//     }

    //--- Check selection cut. Being done here, flags are not available; but this way we 
    //    avoid wasting time on rejected leptons.
    if (!cut(l)) continue;

    //--- Embed flags (ie flags specified in the "flags" pset)
    for(CutSet<pat::Electron>::const_iterator flag = flags.begin(); flag != flags.end(); ++flag) {
      l.addUserFloat(flag->first,int((*(flag->second))(l)));
    }

    result->push_back(l);
    //i++;
  }
  iEvent.put(result);
}


#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(EleFiller);

