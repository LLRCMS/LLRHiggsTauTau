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
  int lep_setup;
  const StringCutObjectSelector<pat::Electron, true> cut;
  const CutSet<pat::Electron> flags;
  //EGammaMvaEleEstimatorCSA14* myMVATrig;
  edm::EDGetTokenT<edm::View<pat::Electron> > electronCollectionToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > miniRelIsoChargedNanoAODToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > miniRelIsoAllNanoAODToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > PFRelIsoChargedNanoAODToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > PFRelIsoAllNanoAODToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > PFRelIsoAll04NanoAODToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > ptRelNanoAODToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > ptRatioNanoAODToken_;  
  edm::EDGetTokenT<edm::ValueMap<float> > jetNDauChargedMVASelNanoAODToken_;
  
  //BDTId* bdt;
  };


EleFiller::EleFiller(const edm::ParameterSet& iConfig) :
  //theCandidateTag(iConfig.getParameter<InputTag>("src")),
  theGenTag(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genCollection"))),
  theRhoTag(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoCollection"))),
  theVtxTag(consumes<vector<Vertex>>(iConfig.getParameter<edm::InputTag>("vtxCollection"))),
  sampleType(iConfig.getParameter<int>("sampleType")),
  setup(iConfig.getParameter<int>("setup")),
  lep_setup(iConfig.getParameter<int>("lep_setup")),
  cut(iConfig.getParameter<std::string>("cut")),
  flags(iConfig.getParameter<ParameterSet>("flags")),//,
  //myMVATrig(0),
  electronCollectionToken_(consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("src"))),
  miniRelIsoChargedNanoAODToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("miniRelIsoChgCollection"))),
  miniRelIsoAllNanoAODToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("miniRelIsoAllCollection"))),
  PFRelIsoChargedNanoAODToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("PFRelIsoChgCollection"))),
  PFRelIsoAllNanoAODToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("PFRelIsoAllCollection"))),
  PFRelIsoAll04NanoAODToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("PFRelIsoAll04Collection"))),
  ptRelNanoAODToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("ptRelCollection"))),
  ptRatioNanoAODToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("ptRatioCollection"))),
  jetNDauChargedMVASelNanoAODToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("jetNDauChargedMVASelCollection")))
 
 //****************
  //bdt(0)
{
  
  produces<pat::ElectronCollection>();
}

  void EleFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
  {  
  edm::Handle<edm::View<pat::Electron> > electrons;
  iEvent.getByToken(electronCollectionToken_, electrons);

  //InputTag theRhoTag = LeptonIsoHelper::getEleRhoTag(sampleType,setup);
  edm::Handle<double> rhoHandle;
  iEvent.getByToken(theRhoTag, rhoHandle); 
  double rho = *rhoHandle;

  edm::Handle<vector<Vertex> >  vertexs;
  //iEvent.getByLabel("offlineSlimmedPrimaryVertices",vertexs);
  iEvent.getByToken(theVtxTag,vertexs);
  
  edm::Handle<edm::ValueMap<float> > miniRelIsoChg_2;
  edm::Handle<edm::ValueMap<float> > miniRelIsoAll_2;
  edm::Handle<edm::ValueMap<float> > PFRelIsoChg_2;
  edm::Handle<edm::ValueMap<float> > PFRelIsoAll_2;
  edm::Handle<edm::ValueMap<float> > PFRelIsoAll04_2;
  iEvent.getByToken(miniRelIsoChargedNanoAODToken_,miniRelIsoChg_2);
  iEvent.getByToken(miniRelIsoAllNanoAODToken_,miniRelIsoAll_2);
  iEvent.getByToken(PFRelIsoChargedNanoAODToken_,PFRelIsoChg_2);
  iEvent.getByToken(PFRelIsoAllNanoAODToken_,PFRelIsoAll_2);
  iEvent.getByToken(PFRelIsoAll04NanoAODToken_,PFRelIsoAll04_2);

  edm::Handle<edm::ValueMap<float> > ptRel_2;
  edm::Handle<edm::ValueMap<float> > ptRatio_2;
  edm::Handle<edm::ValueMap<float> > jetNDauChargedMVASel_2;
  iEvent.getByToken(ptRelNanoAODToken_,ptRel_2);
  iEvent.getByToken(ptRatioNanoAODToken_,ptRatio_2);
  iEvent.getByToken(jetNDauChargedMVASelNanoAODToken_,jetNDauChargedMVASel_2);

//**********************


  // Output collection
  //auto_ptr<pat::ElectronCollection> result( new pat::ElectronCollection() );
  std::unique_ptr<pat::ElectronCollection> result( new pat::ElectronCollection() );

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

    /*bool passVetoId = l.userInt("cutBasedElectronID-Fall17-94X-V2-veto");
    bool passLooseId = l.userInt("cutBasedElectronID-Fall17-94X-V2-loose");
    bool passMediumId = l.userInt("cutBasedElectronID-Fall17-94X-V2-medium");
    bool passTightId = l.userInt("cutBasedElectronID-Fall17-94X-V2-tight");*/
    
    bool passMvaIsowp80Id = l.electronID("mvaEleID-Fall17-iso-V2-wp80");
    bool passMvanonIsowp80Id = l.electronID("mvaEleID-Fall17-noIso-V2-wp80");
    bool passMvaIsowp90Id = l.electronID("mvaEleID-Fall17-iso-V2-wp90");
    bool passMvanonIsowp90Id = l.electronID("mvaEleID-Fall17-noIso-V2-wp90");
    bool passMvaIsowpLooseId = l.electronID("mvaEleID-Fall17-iso-V2-wpLoose");
    bool passMvanonIsowpLooseId = l.electronID("mvaEleID-Fall17-noIso-V2-wpLoose");
 
    float mvaValue_Iso = l.userFloat("ElectronMVAEstimatorRun2Fall17IsoV2Values");
    float mvaValue_noIso = l.userFloat("ElectronMVAEstimatorRun2Fall17NoIsoV2Values");
    float mvaCategory_Iso = l.userInt("ElectronMVAEstimatorRun2Fall17IsoV2Categories");
    float mvaCategory_noIso = l.userInt("ElectronMVAEstimatorRun2Fall17NoIsoV2Categories");

    bool passMvaHZZwpLooseId = l.electronID("mvaEleID-Spring16-HZZ-V1-wpLoose");
    float mvaValue_HZZ = l.userFloat("ElectronMVAEstimatorRun2Spring16HZZV1Values");
    float mvaCategory_HZZ = l.userInt("ElectronMVAEstimatorRun2Spring16HZZV1Categories");

    bool isconversionveto=l.passConversionVeto();
    
    //-- Missing hit  
    //int missinghit = l.gsfTrack()->hitPattern().numberOfHits(HitPattern::MISSING_INNER_HITS);
    int missinghit = l.gsfTrack()->hitPattern().numberOfAllHits(HitPattern::MISSING_INNER_HITS);
    int missinglosthit = l.gsfTrack()->hitPattern().numberOfLostHits(HitPattern::MISSING_INNER_HITS);
	  
    //--- 3 charge assignement
    bool isGsfCtfScPixChargeConsistent=l.isGsfCtfScPixChargeConsistent();
    bool isGsfScPixChargeConsistent=l.isGsfScPixChargeConsistent();

    float IoEmIoP_ttH = 9e9;
    if(l.ecalEnergy()>0)
      IoEmIoP_ttH = (1.0/l.ecalEnergy() - l.eSuperClusterOverP()/l.ecalEnergy());

    //--- scale and smearing corrections and systematics
    float ecalEnergyErrPreCorr = l.userFloat("ecalEnergyErrPreCorr");
    float ecalEnergyErrPostCorr = l.userFloat("ecalEnergyErrPostCorr");
    float ecalTrkEnergyPreCorr = l.userFloat("ecalTrkEnergyPreCorr");
    float ecalTrkEnergyErrPreCorr = l.userFloat("ecalTrkEnergyErrPreCorr");
    float ecalTrkEnergyPostCorr = l.userFloat("ecalTrkEnergyPostCorr");
    float ecalTrkEnergyErrPostCorr = l.userFloat("ecalTrkEnergyErrPostCorr");
    float energyScaleValue = l.userFloat("energyScaleValue");
    float energySigmaValue = l.userFloat("energySigmaValue");
    float energySmearNrSigma = l.userFloat("energySmearNrSigma");
    float energyScaleUp = l.userFloat("energyScaleUp");
    float energyScaleDown = l.userFloat("energyScaleDown");
    float energyScaleStatUp = l.userFloat("energyScaleStatUp");
    float energyScaleStatDown = l.userFloat("energyScaleStatDown");
    float energyScaleSystUp = l.userFloat("energyScaleSystUp");
    float energyScaleSystDown = l.userFloat("energyScaleSystDown");
    float energyScaleGainUp = l.userFloat("energyScaleGainUp");
    float energyScaleGainDown = l.userFloat("energyScaleGainDown");
    float energyScaleEtUp = l.userFloat("energyScaleEtUp");
    float energyScaleEtDown = l.userFloat("energyScaleEtDown");
    float energySigmaUp = l.userFloat("energySigmaUp");
    float energySigmaDown = l.userFloat("energySigmaDown");
    float energySigmaPhiUp = l.userFloat("energySigmaPhiUp");
    float energySigmaPhiDown = l.userFloat("energySigmaPhiDown");
    float energySigmaRhoUp = l.userFloat("energySigmaRhoUp");
    float energySigmaRhoDown = l.userFloat("energySigmaRhoDown");

    float corr_ele = ecalTrkEnergyPostCorr/ecalTrkEnergyPreCorr; 
    float pt_corr = corr_ele*(l.pt());
    float E_corr = corr_ele*(l.energy());
    //std::cout<<"Eta "<<l.eta()<<",post: "<<corr_ele_ecalTrkEnergyPostCorr<<",pre: "<<corr_ele_ecalTrkEnergyPreCorr<<", corr_ele "<<corr_ele<<",pt_corr "<<pt_corr<<std::endl;

    //--- Embed user variables
    l.addUserFloat("PFChargedHadIso",PFChargedHadIso);
    l.addUserFloat("PFNeutralHadIso",PFNeutralHadIso);
    l.addUserFloat("PFPhotonIso",PFPhotonIso);
    l.addUserFloat("combRelIsoPF",combRelIsoPF);
    l.addUserFloat("rho",rho);
    l.addUserFloat("SIP",SIP);
    l.addUserFloat("dxy",dxy);
    l.addUserFloat("dz",dz);
    l.addUserFloat("passMvaIsowp80Id",passMvaIsowp80Id);
    l.addUserFloat("passMvanonIsowp80Id",passMvanonIsowp80Id);
    l.addUserFloat("passMvaIsowp90Id",passMvaIsowp90Id);
    l.addUserFloat("passMvanonIsowp90Id",passMvanonIsowp90Id);
    l.addUserFloat("passMvaIsowpLooseId",passMvaIsowpLooseId);
    l.addUserFloat("passMvanonIsowpLooseId",passMvanonIsowpLooseId);
    l.addUserFloat("mvaValue_Iso",mvaValue_Iso);
    l.addUserFloat("mvaValue_noIso",mvaValue_noIso);
    l.addUserInt("mvaCategory_Iso",mvaCategory_Iso);
    l.addUserInt("mvaCategory_noIso",mvaCategory_noIso);
    l.addUserFloat("passMvaHZZwpLooseId",passMvaHZZwpLooseId);
    l.addUserFloat("mvaValue_HZZ",mvaValue_HZZ);
    l.addUserInt("mvaCategory_HZZ",mvaCategory_HZZ);
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
    l.addUserFloat("ecalEnergyErrPreCorr",ecalEnergyErrPreCorr);
    l.addUserFloat("ecalEnergyErrPostCorr",ecalEnergyErrPostCorr);
    l.addUserFloat("ecalTrkEnergyPreCorr",ecalTrkEnergyPreCorr);
    l.addUserFloat("ecalTrkEnergyErrPreCorr",ecalTrkEnergyErrPreCorr);
    l.addUserFloat("ecalTrkEnergyPostCorr",ecalTrkEnergyPostCorr);
    l.addUserFloat("ecalTrkEnergyErrPostCorr",ecalTrkEnergyErrPostCorr);
    l.addUserFloat("energyScaleValue",energyScaleValue);
    l.addUserFloat("energySigmaValue",energySigmaValue);
    l.addUserFloat("energySmearNrSigma",energySmearNrSigma);
    l.addUserFloat("energyScaleUp",energyScaleUp);
    l.addUserFloat("energyScaleDown",energyScaleDown);
    l.addUserFloat("energyScaleStatUp",energyScaleStatUp);
    l.addUserFloat("energyScaleStatDown",energyScaleStatDown);
    l.addUserFloat("energyScaleSystUp",energyScaleSystUp);
    l.addUserFloat("energyScaleSystDown",energyScaleSystDown);
    l.addUserFloat("energyScaleGainUp",energyScaleGainUp);
    l.addUserFloat("energyScaleGainDown",energyScaleGainDown);
    l.addUserFloat("energyScaleEtUp",energyScaleEtUp);
    l.addUserFloat("energyScaleEtDown",energyScaleEtDown);
    l.addUserFloat("energySigmaUp",energySigmaUp);
    l.addUserFloat("energySigmaDown",energySigmaDown);
    l.addUserFloat("energySigmaPhiUp",energySigmaPhiUp);
    l.addUserFloat("energySigmaPhiDown",energySigmaPhiDown);
    l.addUserFloat("energySigmaRhoUp",energySigmaRhoUp);
    l.addUserFloat("energySigmaRhoDown",energySigmaRhoDown);
    l.addUserFloat("corr_ele",corr_ele);
    l.addUserFloat("pt_corr",pt_corr);
    l.addUserFloat("E_corr",E_corr);
    l.addUserInt("isEB", int(l.isEB()));
    const Ptr<pat::Electron> elPtr(electrons, el - electrons->begin() );
    float eleMVAvalue = mvaValue_Iso;
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

    float miniRelIsoChg2Value = (*miniRelIsoChg_2)[ele];
    float miniRelIsoAll2Value = (*miniRelIsoAll_2)[ele];
    float PFRelIsoChg2Value = (*PFRelIsoChg_2)[ele];
    float PFRelIsoAll2Value = (*PFRelIsoAll_2)[ele];
    float PFRelIsoAll042Value = (*PFRelIsoAll04_2)[ele];
    l.addUserFloat("miniRelIsoChg2Value",miniRelIsoChg2Value);
    l.addUserFloat("miniRelIsoAll2Value",miniRelIsoAll2Value);
    l.addUserFloat("PFRelIsoChg2Value",PFRelIsoChg2Value);
    l.addUserFloat("PFRelIsoAll2Value",PFRelIsoAll2Value);
    l.addUserFloat("PFRelIsoAll042Value",PFRelIsoAll042Value);
	  
    float ptRel2Value = (*ptRel_2)[ele];
    float ptRatio2Value = (*ptRatio_2)[ele];
    float jetNDauChargedMVASel2Value = (*jetNDauChargedMVASel_2)[ele];
    l.addUserFloat("ptRel2Value",ptRel2Value);
    l.addUserFloat("ptRatio2Value",ptRatio2Value);
    l.addUserFloat("jetNDauChargedMVASel2Value",jetNDauChargedMVASel2Value);

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
  //iEvent.put(result);
  iEvent.put(std::move(result));
}


#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(EleFiller);

