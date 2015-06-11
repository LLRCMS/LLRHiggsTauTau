/** \class TauFiller
 *
 *  No description available.
 *
 *  $Date: 2013/05/24 15:42:42 $
 *  $Revision: 1.28 $
 *  \author G. Ortona (LLR)
 */

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <DataFormats/PatCandidates/interface/Tau.h>
//#include <DataFormats/TauReco/interface/PFTauDiscriminator.h>

//#include "DataFormats/VertexReco/interface/Vertex.h"
#include <LLRHiggsTauTau/NtupleProducer/interface/DaughterDataHelpers.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/CutSet.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/LeptonIsoHelper.h>
//#include "BDTId.h"

#include <vector>
#include <string>

using namespace edm;
using namespace std;
using namespace reco;

//bool recomputeBDT = false;

class TauFiller : public edm::EDProducer {
 public:
  /// Constructor
  explicit TauFiller(const edm::ParameterSet&);
    
  /// Destructor
  ~TauFiller(){
    
  };  
  //ByIsolationMVA3oldDMwoLTraw
 private:
  virtual void beginJob(){};  
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};

  const edm::InputTag theCandidateTag;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > theGenTag ;
  //const edm::InputTag theDiscriminatorTag;
  const std::string theDiscriminatorTag;
  const StringCutObjectSelector<pat::Tau, true> cut;
  const CutSet<pat::Tau> flags;
  //BDTId* bdt;
};


TauFiller::TauFiller(const edm::ParameterSet& iConfig) :
  theCandidateTag(iConfig.getParameter<InputTag>("src")),
  theGenTag(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genCollection"))),
  theDiscriminatorTag(iConfig.getParameter<std::string>("discriminator")),
  cut(iConfig.getParameter<std::string>("cut")),
  flags(iConfig.getParameter<ParameterSet>("flags"))//, 
  //bdt(0)
{
  //if (recomputeBDT) bdt = new BDTId;
  produces<pat::TauCollection>();
}


void
TauFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{  

//read one PFTauDiscriminator (set discriminatorSrc_ in to an edm::InputTag before)

// loop over taus
//    for ( unsigned iTau = 0; iTau < taus->size(); ++iTau ) {
//        reco::PFTauRef tauCandidate(taus, iTau);
// check if tau candidate has passed discriminator
//        if( (*discriminator)[tauCandidate] > 0.5 ){
//        // do something with your candidate
//        }
//    }


  // Get leptons and discriminators
  //edm::Handle<pat::PFTauCollection> tauHandle;
  edm::Handle<pat::TauRefVector> tauHandle;
  iEvent.getByLabel(theCandidateTag, tauHandle);
  
  //edm::Handle<reco::PFTauDiscriminator> discriminator;
  //iEvent.getByLabel(theDiscriminatorTag, discriminator);

  //hpsPFTauDiscrimination
  
  edm::Handle<vector<Vertex> >  vertexs;
  iEvent.getByLabel("goodPrimaryVertices",vertexs);

  // Output collection
  auto_ptr<pat::TauCollection> result( new pat::TauCollection() );

  for (unsigned int i = 0; i< tauHandle->size(); ++i){
    //---Clone the pat::Tau
    pat::Tau l(*((*tauHandle)[i].get()));
    
    
    //--- PF ISO
    float PFChargedHadIso   = l.chargedHadronIso();
    float PFNeutralHadIso   = l.neutralHadronIso();
    float PFPhotonIso       = l.photonIso();

    float combRelIsoPF = LeptonIsoHelper::combRelIsoPF(l);

    // NOTE: tauID returns a float but I checked that the returned value is always 0 or 1 for "booleans" --> safe implicit cast to an int
    //float decayModeFindingOldDMs = l.tauID ("decayModeFindingOldDMs");
    int decayModeFindingOldDMs = l.tauID ("decayModeFinding");
    int decayModeFindingNewDMs = l.tauID ("decayModeFindingNewDMs");
    
    int byLooseCombinedIsolationDeltaBetaCorr3Hits = l.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
    int byMediumCombinedIsolationDeltaBetaCorr3Hits = l.tauID ("byMediumCombinedIsolationDeltaBetaCorr3Hits");
    int byTightCombinedIsolationDeltaBetaCorr3Hits = l.tauID ("byTightCombinedIsolationDeltaBetaCorr3Hits");
    
    float byCombinedIsolationDeltaBetaCorrRaw3Hits = l.tauID ("byCombinedIsolationDeltaBetaCorrRaw3Hits");
    float chargedIsoPtSum = l.tauID ("chargedIsoPtSum");
    float neutralIsoPtSum = l.tauID ("neutralIsoPtSum");
    float puCorrPtSum = l.tauID ("puCorrPtSum");
    
    int againstMuonLoose3 = l.tauID ("againstMuonLoose3");
    int againstMuonTight3 = l.tauID ("againstMuonTight3");
    
    int againstElectronVLooseMVA5 = l.tauID ("againstElectronVLooseMVA5");
    int againstElectronLooseMVA5 = l.tauID ("againstElectronLooseMVA5");
    int againstElectronMediumMVA5 = l.tauID ("againstElectronMediumMVA5");
    int againstElectronTightMVA5 = l.tauID ("againstElectronTightMVA5");
    int againstElectronVTightMVA5 = l.tauID ("againstElectronVTightMVA5");
    
    //Decay mode
    //int decayMode = -1;
    //int A = l.signalPFChargedHadrCands().size();
    //int B = l.signalPFGammaCands().size();
    //if(A==1&&B==0)decayMode=1;
    //else if(A==1&&B>0)decayMode=2;
    //else if (A==3&&B==0)decayMode=3;
    float tauid=l.tauID(theDiscriminatorTag);
    //printf("A, B, tau %d %d %f \n",A,B,tauid);

    //if(decayMode<0&&tauid==0)edm::LogWarning("TauFiller: Unrecognized decay mode");
    /*

    //--- SIP, dxy, dz
    float IP      = fabs(l.dB(pat::Electron::PV3D));
    float IPError = l.edB(pat::Electron::PV3D);
    float SIP     = IP/IPError;
    */
    
    float dxy = 999.;
    float dz  = 999.;
    if (vertexs->size()>0) {
      dxy = l.dxy();
      const Vertex* vertex = &(vertexs->front());          
      dz = l.vertex().z() - vertex[0].z();
      //For some reasons, the reference secondaryVertex() is empty EVEN if hasSecondaryVertex is true
      //To be asked to miniAOD people
      //if(l.hasSecondaryVertex()) {
        //dz  = l.secondaryVertex().get()->z()-vertex->z();      
      //}
    } 

   
    //--- Embed user variables
    l.addUserFloat("HPSDiscriminator",tauid);
    l.addUserFloat("decayMode",l.decayMode());
    l.addUserFloat("dxy",dxy);
    l.addUserFloat("dz",dz);
    l.addUserFloat("PFChargedHadIso",PFChargedHadIso);
    l.addUserFloat("PFNeutralHadIso",PFNeutralHadIso);
    l.addUserFloat("PFPhotonIso",PFPhotonIso);
    l.addUserFloat("combRelIsoPF",combRelIsoPF);
    l.addUserInt("decayModeFindingOldDMs", decayModeFindingOldDMs);
    l.addUserInt("decayModeFindingNewDMs", decayModeFindingNewDMs);
    l.addUserInt("byLooseCombinedIsolationDeltaBetaCorr3Hits", byLooseCombinedIsolationDeltaBetaCorr3Hits);
    l.addUserInt("byMediumCombinedIsolationDeltaBetaCorr3Hits", byMediumCombinedIsolationDeltaBetaCorr3Hits);
    l.addUserInt("byTightCombinedIsolationDeltaBetaCorr3Hits", byTightCombinedIsolationDeltaBetaCorr3Hits);
    l.addUserFloat("byCombinedIsolationDeltaBetaCorrRaw3Hits", byCombinedIsolationDeltaBetaCorrRaw3Hits);
    l.addUserFloat("chargedIsoPtSum", chargedIsoPtSum);
    l.addUserFloat("neutralIsoPtSum", neutralIsoPtSum);
    l.addUserFloat("puCorrPtSum", puCorrPtSum);
    l.addUserInt("againstMuonLoose3", againstMuonLoose3);
    l.addUserInt("againstMuonTight3", againstMuonTight3);
    l.addUserInt("againstElectronVLooseMVA5", againstElectronVLooseMVA5);
    l.addUserInt("againstElectronLooseMVA5", againstElectronLooseMVA5);
    l.addUserInt("againstElectronMediumMVA5", againstElectronMediumMVA5);
    l.addUserInt("againstElectronTightMVA5", againstElectronTightMVA5);
    l.addUserInt("againstElectronVTightMVA5", againstElectronVTightMVA5);
    
    

    //--- MC parent code 
    const reco::GenParticle* genL= l.genParticleRef().get();
    float px=0,py=0,pz=0,E=0,fromH=0;
    float pxHad=0, pyHad=0, pzHad=0, EHad=0; // hadronic gen tau
    int status=99999, id=99999;

    if(genL){
      px =genL->p4().Px();
      py =genL->p4().Py();
      pz =genL->p4().Pz();
      E =genL->p4().E();
      status =genL->status();
      id =genL->pdgId();

      //cout << "Tau filler: " << i << " [px, id] = " << l.px() << " , " << l.pdgId() << " | (px, py, pz, e) " << px << " " << py << " " << pz << " " << E << " | ID: " << genL->pdgId() << " | status: " << genL->status() << endl;

      // build hadronic gen tau (all visible sons)
      for (unsigned int iDau = 0; iDau < genL->numberOfDaughters(); iDau++)
      {
          const Candidate * Dau = genL->daughter( iDau );
          int dauId = Dau->pdgId();
          if (abs(dauId) != 12 && abs(dauId) != 14 && abs(dauId) != 16)
          {
              pxHad += Dau->p4().Px();
              pyHad += Dau->p4().Py();
              pzHad += Dau->p4().Pz();
              EHad += Dau->p4().E();
          }
      }
      
      //search if it comes from H
      Handle<edm::View<reco::GenParticle> > genHandle;
      iEvent.getByToken(theGenTag, genHandle);
      for(unsigned int ipruned = 0; ipruned< genHandle->size(); ++ipruned){
        int pdgmot = (&(*genHandle)[ipruned])->pdgId();
        if(abs(pdgmot)==25){
          if(userdatahelpers::isAncestor(&(*genHandle)[ipruned],genL)){
            fromH=1;
            break;
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

    l.addUserFloat("genHadPx",px);
    l.addUserFloat("genHadPy",py);
    l.addUserFloat("genHadPz",pz);
    l.addUserFloat("genHadE",E);

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
    for(CutSet<pat::Tau>::const_iterator flag = flags.begin(); flag != flags.end(); ++flag) {
      l.addUserFloat(flag->first,int((*(flag->second))(l)));
    }
    
    result->push_back(l);
  }
  iEvent.put(result);
}


#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(TauFiller);

