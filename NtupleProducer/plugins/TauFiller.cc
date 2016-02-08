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

  edm::EDGetTokenT<pat::TauRefVector> theCandidateTag;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > theGenTag ;
  edm::EDGetTokenT<vector<Vertex> > theVtxTag ;
  const std::string theDiscriminatorTag;
  const StringCutObjectSelector<pat::Tau, true> cut;
  const CutSet<pat::Tau> flags;
  const std::string NominalUpOrDown;
};


TauFiller::TauFiller(const edm::ParameterSet& iConfig) :
  theCandidateTag(consumes<pat::TauRefVector>(iConfig.getParameter<InputTag>("src"))),
  theGenTag(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genCollection"))),
  theVtxTag(consumes<vector<Vertex>>(iConfig.getParameter<edm::InputTag>("vtxCollection"))),
  theDiscriminatorTag(iConfig.getParameter<std::string>("discriminator")),
  cut(iConfig.getParameter<std::string>("cut")),
  flags(iConfig.getParameter<ParameterSet>("flags")), 
  NominalUpOrDown(iConfig.getParameter<std::string>("NominalUpOrDown"))//,
{
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
  iEvent.getByToken(theCandidateTag, tauHandle);
  
  //edm::Handle<reco::PFTauDiscriminator> discriminator;
  //iEvent.getByLabel(theDiscriminatorTag, discriminator);

  //hpsPFTauDiscrimination
  
  edm::Handle<vector<Vertex> >  vertexs;
  //iEvent.getByLabel("goodPrimaryVertices",vertexs);
  iEvent.getByToken(theVtxTag, vertexs);

  // Output collection
  auto_ptr<pat::TauCollection> result( new pat::TauCollection() );

  for (unsigned int i = 0; i< tauHandle->size(); ++i){
    //---Clone the pat::Tau
    pat::Tau l(*((*tauHandle)[i].get()));
    
    double Shift = 1.;
    if(NominalUpOrDown=="Nominal") Shift = 1.;
    if(NominalUpOrDown=="Up") Shift = 1.03;
    if(NominalUpOrDown=="Down") Shift = 0.97;
    // cout<<"tau #"<<i<<endl;
    // cout<<"NominalUpOrDown = "<<NominalUpOrDown<<endl;
    double shiftP = 1.;
    double shiftMass = 1.;
    int isTESShifted = 0;
    if ( l.genJet() && deltaR(l.p4(), l.genJet()->p4()) < 0.5 && l.genJet()->pt() > 8. ) {

      // cout<<"I have pT = "<<l.pt()<<" and will be shifted"<<endl;
      isTESShifted = 1;
      // cout<<"l.decayMode() = "<< l.decayMode()<<endl;
      // cout<<"l.signalPFChargedHadrCands().size() = "<<(l.signalPFChargedHadrCands()).size()<<endl;
      // cout<<"l.signalPFGammaCands().size() = "<<(l.signalPFGammaCands()).size()<<endl;

      if(l.decayMode()>=1 && l.decayMode()<10){
      // if((l.signalPFChargedHadrCands()).size()==1 && (l.signalPFGammaCands()).size()>0){
	//cout<<"1-prong + pi0"<<endl;
	shiftP = Shift;
	shiftMass = Shift;
      }

      else if(l.decayMode()==0){
      // else if((l.signalPFChargedHadrCands()).size()==1 && (l.signalPFGammaCands()).size()==0){
	//cout<<"1-prong"<<endl;
	shiftP = Shift;
	shiftMass = 1.;
      }
      else if(l.decayMode()==10){
      // else if((l.signalPFChargedHadrCands()).size()==3) {
	//cout<<"3-prong"<<endl;
	shiftP = Shift;
	shiftMass = Shift;
      }
      else isTESShifted = 0;
    }
    
    double pxS = l.px()*shiftP;
    double pyS = l.py()*shiftP;
    double pzS = l.pz()*shiftP;
    double massS = l.mass()*shiftMass;
    double enS = TMath::Sqrt(pxS*pxS + pyS*pyS + pzS*pzS + massS*massS);
    math::XYZTLorentzVectorD p4S( pxS, pyS, pzS, enS );

    //Dump to check the T-ES
    // if((l.decayMode()==0 || l.decayMode()==1 || l.decayMode()==2 || l.decayMode()==10) && l.genJet() && deltaR(l.p4(), l.genJet()->p4()) < 0.5 && l.genJet()->pt() > 8. && NominalUpOrDown=="Up")
    //   {
    // 	cout<<"NominalUpOrDown = "<<NominalUpOrDown<<endl;
    // 	cout<<"l.decayMode() = "<<l.decayMode()<<endl;
    // 	cout<<"pt before = "<<l.pt()<<endl;
    // 	cout<<"pt after = "<<p4S.Pt()<<endl;
    //   }
    l.setP4( p4S );
    // cout<<"after, I have pT = "<<l.pt()<<" and my decayMode is = "<<l.decayMode()<<endl;
    // if((l.decayMode()==0 || l.decayMode()==1 || l.decayMode()==2 || l.decayMode()==10) && l.genJet() && deltaR(l.p4(), l.genJet()->p4()) < 0.5 && l.genJet()->pt() > 8. && NominalUpOrDown=="Up") cout<<"pt after in pat = "<<l.pt()<<endl;
    
    
    
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

    int byLooseIsolationMVArun2v1DBoldDMwLT = l.tauID("byLooseIsolationMVArun2v1DBoldDMwLT");
    int byMediumIsolationMVArun2v1DBoldDMwLT = l.tauID("byMediumIsolationMVArun2v1DBoldDMwLT");
    int byTightIsolationMVArun2v1DBoldDMwLT = l.tauID("byTightIsolationMVArun2v1DBoldDMwLT");
    int byVTightIsolationMVArun2v1DBoldDMwLT = l.tauID("byVTightIsolationMVArun2v1DBoldDMwLT");
    int byLooseIsolationMVArun2v1DBnewDMwLT = l.tauID("byLooseIsolationMVArun2v1DBnewDMwLT");
    int byMediumIsolationMVArun2v1DBnewDMwLT = l.tauID("byMediumIsolationMVArun2v1DBnewDMwLT");
    int byTightIsolationMVArun2v1DBnewDMwLT = l.tauID("byTightIsolationMVArun2v1DBnewDMwLT");
    int byVTightIsolationMVArun2v1DBnewDMwLT = l.tauID("byVTightIsolationMVArun2v1DBnewDMwLT");
    int byLooseCombinedIsolationDeltaBetaCorr3HitsdR03 = l.tauID("byLooseCombinedIsolationDeltaBetaCorr3HitsdR03");
    int byMediumCombinedIsolationDeltaBetaCorr3HitsdR03 = l.tauID("byMediumCombinedIsolationDeltaBetaCorr3HitsdR03");
    int byTightCombinedIsolationDeltaBetaCorr3HitsdR03 = l.tauID("byTightCombinedIsolationDeltaBetaCorr3HitsdR03");
    int byLooseIsolationMVArun2v1DBdR03oldDMwLT = l.tauID("byLooseIsolationMVArun2v1DBdR03oldDMwLT");
    int byMediumIsolationMVArun2v1DBdR03oldDMwLT = l.tauID("byMediumIsolationMVArun2v1DBdR03oldDMwLT");
    int byTightIsolationMVArun2v1DBdR03oldDMwLT = l.tauID("byTightIsolationMVArun2v1DBdR03oldDMwLT");
    int byVTightIsolationMVArun2v1DBdR03oldDMwLT = l.tauID("byVTightIsolationMVArun2v1DBdR03oldDMwLT");

    float byCombinedIsolationDeltaBetaCorrRaw3Hits = l.tauID ("byCombinedIsolationDeltaBetaCorrRaw3Hits");
    float chargedIsoPtSum = l.tauID ("chargedIsoPtSum");
    float neutralIsoPtSum = l.tauID ("neutralIsoPtSum");
    float puCorrPtSum = l.tauID ("puCorrPtSum");
    
    int againstMuonLoose3 = l.tauID ("againstMuonLoose3");
    int againstMuonTight3 = l.tauID ("againstMuonTight3");
    
    int againstElectronVLooseMVA6 = l.tauID ("againstElectronVLooseMVA6");
    int againstElectronLooseMVA6 = l.tauID ("againstElectronLooseMVA6");
    int againstElectronMediumMVA6 = l.tauID ("againstElectronMediumMVA6");
    int againstElectronTightMVA6 = l.tauID ("againstElectronTightMVA6");
    int againstElectronVTightMVA6 = l.tauID ("againstElectronVTightMVA6");

    int numChargedParticlesSignalCone = l.signalChargedHadrCands().size();
    int numNeutralHadronsSignalCone = l.signalNeutrHadrCands().size();
    int numPhotonsSignalCone = l.signalGammaCands().size();
    int numParticlesSignalCone = l.signalCands().size();
    int numChargedParticlesIsoCone = l.isolationChargedHadrCands().size();
    int numNeutralHadronsIsoCone = l.isolationNeutrHadrCands().size();
    int numPhotonsIsoCone = l.isolationGammaCands().size();
    int numParticlesIsoCone = l.isolationCands().size();
    float leadChargedParticlePt=l.leadCand()->pt();
    float trackRefPt = (l.leadChargedHadrCand().isNonnull() ? l.leadChargedHadrCand()->pt() : 0.);

    
    //Decay mode
    //int decayMode = -1;
    //int A = l.signalPFChargedHadrCands().size();
    //int B = l.signalPFGammaCands().size();
    //if(A==1&&B==0)decayMode=1;
    //else if(A==1&&B>0)decayMode=2;
    //else if (A==3&&B==0)decayMode=3;
    float tauid = (l.isTauIDAvailable(theDiscriminatorTag) ? l.tauID(theDiscriminatorTag) : -999);
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
      //dxy = l.dxy();
      //const Vertex* vertex = &(vertexs->front());          
      //dz = l.vertex().z() - vertex[0].z();

      pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(l.leadChargedHadrCand().get());
      dz=packedLeadTauCand->dz();
      dxy=packedLeadTauCand->dxy();

      //For some reasons, the reference secondaryVertex() is empty EVEN if hasSecondaryVertex is true
      //To be asked to miniAOD people
      //if(l.hasSecondaryVertex()) {
        //dz  = l.secondaryVertex().get()->z()-vertex->z();      
      //}
    } 
   
    //--- Embed user variables
    l.addUserInt("isTESShifted",isTESShifted);
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

    l.addUserInt("byLooseIsolationMVArun2v1DBoldDMwLT",byLooseIsolationMVArun2v1DBoldDMwLT);
    l.addUserInt("byMediumIsolationMVArun2v1DBoldDMwLT",byMediumIsolationMVArun2v1DBoldDMwLT);
    l.addUserInt("byTightIsolationMVArun2v1DBoldDMwLT",byTightIsolationMVArun2v1DBoldDMwLT);
    l.addUserInt("byVTightIsolationMVArun2v1DBoldDMwLT",byVTightIsolationMVArun2v1DBoldDMwLT);
    l.addUserInt("byLooseIsolationMVArun2v1DBnewDMwLT",byLooseIsolationMVArun2v1DBnewDMwLT);
    l.addUserInt("byMediumIsolationMVArun2v1DBnewDMwLT",byMediumIsolationMVArun2v1DBnewDMwLT);
    l.addUserInt("byTightIsolationMVArun2v1DBnewDMwLT",byTightIsolationMVArun2v1DBnewDMwLT);
    l.addUserInt("byVTightIsolationMVArun2v1DBnewDMwLT",byVTightIsolationMVArun2v1DBnewDMwLT);
    l.addUserInt("byLooseCombinedIsolationDeltaBetaCorr3HitsdR03",byLooseCombinedIsolationDeltaBetaCorr3HitsdR03);
    l.addUserInt("byMediumCombinedIsolationDeltaBetaCorr3HitsdR03",byMediumCombinedIsolationDeltaBetaCorr3HitsdR03);
    l.addUserInt("byTightCombinedIsolationDeltaBetaCorr3HitsdR03",byTightCombinedIsolationDeltaBetaCorr3HitsdR03);
    l.addUserInt("byLooseIsolationMVArun2v1DBdR03oldDMwLT",byLooseIsolationMVArun2v1DBdR03oldDMwLT);
    l.addUserInt("byMediumIsolationMVArun2v1DBdR03oldDMwLT",byMediumIsolationMVArun2v1DBdR03oldDMwLT);
    l.addUserInt("byTightIsolationMVArun2v1DBdR03oldDMwLT",byTightIsolationMVArun2v1DBdR03oldDMwLT);
    l.addUserInt("byVTightIsolationMVArun2v1DBdR03oldDMwLT",byVTightIsolationMVArun2v1DBdR03oldDMwLT);

    l.addUserFloat("byCombinedIsolationDeltaBetaCorrRaw3Hits", byCombinedIsolationDeltaBetaCorrRaw3Hits);
    if (l.isTauIDAvailable("byIsolationMVA3oldDMwoLTraw")) l.addUserFloat("byIsolationMVA3oldDMwoLTraw",l.tauID("byIsolationMVA3oldDMwoLTraw"));
    else                                                   l.addUserFloat("byIsolationMVA3oldDMwoLTraw",-999.);    
    l.addUserFloat("byIsolationMVA3oldDMwLTraw",l.tauID("byIsolationMVA3oldDMwLTraw"));
    if (l.isTauIDAvailable("byIsolationMVA3newDMwoLTraw")) l.addUserFloat("byIsolationMVA3newDMwoLTraw",l.tauID("byIsolationMVA3newDMwoLTraw"));
    else                                                   l.addUserFloat("byIsolationMVA3newDMwoLTraw",-999.);    
    if (l.isTauIDAvailable("againstElectronMVA5category")) l.addUserInt("againstElectronMVA5category",l.tauID("againstElectronMVA5category"));
    else l.addUserInt("againstElectronMVA5category",0);
    if (l.isTauIDAvailable("againstElectronMVA5raw")) l.addUserFloat("againstElectronMVA5raw",l.tauID("againstElectronMVA5raw"));
    else l.addUserFloat("againstElectronMVA5raw",-999);
    if (l.isTauIDAvailable("byLooseIsolationMVA3newDMwLT")) l.addUserInt("byLooseIsolationMVA3newDMwLT",l.tauID("byLooseIsolationMVA3newDMwLT"));
    else l.addUserInt("byLooseIsolationMVA3newDMwLT",0);
    if (l.isTauIDAvailable("byLooseIsolationMVA3oldDMwLT")) l.addUserInt("byLooseIsolationMVA3oldDMwLT",l.tauID("byLooseIsolationMVA3oldDMwLT"));
    else l.addUserFloat("byLooseIsolationMVA3oldDMwLT",0);
    if (l.isTauIDAvailable("byLoosePileupWeightedIsolation3Hits")) l.addUserInt("byLoosePileupWeightedIsolation3Hits",l.tauID("byLoosePileupWeightedIsolation3Hits"));
    else l.addUserInt("byLoosePileupWeightedIsolation3Hits",0);
    if (l.isTauIDAvailable("byMediumIsolationMVA3newDMwLT")) l.addUserInt("byMediumIsolationMVA3newDMwLT",l.tauID("byMediumIsolationMVA3newDMwLT"));
    else l.addUserInt("byMediumIsolationMVA3newDMwLT",0);
    if (l.isTauIDAvailable("byMediumIsolationMVA3oldDMwLT")) l.addUserInt("byMediumIsolationMVA3oldDMwLT",l.tauID("byMediumIsolationMVA3oldDMwLT"));
    else l.addUserInt("byMediumIsolationMVA3oldDMwLT",0);
    if (l.isTauIDAvailable("byMediumPileupWeightedIsolation3Hits")) l.addUserInt("byMediumPileupWeightedIsolation3Hits",l.tauID("byMediumPileupWeightedIsolation3Hits"));
    else l.addUserInt("byMediumPileupWeightedIsolation3Hits",0);
    if (l.isTauIDAvailable("byPhotonPtSumOutsideSignalCone")) l.addUserFloat("byPhotonPtSumOutsideSignalCone",l.tauID("byPhotonPtSumOutsideSignalCone"));
    else l.addUserFloat("byPhotonPtSumOutsideSignalCone",-999);
    if (l.isTauIDAvailable("byPileupWeightedIsolationRaw3Hits")) l.addUserFloat("byPileupWeightedIsolationRaw3Hits",l.tauID("byPileupWeightedIsolationRaw3Hits"));
    else l.addUserFloat("byPileupWeightedIsolationRaw3Hits",-999);
    if (l.isTauIDAvailable("byTightIsolationMVA3newDMwLT")) l.addUserInt("byTightIsolationMVA3newDMwLT",l.tauID("byTightIsolationMVA3newDMwLT"));
    else l.addUserInt("byTightIsolationMVA3newDMwLT",-999);
    if (l.isTauIDAvailable("byTightIsolationMVA3oldDMwLT")) l.addUserInt("byTightIsolationMVA3oldDMwLT",l.tauID("byTightIsolationMVA3oldDMwLT"));
    else l.addUserInt("byTightIsolationMVA3oldDMwLT",-999);
    if (l.isTauIDAvailable("byTightPileupWeightedIsolation3Hits")) l.addUserInt("byTightPileupWeightedIsolation3Hits",l.tauID("byTightPileupWeightedIsolation3Hits"));
    else l.addUserInt("byTightPileupWeightedIsolation3Hits",-999);
    if (l.isTauIDAvailable("byVLooseIsolationMVA3newDMwLT")) l.addUserInt("byVLooseIsolationMVA3newDMwLT",l.tauID("byVLooseIsolationMVA3newDMwLT"));
    else l.addUserInt("byVLooseIsolationMVA3newDMwLT",0);
    if (l.isTauIDAvailable("byVLooseIsolationMVA3oldDMwLT")) l.addUserInt("byVLooseIsolationMVA3oldDMwLT",l.tauID("byVLooseIsolationMVA3oldDMwLT"));
    else l.addUserInt("byVLooseIsolationMVA3oldDMwLT",0);
    if (l.isTauIDAvailable("byVTightIsolationMVA3newDMwLT")) l.addUserInt("byVTightIsolationMVA3newDMwLT",l.tauID("byVTightIsolationMVA3newDMwLT"));
    else l.addUserInt("byVTightIsolationMVA3newDMwLT",0);
    if (l.isTauIDAvailable("byVTightIsolationMVA3oldDMwLT")) l.addUserInt("byVTightIsolationMVA3oldDMwLT",l.tauID("byVTightIsolationMVA3oldDMwLT"));
    else l.addUserInt("byVTightIsolationMVA3oldDMwLT",0);
    if (l.isTauIDAvailable("byVVTightIsolationMVA3newDMwLT")) l.addUserInt("byVVTightIsolationMVA3newDMwLT",l.tauID("byVVTightIsolationMVA3newDMwLT"));
    else l.addUserInt("byVVTightIsolationMVA3newDMwLT",0);
    if (l.isTauIDAvailable("byVVTightIsolationMVA3oldDMwLT")) l.addUserInt("byVVTightIsolationMVA3oldDMwLT",l.tauID("byVVTightIsolationMVA3oldDMwLT"));
    else l.addUserInt("byVVTightIsolationMVA3oldDMwLT",0);
    if (l.isTauIDAvailable("footprintCorrection")) l.addUserFloat("footprintCorrection",l.tauID("footprintCorrection"));
    else l.addUserFloat("footprintCorrection",-999);
    if (l.isTauIDAvailable("neutralIsoPtSumWeight")) l.addUserFloat("neutralIsoPtSumWeight",l.tauID("neutralIsoPtSumWeight"));
    else l.addUserFloat("neutralIsoPtSumWeight",-999);
    if (l.isTauIDAvailable("photonPtSumOutsideSignalCone")) l.addUserFloat("photonPtSumOutsideSignalCone",l.tauID("photonPtSumOutsideSignalCone"));
    else l.addUserFloat("photonPtSumOutsideSignalCone",-999);
    l.addUserFloat("byIsolationMVA3newDMwLTraw",l.tauID("byIsolationMVA3newDMwLTraw"));


    l.addUserFloat("chargedIsoPtSum", chargedIsoPtSum);
    l.addUserFloat("neutralIsoPtSum", neutralIsoPtSum);
    l.addUserFloat("puCorrPtSum", puCorrPtSum);
    l.addUserInt("againstMuonLoose3", againstMuonLoose3);
    l.addUserInt("againstMuonTight3", againstMuonTight3);
    l.addUserInt("againstElectronVLooseMVA6", againstElectronVLooseMVA6);
    l.addUserInt("againstElectronLooseMVA6", againstElectronLooseMVA6);
    l.addUserInt("againstElectronMediumMVA6", againstElectronMediumMVA6);
    l.addUserInt("againstElectronTightMVA6", againstElectronTightMVA6);
    l.addUserInt("againstElectronVTightMVA6", againstElectronVTightMVA6);
    l.addUserInt("numChargedParticlesSignalCone",numChargedParticlesSignalCone);
    l.addUserInt("numNeutralHadronsSignalCone",numNeutralHadronsSignalCone);
    l.addUserInt("numPhotonsSignalCone",numPhotonsSignalCone);
    l.addUserInt("numParticlesSignalCone",numParticlesSignalCone);
    l.addUserInt("numChargedParticlesIsoCone",numChargedParticlesIsoCone);
    l.addUserInt("numNeutralHadronsIsoCone",numNeutralHadronsIsoCone);
    l.addUserInt("numPhotonsIsoCone",numPhotonsIsoCone);
    l.addUserInt("numParticlesIsoCone",numParticlesIsoCone);
    l.addUserFloat("leadChargedParticlePt",leadChargedParticlePt);
    l.addUserFloat("trackRefPt",trackRefPt);

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

