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
  const bool ApplyTESCentralCorr; // shift the central TES value
  //const double NominalTESCorrection;     // value of correction of centrale TES value used in HIG-17-002, same for all DMs
  const double NominalTESCorrection1Pr;    // DM==0  - correction of central TES
  const double NominalTESCorrection1PrPi0; // DM==1  - correction of central TES
  const double NominalTESCorrection3Pr;    // DM==10 - correction of central TES
  const double NominalTESCorrection3PrPi0;    // DM==11 - correction of central TES
  const double NominalTESUncertainty1Pr;      // Up/Down uncertainty for TES 
  const double NominalTESUncertainty1PrPi0;      // Up/Down uncertainty for TES 
  const double NominalTESUncertainty3Pr;      // Up/Down uncertainty for TES 
  const double NominalTESUncertainty3PrPi0;      // Up/Down uncertainty for TES 

  vector<string> tauIntDiscrims_; // tau discrims to be added as userInt
  vector<string> tauFloatDiscrims_; // tau discrims to be added as userFloats
};


TauFiller::TauFiller(const edm::ParameterSet& iConfig) :
  theCandidateTag(consumes<pat::TauRefVector>(iConfig.getParameter<InputTag>("src"))),
  theGenTag(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genCollection"))),
  theVtxTag(consumes<vector<Vertex>>(iConfig.getParameter<edm::InputTag>("vtxCollection"))),
  theDiscriminatorTag(iConfig.getParameter<std::string>("discriminator")),
  cut(iConfig.getParameter<std::string>("cut")),
  flags(iConfig.getParameter<ParameterSet>("flags")), 
  ApplyTESCentralCorr(iConfig.getParameter<bool>("ApplyTESCentralCorr")),
  //NominalTESCorrection(iConfig.getParameter<double>("NominalTESCorrection")), // used for HIG-17-002, same for all values
  NominalTESCorrection1Pr(iConfig.getParameter<double>("NominalTESCorrection1Pr")),
  NominalTESCorrection1PrPi0(iConfig.getParameter<double>("NominalTESCorrection1PrPi0")),
  NominalTESCorrection3Pr(iConfig.getParameter<double>("NominalTESCorrection3Pr")),
  NominalTESCorrection3PrPi0(iConfig.getParameter<double>("NominalTESCorrection3PrPi0")),
  NominalTESUncertainty1Pr(iConfig.getParameter<double>("NominalTESUncertainty1Pr")),
  NominalTESUncertainty1PrPi0(iConfig.getParameter<double>("NominalTESUncertainty1PrPi0")),
  NominalTESUncertainty3Pr(iConfig.getParameter<double>("NominalTESUncertainty3Pr")),
  NominalTESUncertainty3PrPi0(iConfig.getParameter<double>("NominalTESUncertainty3PrPi0"))
{
  produces<pat::TauCollection>();

  tauIntDiscrims_ = 
  {
    "decayModeFinding", // it is decayModeFindingOldDMs
    "decayModeFindingNewDMs",
    
    "byLooseCombinedIsolationDeltaBetaCorr3Hits",
    "byMediumCombinedIsolationDeltaBetaCorr3Hits",
    "byTightCombinedIsolationDeltaBetaCorr3Hits",
    
    "byVLooseIsolationMVArun2v1DBoldDMwLT",
    "byLooseIsolationMVArun2v1DBoldDMwLT",
    "byMediumIsolationMVArun2v1DBoldDMwLT",
    "byTightIsolationMVArun2v1DBoldDMwLT",
    "byVTightIsolationMVArun2v1DBoldDMwLT",

    "byVLooseIsolationMVArun2v1DBnewDMwLT",    
    "byLooseIsolationMVArun2v1DBnewDMwLT",
    "byMediumIsolationMVArun2v1DBnewDMwLT",
    "byTightIsolationMVArun2v1DBnewDMwLT",
    "byVTightIsolationMVArun2v1DBnewDMwLT",

    "byLooseIsolationMVArun2v1DBdR03oldDMwLT",
    "byMediumIsolationMVArun2v1DBdR03oldDMwLT",
    "byTightIsolationMVArun2v1DBdR03oldDMwLT",
    "byVTightIsolationMVArun2v1DBdR03oldDMwLT",
    
    "againstMuonLoose3",
    "againstMuonTight3",

    "againstElectronVLooseMVA6",
    "againstElectronLooseMVA6",
    "againstElectronMediumMVA6",
    "againstElectronTightMVA6",
    "againstElectronVTightMVA6",

    "byVVLooseIsolationMVArun2017v2DBoldDMwLT2017", //FRA syncApr2018
    "byVLooseIsolationMVArun2017v2DBoldDMwLT2017", //FRA syncApr2018
    "byLooseIsolationMVArun2017v2DBoldDMwLT2017",  //FRA syncApr2018
    "byMediumIsolationMVArun2017v2DBoldDMwLT2017", //FRA syncApr2018
    "byTightIsolationMVArun2017v2DBoldDMwLT2017",  //FRA syncApr2018
    "byVTightIsolationMVArun2017v2DBoldDMwLT2017", //FRA syncApr2018
    "byVVTightIsolationMVArun2017v2DBoldDMwLT2017", //FRA syncApr2018
    
    "byVLooseIsolationMVArun2017v1DBoldDMwLT2017", //FRA syncApr2018
    "byLooseIsolationMVArun2017v1DBoldDMwLT2017",  //FRA syncApr2018
    "byMediumIsolationMVArun2017v1DBoldDMwLT2017", //FRA syncApr2018
    "byTightIsolationMVArun2017v1DBoldDMwLT2017",  //FRA syncApr2018
    "byVTightIsolationMVArun2017v1DBoldDMwLT2017", //FRA syncApr2018
    
    "byVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017", //FRA syncApr2018
    "byLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017",  //FRA syncApr2018
    "byMediumIsolationMVArun2017v2DBoldDMdR0p3wLT2017", //FRA syncApr2018
    "byTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017",  //FRA syncApr2018
    "byVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017", //FRA syncApr2018
  };

  tauFloatDiscrims_ =
  {
    "byCombinedIsolationDeltaBetaCorrRaw3Hits",
    "byIsolationMVArun2v1DBoldDMwLTraw",
    "byPhotonPtSumOutsideSignalCone",
    "footprintCorrection",
    "neutralIsoPtSumWeight",
    "photonPtSumOutsideSignalCone",
    "chargedIsoPtSum",
    "neutralIsoPtSum",
    "puCorrPtSum",
    "byIsolationMVArun2017v1DBoldDMwLTraw2017",      //FRA syncApr2018
    "byIsolationMVArun2017v2DBoldDMwLTraw2017",      //FRA syncApr2018
    "byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017", //FRA syncApr2018
  };


}

using LorentzVectorE = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<double>>;

void
TauFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{  

//read one PFTauDiscriminator (set discriminatorSrc_ in to an edm::InputTag before)

  // Get leptons and discriminators
  edm::Handle<pat::TauRefVector> tauHandle;
  iEvent.getByToken(theCandidateTag, tauHandle);
    
  edm::Handle<vector<Vertex> >  vertexs;
  iEvent.getByToken(theVtxTag, vertexs);

  // Output collection
  //auto_ptr<pat::TauCollection> result( new pat::TauCollection() );
  std::unique_ptr<pat::TauCollection> result( new pat::TauCollection() );

  for (unsigned int itau = 0; itau < tauHandle->size(); ++itau){
    //---Clone the pat::Tau
    pat::Tau l(*((*tauHandle)[itau].get()));
    
    // Nominal TES Correction
    double Shift1Pr    = 1. + NominalTESCorrection1Pr/100.;
    double Shift1PrPi0 = 1. + NominalTESCorrection1PrPi0/100.;
    double Shift3Pr    = 1. + NominalTESCorrection3Pr/100.;
    double Shift3PrPi0  = 1. + NominalTESCorrection3PrPi0/100.;
    double shiftP = 1.;
    double shiftMass = 1.;
    
    //if ( l.genJet() && deltaR(l.p4(), l.genJet()->p4()) < 0.5 && l.genJet()->pt() > 8. && ApplyTESCentralCorr) // 2016 data
    if ( l.genJet() && deltaR(l.p4(), l.genJet()->p4()) < 0.2 && l.genJet()->pt() > 15. && ApplyTESCentralCorr)   // 2017 data
    {

      //cout << "---- gen get pt: " << l.genJet()->pt() << endl;
      if (l.decayMode()==0)       // 1prong
      {
        shiftP    = Shift1Pr;
        shiftMass = 1.;
      }
      else if (l.decayMode()==1)  // 1prong+pi0
      {
        shiftP    = Shift1PrPi0;
        shiftMass = Shift1PrPi0;
      }
      else if (l.decayMode()==10) // 3prong
      {
        shiftP    = Shift3Pr;
        shiftMass = Shift3Pr;
      }
      else if (l.decayMode()==11) // 3prong+pi0
      {
        shiftP    = Shift3PrPi0;
        shiftMass = Shift3PrPi0;
      }
      else  // these are not real taus and will be rejected --> we don't care about the shift and just put 1
      {
        shiftP    = 1.;
        shiftMass = 1.;
      }
      //if(l.decayMode()>=1 && l.decayMode()<=10){
      //  shiftP = Shift;
      //  shiftMass = Shift;
      //}
      //else if(l.decayMode()==0){
      //  shiftP = Shift;
      //  shiftMass = 1.;
      //}
    }

    double pxS_Nominal = l.px()*shiftP;
    double pyS_Nominal = l.py()*shiftP;
    double pzS_Nominal = l.pz()*shiftP;
    double massS_Nominal = l.mass()*shiftMass;
    double enS_Nominal = TMath::Sqrt(pxS_Nominal*pxS_Nominal + pyS_Nominal*pyS_Nominal + pzS_Nominal*pzS_Nominal + massS_Nominal*massS_Nominal);
    math::XYZTLorentzVectorD p4S_Nominal( pxS_Nominal, pyS_Nominal, pzS_Nominal, enS_Nominal );

    //l.setP4( p4S_Nominal ); // was here for 2016 data, now it's moved almost at the end of taufiller

    //Up and Down variations: NominalTESUncertainty from python cfg
    //const float udShift[2] = {1.03, 0.97}; // 0: UP, 1: DOWN
    //const double udShift[2] = {1. + (NominalTESUncertainty/100.), 1. - (NominalTESUncertainty/100.)}; // 0: UP, 1: DOWN // unused

    float udshiftP[2] = {1., 1.};
    float udshiftMass[2] = {1., 1.};
    bool isTESShifted = false;
    //if ( l.genJet() && deltaR(l.p4(), l.genJet()->p4()) < 0.5 && l.genJet()->pt() > 8.) { // 2016 data
    if ( l.genJet() && deltaR(l.p4(), l.genJet()->p4()) < 0.2 && l.genJet()->pt() > 15.) { // 2017 data

      isTESShifted = true;

      if (l.decayMode()==0)       // 1prong
      {
        udshiftP[0]    =  Shift1Pr + (NominalTESUncertainty1Pr/100.); //udShift[0]; // up
        udshiftP[1]    =  Shift1Pr - (NominalTESUncertainty1Pr/100.); //udShift[1]; // down
        udshiftMass[0] = udshiftMass[1] = 1.; // no mass shift for pi0
      }
      else if (l.decayMode()==1)  // 1prong+pi0
      {
        udshiftP[0]    = Shift1PrPi0 + (NominalTESUncertainty1PrPi0/100.); //udShift[0]; // up
        udshiftP[1]    = Shift1PrPi0 - (NominalTESUncertainty1PrPi0/100.); //udShift[1]; // down
        udshiftMass[0] = Shift1PrPi0 + (NominalTESUncertainty1PrPi0/100.); //udShift[0]; // up
        udshiftMass[1] = Shift1PrPi0 - (NominalTESUncertainty1PrPi0/100.); //udShift[1]; // down
      }
      else if (l.decayMode()==10) // 3prong
      {
        udshiftP[0]    = Shift3Pr + (NominalTESUncertainty3Pr/100.); //udShift[0]; // up
        udshiftP[1]    = Shift3Pr - (NominalTESUncertainty3Pr/100.); //udShift[1]; // down
        udshiftMass[0] = Shift3Pr + (NominalTESUncertainty3Pr/100.); //udShift[0]; // up
        udshiftMass[1] = Shift3Pr - (NominalTESUncertainty3Pr/100.); //udShift[1]; // down
      }
      else if (l.decayMode()==11) // 3prong
      {
        udshiftP[0]    = Shift3PrPi0 + (NominalTESUncertainty3PrPi0/100.); //udShift[0]; // up
        udshiftP[1]    = Shift3PrPi0 - (NominalTESUncertainty3PrPi0/100.); //udShift[1]; // down
        udshiftMass[0] = Shift3PrPi0 + (NominalTESUncertainty3PrPi0/100.); //udShift[0]; // up
        udshiftMass[1] = Shift3PrPi0 - (NominalTESUncertainty3PrPi0/100.); //udShift[1]; // down
      }
      else  // these are not real taus and will be rejected --> we don't care about the shift and just put 1
      {
        isTESShifted = false;
        udshiftP[0]    = udshiftP[1]    = 1.;
        udshiftMass[0] = udshiftMass[1] = 1.;
      }
      //if(l.decayMode()>=1 && l.decayMode()<=10){
      //  udshiftP[0] = udShift[0]; // up
      //  udshiftP[1] = udShift[1]; // down
      //  udshiftMass[0] = udShift[0]; // up
      //  udshiftMass[1] = udShift[1]; // down
      //}
      //else if(l.decayMode()==0){
      //  udshiftP[0] = udShift[0]; // up
      //  udshiftP[1] = udShift[1]; // down
      //  udshiftMass[0] = udshiftMass[1] = 1.; // no mass shift for pi0
      //}
      //else isTESShifted = false;
    }

    if (isTESShifted)
    {
      // up shift
      double pxS = l.px()*udshiftP[0];
      double pyS = l.py()*udshiftP[0];
      double pzS = l.pz()*udshiftP[0];
      double massS = l.mass()*udshiftMass[0];
      double enS = TMath::Sqrt(pxS*pxS + pyS*pyS + pzS*pzS + massS*massS);
      math::XYZTLorentzVectorD p4SUP( pxS, pyS, pzS, enS );
      // set userfloats
      l.addUserFloat("px_TauUp",p4SUP.px());
      l.addUserFloat("py_TauUp",p4SUP.py());
      l.addUserFloat("pz_TauUp",p4SUP.pz());
      l.addUserFloat("e_TauUp",p4SUP.energy());
      l.addUserFloat("m_TauUp",p4SUP.mass());

      // down shift
      pxS = l.px()*udshiftP[1];
      pyS = l.py()*udshiftP[1];
      pzS = l.pz()*udshiftP[1];
      massS = l.mass()*udshiftMass[1];
      enS = TMath::Sqrt(pxS*pxS + pyS*pyS + pzS*pzS + massS*massS);
      math::XYZTLorentzVectorD p4SDOWN( pxS, pyS, pzS, enS );
      // set userfloats
      l.addUserFloat("px_TauDown",p4SDOWN.px());
      l.addUserFloat("py_TauDown",p4SDOWN.py());
      l.addUserFloat("pz_TauDown",p4SDOWN.pz());
      l.addUserFloat("e_TauDown",p4SDOWN.energy());
      l.addUserFloat("m_TauDown",p4SDOWN.mass());
    }
    l.addUserInt("TauUpExists", isTESShifted ? 1 : 0);
    l.addUserInt("TauDownExists", isTESShifted ? 1 : 0);

    //--- PF ISO
    float PFChargedHadIso   = l.chargedHadronIso();
    float PFNeutralHadIso   = l.neutralHadronIso();
    float PFPhotonIso       = l.photonIso();

    float combRelIsoPF = LeptonIsoHelper::combRelIsoPF(l);

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

    // fill all userfloats
    for (unsigned int iuf = 0; iuf < tauFloatDiscrims_.size(); iuf++)
    {
      string ID = tauFloatDiscrims_.at(iuf);
      l.addUserFloat (ID.c_str(), l.isTauIDAvailable(ID.c_str()) ? l.tauID (ID.c_str()) : -999.);
    }

    // fill all userints
    for (unsigned int iui = 0; iui < tauIntDiscrims_.size(); iui++)
    {
      string ID = tauIntDiscrims_.at(iui);
      int ui = -999;
      if (l.isTauIDAvailable(ID.c_str()))
      {
        ui = ( (l.tauID (ID.c_str()) > 0.5) ? 1 : 0);
      }
      l.addUserInt (ID.c_str(), ui);
    }


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
    
    // apply the actual shift of the central value here
    l.setP4( p4S_Nominal );

    //--- Check selection cut. Being done here, flags are not available; but this way we 
    //    avoid wasting time on rejected leptons.
    if (!cut(l)) continue;

    //--- Embed flags (ie flags specified in the "flags" pset)
    for(CutSet<pat::Tau>::const_iterator flag = flags.begin(); flag != flags.end(); ++flag) {
      l.addUserFloat(flag->first,int((*(flag->second))(l)));
    }
    //cout << " ---------> DOPO pattau - p4: " << LorentzVectorE(l.p4()) << endl;
    result->push_back(l);
  }
  //iEvent.put(result);
  iEvent.put(std::move(result));
}


#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(TauFiller);

