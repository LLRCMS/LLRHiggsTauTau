/** \file
 *
 *  $Date: 2013/05/13 17:10:20 $
 *  $Revision: 1.6 $
 *  \author N. Amapane
 */

#include <LLRHiggsTauTau/NtupleProducer/interface/LeptonIsoHelper.h>
#include <Muon/MuonAnalysisTools/interface/MuonEffectiveArea.h>
//#include <EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/CustomElectronEffectiveArea.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/DaughterDataHelpers.h>
#include <DataFormats/PatCandidates/interface/PackedCandidate.h>
#include "DataFormats/Math/interface/Point3D.h"

#include <TLorentzVector.h>

#include <iostream>

using namespace std;
using namespace edm;
using namespace pat;
using namespace reco;


namespace
{
  struct ByEta
  {
    bool operator()(const pat::PackedCandidate * c1, const pat::PackedCandidate * c2) const
    {
      return c1->eta() < c2->eta();
    }
    bool operator()(float c1eta, const pat::PackedCandidate *c2) const
    {
      return c1eta < c2->eta();
    }
    bool operator()(const pat::PackedCandidate *c1, float c2eta) const
    {
      return c1->eta() < c2eta;
    }
  };
}



int correctionType = 2; //1 = rho; 2 = dbeta;

InputTag LeptonIsoHelper::getMuRhoTag(int sampleType, int setup) {
  InputTag rhoTag;
  if (sampleType ==2011) {
    //rhoTag = InputTag("kt6PFJetsForIso","rho");//RH
    rhoTag = InputTag("fixedGridRhoFastjetAll","");
  } else if (sampleType ==2012) { 
    //rhoTag = InputTag("kt6PFJetsCentralNeutral","rho");//RH
    rhoTag = InputTag("fixedGridRhoFastjetAll","");
  } else {
    cout << "LeptonIsoHelper: Incorrect setup: " << sampleType << " " << setup << endl;
    abort();
  }
  return rhoTag;
}

InputTag LeptonIsoHelper::getEleRhoTag(int sampleType, int setup) {
  InputTag rhoTag;
  if (sampleType ==2011) {
    //rhoTag = InputTag("kt6PFJetsForIso","rho");
    rhoTag = InputTag("fixedGridRhoFastjetAll","");
  } else if (sampleType ==2012) {
    //rhoTag = InputTag("kt6PFJets","rho","RECO");
    rhoTag = InputTag("fixedGridRhoFastjetAll",""); // or "fixedGridRhoFastjetCentralNeutral"? 
  } else {
    cout << "LeptonIsoHelper: Incorect setup: " << sampleType << endl;
    abort();
  }
  return rhoTag;
}


float LeptonIsoHelper::combRelIsoPF(int sampleType, int setup, double rho, const pat::Muon& l, bool dr03, float fsr) {
  float PFChargedHadIso   = l.chargedHadronIso();
  float PFNeutralHadIso   = l.neutralHadronIso();
  float PFPhotonIso       = l.photonIso();
  //float PFPUChargedHadIso = l.puChargedHadronIso();
    
  MuonEffectiveArea::MuonEffectiveAreaTarget EAsetup;
  if (sampleType==2011) {
    EAsetup = MuonEffectiveArea::kMuEAData2011;
  } else if (sampleType ==2012) { 
    EAsetup = MuonEffectiveArea::kMuEAData2012;
  } else abort();

  if (correctionType==1) {
    float EA = MuonEffectiveArea::GetMuonEffectiveArea(MuonEffectiveArea::kMuGammaAndNeutralHadronIso04, 
						       l.eta(), EAsetup);
    return  (PFChargedHadIso + max(0., PFNeutralHadIso + PFPhotonIso - fsr - rho * EA))/l.pt();

  } else if (correctionType==2) {
    //return  (PFChargedHadIso + max(0., PFNeutralHadIso + PFPhotonIso - fsr - 0.5*PFPUChargedHadIso))/l.pt();
    if (dr03)
    {
      return (l.pfIsolationR03().sumChargedHadronPt + max(
             l.pfIsolationR03().sumNeutralHadronEt +
             l.pfIsolationR03().sumPhotonEt - 
             0.5 * l.pfIsolationR03().sumPUPt, 0.0)) / l.pt();
    }
    else
    {
      return (l.pfIsolationR04().sumChargedHadronPt + max(
             l.pfIsolationR04().sumNeutralHadronEt +
             l.pfIsolationR04().sumPhotonEt - 
             0.5 * l.pfIsolationR04().sumPUPt, 0.0)) / l.pt();

    }
  }
  return 0;
}

float LeptonIsoHelper::combRelIsoPF(const pat::Tau& l) {

  float PFChargedHadIso   = l.tauID ("chargedIsoPtSum");
  float PFNeutralHadIso   = l.tauID ("neutralIsoPtSum");
  float PFPhotonIso       = 0;//l.photonIso();
  float PFPUChargedHadIso = l.tauID ("puCorrPtSum");

  return  (PFChargedHadIso + max(0., PFNeutralHadIso + PFPhotonIso - 0.5*PFPUChargedHadIso))/l.pt();

//    return (l.pfIsolationVariables().sumChargedHadronPt + max(
//           l.pfIsolationVariables().sumNeutralHadronEt +
//           l.pfIsolationVariables().sumPhotonEt -
//           0.5 * l.pfIsolationVariables().sumPUPt, 0.0)) / l.pt();

}

float LeptonIsoHelper::combRelIsoPF(int sampleType, int setup, double rho, const pat::Electron& l, float fsr) {
  float PFChargedHadIso   = l.chargedHadronIso();
  float PFNeutralHadIso   = l.neutralHadronIso();
  float PFPhotonIso       = l.photonIso();
  if(correctionType==1){
  ElectronEffectiveArea::ElectronEffectiveAreaTarget EAsetup;
  if (sampleType ==2011) {
    EAsetup = ElectronEffectiveArea::kEleEAData2011;
  } else if (sampleType ==2012) { 
    EAsetup = ElectronEffectiveArea::kEleEAData2012; // Legacy
    // EAsetup = ElectronEffectiveArea::kEleEASpring14MC_PU20bx25 // retuned by Simon
  } else {
    abort();
  }

  float EA = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso04,
							     l.superCluster()->eta(), EAsetup);
  return  (PFChargedHadIso + max(0., PFNeutralHadIso + PFPhotonIso - fsr - rho * EA))/l.pt();
  }else{
    return (l.pfIsolationVariables().sumChargedHadronPt + max(
           l.pfIsolationVariables().sumNeutralHadronEt +
           l.pfIsolationVariables().sumPhotonEt - 
           0.5 * l.pfIsolationVariables().sumPUPt, 0.0)) / l.pt();
  }

}


float LeptonIsoHelper::combRelIsoPF(int sampleType, int setup, double rho, const Candidate* lep, float fsr) {
  // should check if lep->hasMasterClone()?  
  if (lep->isMuon()) {
    const pat::Muon* mu = dynamic_cast<const pat::Muon*>(lep->masterClone().get());
    return combRelIsoPF(sampleType, setup, rho, *mu, false, fsr); //NB: defaults to 2016 0.4 mu cone
  } else if (lep->isElectron()) {
    const pat::Electron* ele = dynamic_cast<const pat::Electron*>(lep->masterClone().get());
    return combRelIsoPF(sampleType, setup, rho, *ele, fsr);    
  }else {
    cout << "ERROR: LeptonIsoHelper: unknown type" << endl;
    abort();
  }
  return 0;
}





int LeptonIsoHelper::jetNDauChargedMVASel(const reco::Candidate* cand, pat::Jet jet){

  int jetNDauCharged = 0;

  for (unsigned int i = 0, n = jet.numberOfSourceCandidatePtrs(); i < n; ++i) {
      
    const pat::PackedCandidate &dau_jet = dynamic_cast<const pat::PackedCandidate &>(*(jet.sourceCandidatePtr(i)));
    
    if (!dau_jet.bestTrack()) //FRA 
    {
      //std::cout << " ------ Skipping.." << std::endl;
      continue;
    }

    float dR = deltaR(jet,dau_jet);
    
    bool isgoodtrk = false;
    const reco::Track trk = dau_jet.pseudoTrack();
    const math::XYZPoint vtx_position = cand->vertex();
    
    if(trk.pt()>1 &&
       trk.hitPattern().numberOfValidHits()>=8 &&
       trk.hitPattern().numberOfValidPixelHits()>=2 &&
       trk.normalizedChi2()<5 &&
       std::fabs(trk.dxy(vtx_position))<0.2 &&
       std::fabs(trk.dz(vtx_position))<17
       ) isgoodtrk = true;
    
    if( dR<=0.4 && dau_jet.charge()!=0 && dau_jet.fromPV()>1 && isgoodtrk)
      jetNDauCharged++;
    
  }
  
  return jetNDauCharged;

}




void LeptonIsoHelper::PFIso_particles(const edm::View<pat::PackedCandidate>* pfCands, std::vector<const pat::PackedCandidate *> & pfCands_charged, std::vector<const pat::PackedCandidate *> & pfCands_neutral){

  pfCands_charged.clear(); pfCands_neutral.clear();

  for(edm::View<pat::PackedCandidate>::const_iterator pfCandi = pfCands->begin(); pfCandi!=pfCands->end(); ++pfCandi){

    const pat::PackedCandidate* p = &(*pfCandi);

    if (p->charge() == 0) {
      pfCands_neutral.push_back(p);
    } 

    else{
      if( fabs(p->pdgId()) == 211 ){
	if (p->fromPV() > 1 && fabs(p->dz()) < 9999. ){
	  pfCands_charged.push_back(p);
	}
      }
    }

  }

  std::sort(pfCands_charged.begin(), pfCands_charged.end(), ByEta());
  std::sort(pfCands_neutral.begin(), pfCands_neutral.end(), ByEta());

  return;

}


float LeptonIsoHelper::isoSumRaw(const reco::Candidate* cand, const std::vector<const pat::PackedCandidate *> pfCands_Iso, float dR, float innerR, float threshold, SelfVetoPolicy::SelfVetoPolicy selfVeto, int pdgId){


  std::vector<const reco::Candidate *> vetos;
   for( unsigned int i=0,n=cand->numberOfSourceCandidatePtrs();i<n;++i )
     {
        if(selfVeto == SelfVetoPolicy::selfVetoNone) break;
        const reco::CandidatePtr &cp = cand->sourceCandidatePtr(i);
        if( cp.isNonnull() && cp.isAvailable() )
	  {
	     vetos.push_back(&*cp);
	     if (selfVeto == SelfVetoPolicy::selfVetoFirst) break;
        }
    }

   typedef std::vector<const pat::PackedCandidate *>::const_iterator IT;
   IT candsbegin = std::lower_bound(pfCands_Iso.begin(), pfCands_Iso.end(), cand->eta() - dR, ByEta());
   IT candsend = std::upper_bound(candsbegin, pfCands_Iso.end(), cand->eta() + dR, ByEta());

   double isosum = 0;
   for( IT icand=candsbegin;icand<candsend;++icand )
     {
        // pdgId
        if( pdgId > 0 && abs((*icand)->pdgId()) != pdgId ) continue;
        // threshold
        if( threshold > 0 && (*icand)->pt() < threshold ) continue;
        // cone
	float current_dR = deltaR((*icand)->eta(),(*icand)->phi(),cand->eta(),cand->phi());
        if( current_dR >= dR || current_dR < innerR ) continue;
        // veto
        if( std::find(vetos.begin(), vetos.end(), *icand) != vetos.end() )
	  {
	     continue;
        }

        // add to sum
        isosum += (*icand)->pt();
    }
   return isosum;


}





float LeptonIsoHelper::PfIsoCharged(const reco::Candidate* cand, const std::vector<const pat::PackedCandidate *> pfCands_charged, float miniIsoR){

  float result = -1.;

  if(cand->isElectron()){

    float innerR_Ch = .0;

    if( userdatahelpers::getUserInt(cand,"isEB") ) { innerR_Ch = 0.0; }
    else { innerR_Ch = 0.015; }

    result = isoSumRaw(cand,pfCands_charged,miniIsoR,innerR_Ch,0.0,SelfVetoPolicy::selfVetoNone);

  }

  else if(cand->isMuon()){

    result = isoSumRaw(cand,pfCands_charged,miniIsoR,0.0001,0.0,SelfVetoPolicy::selfVetoAll);;

  }
  
  return result;

}



float LeptonIsoHelper::PfIsoNeutral(const reco::Candidate* cand, const std::vector<const pat::PackedCandidate *> pfCands_neutral, float miniIsoR)
{

  float result = -1.;

  if(cand->isElectron()){

    float innerR_N = .0;
    if( userdatahelpers::getUserInt(cand,"isEB") ) { innerR_N = 0.0; }
    else { innerR_N = 0.08; }
    
    float result1 = isoSumRaw(cand,pfCands_neutral,miniIsoR,innerR_N,0.0,SelfVetoPolicy::selfVetoNone,22 );
    float result2 = isoSumRaw(cand,pfCands_neutral,miniIsoR,0.0     ,0.0,SelfVetoPolicy::selfVetoNone,130);
    result = result1 + result2;

  }

  else if(cand->isMuon()){

    result = isoSumRaw(cand,pfCands_neutral,miniIsoR,0.01,0.5,SelfVetoPolicy::selfVetoAll);

  }

  return result;
}



std::pair<float,float> LeptonIsoHelper::miniRelIso_ChargedNeutral(const reco::Candidate* cand, const std::vector<const pat::PackedCandidate *> pfCands_charged, const std::vector<const pat::PackedCandidate *> pfCands_neutral, float rho, int year){


  float miniIsoR = 10.0/std::min(std::max(float(cand->pt()),float(50.)),float(200.));
  float EffArea = 0.;
  float eta = cand->eta();

  if(cand->isElectron()){
  
    // Valid for 2016 data
    if (year == 2016)
    {
      if     ( fabs(eta) >= 0     && fabs(eta) < 1.0   ) EffArea = 0.1752;
      else if( fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffArea = 0.1862;
      else if( fabs(eta) >= 1.479 && fabs(eta) < 2.0   ) EffArea = 0.1411;
      else if( fabs(eta) >= 2.0   && fabs(eta) < 2.2   ) EffArea = 0.1534;
      else if( fabs(eta) >= 2.2   && fabs(eta) < 2.3   ) EffArea = 0.1903;
      else if( fabs(eta) >= 2.3   && fabs(eta) < 2.4   ) EffArea = 0.2243;
      else if( fabs(eta) >= 2.4   && fabs(eta) <= 2.5  ) EffArea = 0.2687;
    }

    // Valid for 2017 and 2018 data
    if (year == 2017 || year == 2018)
    {
      if     ( fabs(eta) >= 0     && fabs(eta) < 1.0   ) EffArea = 0.1440;
      else if( fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffArea = 0.1562;
      else if( fabs(eta) >= 1.479 && fabs(eta) < 2.0   ) EffArea = 0.1032;
      else if( fabs(eta) >= 2.0   && fabs(eta) < 2.2   ) EffArea = 0.0859;
      else if( fabs(eta) >= 2.2   && fabs(eta) < 2.3   ) EffArea = 0.1116;
      else if( fabs(eta) >= 2.3   && fabs(eta) < 2.4   ) EffArea = 0.1321;
      else if( fabs(eta) >= 2.4   && fabs(eta) <= 2.5  ) EffArea = 0.1654;
    }

  }

  else if(cand->isMuon()){

    // Valid for 2016 data
    if (year == 2016)
    {
      if     ( fabs(eta) >= 0   && fabs(eta) < 0.8  ) EffArea = 0.0735;
      else if( fabs(eta) >= 0.8 && fabs(eta) < 1.3  ) EffArea = 0.0619;
      else if( fabs(eta) >= 1.3 && fabs(eta) < 2.0  ) EffArea = 0.0465;
      else if( fabs(eta) >= 2.0 && fabs(eta) < 2.2  ) EffArea = 0.0433;
      else if( fabs(eta) >= 2.2 && fabs(eta) <= 2.5 ) EffArea = 0.0577;
    }

    // Valid for 2017 and 2018 data
    if (year == 2017 || year == 2018)
    {
      if     ( fabs(eta) >= 0   && fabs(eta) < 0.8  ) EffArea = 0.0566;
      else if( fabs(eta) >= 0.8 && fabs(eta) < 1.3  ) EffArea = 0.0562;
      else if( fabs(eta) >= 1.3 && fabs(eta) < 2.0  ) EffArea = 0.0363;
      else if( fabs(eta) >= 2.0 && fabs(eta) < 2.2  ) EffArea = 0.0119;
      else if( fabs(eta) >= 2.2 && fabs(eta) <= 2.5 ) EffArea = 0.0064;
    }

  }

  float correction = rho*EffArea*(miniIsoR/0.3)*(miniIsoR/0.3);
    
  float pfIsoCharged = PfIsoCharged(cand,pfCands_charged,miniIsoR);
  float pfIsoNeutral = PfIsoNeutral(cand,pfCands_neutral,miniIsoR);
  float pfIsoPUSubtracted = std::max(float(0.0),float(pfIsoNeutral-correction));
  
  float miniRelIsoCharged = pfIsoCharged/cand->pt();
  float miniRelIsoNeutral = pfIsoPUSubtracted/cand->pt();

  pair<float,float> miniRelIso = make_pair( miniRelIsoCharged, miniRelIsoNeutral );
  return miniRelIso;
  

}





float LeptonIsoHelper::jetPtRel(const reco::Candidate& cand, const pat::Jet& jet, string JECname)
{

  float PtRel = 0;

  if(jet.numberOfDaughters()>1){
    
    pat::Jet myCorJet;    
    myCorJet.setP4(jet.correctedJet("L1FastJet","none",JECname).p4());
    
    float          SF          = jet.p4().E() / myCorJet.p4().E();
    
    auto lepAwareJetp4 = ( myCorJet.p4() - cand.p4() ) * SF + cand.p4();
    
    TLorentzVector candV = TLorentzVector(cand.px(),cand.py(),cand.pz(),cand.p4().E());
    TLorentzVector jetV = TLorentzVector(lepAwareJetp4.px(),lepAwareJetp4.py(),lepAwareJetp4.pz(),lepAwareJetp4.E());
    
    PtRel = candV.Perp( (jetV - candV).Vect() );
    
  }

  return (PtRel > 0) ? PtRel : 0.0;
}




float LeptonIsoHelper::jetPtRatio(const reco::Candidate& cand, const pat::Jet& jet, string JECname)
{

  float PtRatio = 1.;

  if(jet.numberOfDaughters()>1){

    pat::Jet myCorJet;
    myCorJet.setP4(jet.correctedJet("L1FastJet","none",JECname).p4());

    float          SF          = jet.p4().E() / myCorJet.p4().E();    
    auto lepAwareJetp4 = ( myCorJet.p4() - cand.p4() ) * SF + cand.p4();    

    PtRatio = cand.pt() / lepAwareJetp4.pt();
    
  }

  return (PtRatio > 0) ? PtRatio : 0.0;
}
