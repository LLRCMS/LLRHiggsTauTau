/*
  Adapted from michelif/flashgg implementation:
  https://github.com/michelif/flashgg/blob/bregProducer_90X/Taggers/plugins/flashggbRegressionProducer.cc

 TODO:
 - Find out which is the correct weight file
 - Find out if/how to consider only b-jets (or b-reg applied to ALL jets?)
*/

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "RecoJets/JetProducers/interface/PileupJetIdAlgo.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"

#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include <iostream>
#include <string>
#include <vector>


#define debug 0

using namespace std;
using namespace edm;

class bRegressionProducer : public EDProducer
{

public:
    bRegressionProducer( const ParameterSet & );
    ~bRegressionProducer(){};
    void InitJet();
    void SetNNVectorVar();
    std::vector<float> EvaluateNN();
private:
    void produce( Event &, const EventSetup & ) override;

    edm::InputTag inputTagJets_;
    EDGetTokenT<View<pat::Jet> > jetToken_;
    edm::EDGetTokenT<double> rhoToken_;      
 
    FileInPath bRegressionWeightfile_;

    edm::EDGetTokenT<std::vector<reco::Vertex>> srcVtx_;
    edm::EDGetTokenT<edm::View<reco::VertexCompositePtrCandidate>> srcSV_;  

    double y_mean_;
    double y_std_;

    tensorflow::Session* session;
    std::vector<float> NNvectorVar_; 


    //mva variables
    float Jet_pt ;
    float Jet_eta ;
    float rho ;
    float Jet_mt ;
    float Jet_leadTrackPt ;
    float Jet_leptonPtRel ;
    float Jet_leptonDeltaR ;
    float Jet_neHEF ;
    float Jet_neEmEF ;
    float Jet_vtxPt ;
    float Jet_vtxMass ;
    float Jet_vtx3dL ;
    float Jet_vtxNtrk ;
    float Jet_vtx3deL ;
    float Jet_energyRing_dR0_em_Jet_e ;
    float Jet_energyRing_dR1_em_Jet_e ;
    float Jet_energyRing_dR2_em_Jet_e ;
    float Jet_energyRing_dR3_em_Jet_e ;
    float Jet_energyRing_dR4_em_Jet_e ;
    float Jet_energyRing_dR0_neut_Jet_e ;
    float Jet_energyRing_dR1_neut_Jet_e ;
    float Jet_energyRing_dR2_neut_Jet_e ;
    float Jet_energyRing_dR3_neut_Jet_e ;
    float Jet_energyRing_dR4_neut_Jet_e ;
    float Jet_energyRing_dR0_ch_Jet_e ;
    float Jet_energyRing_dR1_ch_Jet_e ;
    float Jet_energyRing_dR2_ch_Jet_e ;
    float Jet_energyRing_dR3_ch_Jet_e ;
    float Jet_energyRing_dR4_ch_Jet_e ;
    float Jet_energyRing_dR0_mu_Jet_e ;
    float Jet_energyRing_dR1_mu_Jet_e ;
    float Jet_energyRing_dR2_mu_Jet_e ;
    float Jet_energyRing_dR3_mu_Jet_e ;
    float Jet_energyRing_dR4_mu_Jet_e ;
    float Jet_numDaughters_pt03 ;
    float Jet_chHEF;//implement from here
    float Jet_chEmEF;
    float Jet_leptonPtRelInv;
    int isEle;
    int isMu;
    int isOther;
    float Jet_mass;
    float Jet_withPtd;

        
};


bRegressionProducer::bRegressionProducer( const ParameterSet &iConfig ) :
    inputTagJets_( iConfig.getParameter<edm::InputTag>( "JetTag" )) ,
    rhoToken_( consumes<double>(iConfig.getParameter<edm::InputTag>( "rhoFixedGridCollection" ) )),
    bRegressionWeightfile_(iConfig.getParameter<edm::FileInPath>("bRegressionWeightfile")),
    srcVtx_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("pvsrc"))),
    srcSV_(consumes<edm::View<reco::VertexCompositePtrCandidate>>(iConfig.getParameter<edm::InputTag>("svsrc"))),
    y_mean_(iConfig.getParameter<double>("y_mean")),
    y_std_(iConfig.getParameter<double>("y_std"))

{
    jetToken_= consumes<View<pat::Jet> >(inputTagJets_);


    tensorflow::GraphDef* graphDef= tensorflow::loadGraphDef(bRegressionWeightfile_.fullPath()); 

    session = tensorflow::createSession(graphDef);

    //for variables for breg check this PR https://github.com/cms-analysis/flashgg/pull/968
    InitJet();

    produces<vector<pat::Jet> > ();
}



void bRegressionProducer::produce( Event &evt, const EventSetup & )
{
    // input jets
    Handle<View<pat::Jet> > jets;
    evt.getByToken( jetToken_, jets );//just to try get the first one
    unique_ptr<vector<pat::Jet> > jetColl( new vector<pat::Jet> );
  
    const auto& vtxProd = evt.get(srcVtx_);
    const auto& svProd = evt.get(srcSV_);


    for( unsigned int i = 0 ; i < jets->size() ; i++ ) {
        // Reset all input features
        InitJet();
            
        Ptr<pat::Jet> pjet = jets->ptrAt( i );
        pat::Jet fjet = pat::Jet( *pjet );

        if (fjet.pt()<15. || fabs(fjet.eta())>2.5) continue;

        std::array<float, 5> cone_boundaries{{ 0.05, 0.1, 0.2, 0.3, 0.4 }}; // hardcoded boundaries: should be made configurable
        size_t ncone_boundaries = cone_boundaries.size();
        std::vector<float> chEnergies(ncone_boundaries+1,0.);
        std::vector<float> emEnergies(ncone_boundaries+1,0.); 
        std::vector<float> neEnergies(ncone_boundaries+1,0.); 
        std::vector<float> muEnergies(ncone_boundaries+1,0.); 
        
        float leadTrackPt_ = 0, softLepPt = 0, softLepDr = 0;

        float softLepPtRel = 0.;
        float softLepPtRelInv=0.;
        int softLepPDGId=0.;
        int numDaug03 = 0;

        for ( unsigned k = 0; k < fjet.numberOfSourceCandidatePtrs(); ++k ) {
            reco::CandidatePtr pfJetConstituent = fjet.sourceCandidatePtr(k);
                    
            const reco::Candidate* kcand = pfJetConstituent.get();
            const pat::PackedCandidate* lPack = dynamic_cast<const pat::PackedCandidate *>( kcand );
            if ( !lPack ) throw cms::Exception( "NoPackedConstituent" ) << " For jet " << i << " failed to get constituent " << k << std::endl;
            int PDGID = abs(lPack->pdgId()); 
            float candPt = kcand->pt();
            float candDr   = reco::deltaR(*kcand,fjet);

            if( candPt > 0.3 ) { ++numDaug03; }
            if(lPack->charge() != 0 && candPt > leadTrackPt_) leadTrackPt_ = candPt;

            if(PDGID == 11 || PDGID == 13) {
                if(candPt > softLepPt){
                    
                    softLepPt = candPt;
                    softLepDr = candDr;
                    
                    // NanoAOD version here: https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/plugins/BJetEnergyRegressionVarProducer.cc#L248
                    softLepPtRel = ( pjet->px()*lPack->px() + 
                                     pjet->py()*lPack->py() + 
                                     pjet->pz()*lPack->pz() ) / pjet->p();

                    softLepPtRel = sqrt( lPack->p()*lPack->p() - 
                                         softLepPtRel*softLepPtRel );

                    softLepPtRelInv = ( pjet->px()*lPack->px() + 
                                        pjet->py()*lPack->py() +
                                        pjet->pz()*lPack->pz() ) / lPack->p();

                    softLepPtRelInv = sqrt( pjet->p()*pjet->p() - 
                                            softLepPtRelInv*softLepPtRelInv );

                    softLepPDGId = lPack->pdgId();
                }
            }

            size_t icone = std::lower_bound(&cone_boundaries[0],
                                            &cone_boundaries[ncone_boundaries],
                                            candDr) - &cone_boundaries[0];
            float candEnergy = kcand->energy();

            if( PDGID == 22 || PDGID == 11 ) {
                emEnergies[icone] += candEnergy;
            } else if ( PDGID == 13 ) { 
                muEnergies[icone] += candEnergy;
            } else if ( lPack-> charge() != 0 ) {
                chEnergies[icone] += candEnergy;
            } else {
                neEnergies[icone] += candEnergy;
            }
        }
        
        if (abs(softLepPDGId)==13){
            isMu=1; 
        }else if (abs(softLepPDGId)==11){
            isEle=1;
        }else{
            isOther=1;
        }


        float maxFoundSignificance = 0;
        VertexDistance3D vdist;
        int NSV = 0;
        const auto& pv = vtxProd.at(0);
        for (const auto& sv : svProd) {
            NSV++;

            GlobalVector flightDir(sv.vertex().x() - pv.x(), sv.vertex().y() - pv.y(), sv.vertex().z() - pv.z());
            GlobalVector jetDir(pjet->px(), pjet->py(), pjet->pz());
            if (reco::deltaR2(flightDir, jetDir) < 0.09) {
                Measurement1D dl = vdist.distance(
                                                  pv, VertexState(RecoVertex::convertPos(sv.position()), RecoVertex::convertError(sv.error())));
                if (dl.significance() > maxFoundSignificance) {
                    maxFoundSignificance = dl.significance();
                    Jet_vtxPt = sv.pt();
                    Jet_vtxMass = sv.p4().M();
                    Jet_vtx3dL = dl.value();
                    Jet_vtx3deL = dl.error();
                    Jet_vtxNtrk = sv.numberOfSourceCandidatePtrs();
                }
            }
        }
        //variables needed for regression
        Jet_pt = fjet.pt();
        Jet_eta = fjet.eta() ;
        Jet_leadTrackPt = leadTrackPt_; 
        edm::Handle<double> rhoHandle;
        evt.getByToken( rhoToken_, rhoHandle );
        const double rhoFixedGrd = *( rhoHandle.product() );
        rho = rhoFixedGrd;
        Jet_mt = sqrt(fjet.energy()*fjet.energy()-fjet.pz()*fjet.pz());

        //this max probably not needed, it's just heppy
        Jet_leptonPtRel = std::max(float(0.),softLepPtRel);
        Jet_leptonDeltaR = std::max(float(0.),softLepDr);
        Jet_neHEF = fjet.neutralHadronEnergyFraction();
        Jet_neEmEF = fjet.neutralEmEnergyFraction();
        Jet_chHEF = fjet.chargedHadronEnergyFraction(); 
        Jet_chEmEF = fjet.chargedEmEnergyFraction();

        Jet_leptonPtRelInv = std::max(float(0.),softLepPtRelInv);
        Jet_mass=fjet.mass();

        float sumWeight=0;
        float sumPt=0;
        for(const auto & d : pjet->daughterPtrVector()){
            sumWeight+=(d->pt())*(d->pt());
            sumPt+=d->pt();
        }
        Jet_withPtd = (sumWeight > 0 ? sqrt(sumWeight)/sumPt : 0);

    
        

        Jet_energyRing_dR0_em_Jet_e   = emEnergies[0]/fjet.energy();
        Jet_energyRing_dR1_em_Jet_e   = emEnergies[1]/fjet.energy();
        Jet_energyRing_dR2_em_Jet_e   = emEnergies[2]/fjet.energy();
        Jet_energyRing_dR3_em_Jet_e   = emEnergies[3]/fjet.energy();
        Jet_energyRing_dR4_em_Jet_e   = emEnergies[4]/fjet.energy();
        Jet_energyRing_dR0_neut_Jet_e = neEnergies[0]/fjet.energy();
        Jet_energyRing_dR1_neut_Jet_e = neEnergies[1]/fjet.energy();
        Jet_energyRing_dR2_neut_Jet_e = neEnergies[2]/fjet.energy();
        Jet_energyRing_dR3_neut_Jet_e = neEnergies[3]/fjet.energy();
        Jet_energyRing_dR4_neut_Jet_e = neEnergies[4]/fjet.energy();
        Jet_energyRing_dR0_ch_Jet_e   = chEnergies[0]/fjet.energy();
        Jet_energyRing_dR1_ch_Jet_e   = chEnergies[1]/fjet.energy();
        Jet_energyRing_dR2_ch_Jet_e   = chEnergies[2]/fjet.energy();
        Jet_energyRing_dR3_ch_Jet_e   = chEnergies[3]/fjet.energy();
        Jet_energyRing_dR4_ch_Jet_e   = chEnergies[4]/fjet.energy();
        Jet_energyRing_dR0_mu_Jet_e   = muEnergies[0]/fjet.energy();
        Jet_energyRing_dR1_mu_Jet_e   = muEnergies[1]/fjet.energy();
        Jet_energyRing_dR2_mu_Jet_e   = muEnergies[2]/fjet.energy();
        Jet_energyRing_dR3_mu_Jet_e   = muEnergies[3]/fjet.energy();
        Jet_energyRing_dR4_mu_Jet_e   = muEnergies[4]/fjet.energy();

        Jet_numDaughters_pt03 = numDaug03;
            
        std::vector<float> bRegNN(3,-999);
            

        if(debug){
            cout<<"Jet_pt :"<<Jet_pt <<endl;
            cout<<"Jet_eta :"<<Jet_eta <<endl;
            cout<<"rho :"<<rho <<endl;
            cout<<"Jet_mt :"<<Jet_mt <<endl;
            cout<<"Jet_leadTrackPt :"<<Jet_leadTrackPt <<endl;
            cout<<"Jet_leptonPtRel :"<<Jet_leptonPtRel <<endl;
            cout<<"Jet_leptonDeltaR :"<<Jet_leptonDeltaR <<endl;
            cout<<"Jet_neHEF :"<<Jet_neHEF <<endl;
            cout<<"Jet_neEmEF :"<<Jet_neEmEF <<endl;
            cout<<"Jet_vtxPt :"<<Jet_vtxPt <<endl;
            cout<<"Jet_vtxMass :"<<Jet_vtxMass <<endl;
            cout<<"Jet_vtx3dL :"<<Jet_vtx3dL <<endl;
            cout<<"Jet_vtxNtrk :"<<Jet_vtxNtrk <<endl;
            cout<<"Jet_vtx3deL :"<<Jet_vtx3deL <<endl;
            cout<<"Jet_energyRing_dR0_em_Jet_e :"<<Jet_energyRing_dR0_em_Jet_e <<endl;
            cout<<"Jet_energyRing_dR1_em_Jet_e :"<<Jet_energyRing_dR1_em_Jet_e <<endl;
            cout<<"Jet_energyRing_dR2_em_Jet_e :"<<Jet_energyRing_dR2_em_Jet_e <<endl;
            cout<<"Jet_energyRing_dR3_em_Jet_e :"<<Jet_energyRing_dR3_em_Jet_e <<endl;
            cout<<"Jet_energyRing_dR4_em_Jet_e :"<<Jet_energyRing_dR4_em_Jet_e <<endl;
            cout<<"Jet_energyRing_dR0_neut_Jet_e :"<<Jet_energyRing_dR0_neut_Jet_e <<endl;
            cout<<"Jet_energyRing_dR1_neut_Jet_e :"<<Jet_energyRing_dR1_neut_Jet_e <<endl;
            cout<<"Jet_energyRing_dR2_neut_Jet_e :"<<Jet_energyRing_dR2_neut_Jet_e <<endl;
            cout<<"Jet_energyRing_dR3_neut_Jet_e :"<<Jet_energyRing_dR3_neut_Jet_e <<endl;
            cout<<"Jet_energyRing_dR4_neut_Jet_e :"<<Jet_energyRing_dR4_neut_Jet_e <<endl;
            cout<<"Jet_energyRing_dR0_ch_Jet_e :"<<Jet_energyRing_dR0_ch_Jet_e <<endl;
            cout<<"Jet_energyRing_dR1_ch_Jet_e :"<<Jet_energyRing_dR1_ch_Jet_e <<endl;
            cout<<"Jet_energyRing_dR2_ch_Jet_e :"<<Jet_energyRing_dR2_ch_Jet_e <<endl;
            cout<<"Jet_energyRing_dR3_ch_Jet_e :"<<Jet_energyRing_dR3_ch_Jet_e <<endl;
            cout<<"Jet_energyRing_dR4_ch_Jet_e :"<<Jet_energyRing_dR4_ch_Jet_e <<endl;
            cout<<"Jet_energyRing_dR0_mu_Jet_e :"<<Jet_energyRing_dR0_mu_Jet_e <<endl;
            cout<<"Jet_energyRing_dR1_mu_Jet_e :"<<Jet_energyRing_dR1_mu_Jet_e <<endl;
            cout<<"Jet_energyRing_dR2_mu_Jet_e :"<<Jet_energyRing_dR2_mu_Jet_e <<endl;
            cout<<"Jet_energyRing_dR3_mu_Jet_e :"<<Jet_energyRing_dR3_mu_Jet_e <<endl;
            cout<<"Jet_energyRing_dR4_mu_Jet_e :"<<Jet_energyRing_dR4_mu_Jet_e <<endl;
            cout<<"Jet_numDaughters_pt03 :"<<Jet_numDaughters_pt03 <<endl;
            cout<<"Jet_chHEF:"<<Jet_chHEF<<endl;
            cout<<"Jet_chEmEF:"<<Jet_chEmEF<<endl;
            cout<<"Jet_leptonPtRelInv:"<<Jet_leptonPtRelInv<<endl;
            cout<<"isEle:"<<isEle<<endl;
            cout<<"isMu:"<<isMu<<endl;
            cout<<"isOther:"<<isOther<<endl;
            cout<<"Jet_mass:"<<Jet_mass<<endl;
            cout<<"Jet_withPtd:"<<Jet_withPtd<<endl;
            cout << endl;


        }

        //        //..... gen jets info     [ DO WE NEED TO DO SOMETHING LIKE THAT?]                                                                              
        //        int cflav = 0; //~correct flavour definition
        //        if ( !evt.isRealData() ) {
        //            int hflav = fjet.hadronFlavour();//4 if c, 5 if b, 0 if light jets
        //            int pflav = fjet.partonFlavour();
        //
        //            if( hflav != 0 ) {
        //                cflav = hflav;
        //            } else { //not a heavy jet                                              
        //                cflav = std::abs(pflav) == 4 || std::abs(pflav) == 5 ? 0 : pflav;
        //            }
        //            std::cout<< "cflav==" << cflav<<std::endl;
        //            //if (cflav != 5) continue;//i want only bjets
        //        }


        SetNNVectorVar();
        bRegNN = EvaluateNN();
        NNvectorVar_.clear();

        fjet.addUserFloat("bRegNNCorr", bRegNN[0]*y_std_+y_mean_);
        fjet.addUserFloat("bRegNNResolution",0.5*(bRegNN[2]-bRegNN[1])*y_std_);

        if(debug){
            cout<<"bRegNNCorr = "       << bRegNN[0]*y_std_+y_mean_ << endl;
            cout<<"bRegNNResolution = " << 0.5*(bRegNN[2]-bRegNN[1])*y_std_ << endl;
            cout << "===" << endl << endl << "===" << endl;
        }

        jetColl->push_back( fjet );

            

    }
    evt.put( std::move( jetColl ) );
}
    
void bRegressionProducer::InitJet(){
    Jet_pt = 0.;
    Jet_eta = 0.;
    rho = 0.;
    Jet_mt = 0.;
    Jet_leadTrackPt = 0.;
    Jet_leptonPtRel = 0.;
    Jet_leptonDeltaR = 0.;
    Jet_neHEF = 0.;
    Jet_neEmEF = 0.;
    Jet_vtxPt = 0.;
    Jet_vtxMass = 0.;
    Jet_vtx3dL = 0.;
    Jet_vtxNtrk = 0.;
    Jet_vtx3deL = 0.;
    Jet_energyRing_dR0_em_Jet_e = 0.;
    Jet_energyRing_dR1_em_Jet_e = 0.;
    Jet_energyRing_dR2_em_Jet_e = 0.;
    Jet_energyRing_dR3_em_Jet_e = 0.;
    Jet_energyRing_dR4_em_Jet_e = 0.;
    Jet_energyRing_dR0_neut_Jet_e = 0.;
    Jet_energyRing_dR1_neut_Jet_e = 0.;
    Jet_energyRing_dR2_neut_Jet_e = 0.;
    Jet_energyRing_dR3_neut_Jet_e = 0.;
    Jet_energyRing_dR4_neut_Jet_e = 0.;
    Jet_energyRing_dR0_ch_Jet_e = 0.;
    Jet_energyRing_dR1_ch_Jet_e = 0.;
    Jet_energyRing_dR2_ch_Jet_e = 0.;
    Jet_energyRing_dR3_ch_Jet_e = 0.;
    Jet_energyRing_dR4_ch_Jet_e = 0.;
    Jet_energyRing_dR0_mu_Jet_e = 0.;
    Jet_energyRing_dR1_mu_Jet_e = 0.;
    Jet_energyRing_dR2_mu_Jet_e = 0.;
    Jet_energyRing_dR3_mu_Jet_e = 0.;
    Jet_energyRing_dR4_mu_Jet_e = 0.;
    Jet_numDaughters_pt03 = 0;

    Jet_chHEF = 0.;//implement from here
    Jet_chEmEF = 0.;
    Jet_leptonPtRelInv = 0.;
    isEle = 0;
    isMu = 0;
    isOther = 0;
    Jet_mass = 0.;
    Jet_withPtd = 0.;

}//end InitJet

void bRegressionProducer::SetNNVectorVar(){

    NNvectorVar_.push_back(Jet_pt) ;//0
    NNvectorVar_.push_back(Jet_eta) ;
    NNvectorVar_.push_back(rho) ;
    NNvectorVar_.push_back(Jet_mt) ;
    NNvectorVar_.push_back(Jet_leadTrackPt) ;
    NNvectorVar_.push_back(Jet_leptonPtRel) ;//5
    NNvectorVar_.push_back(Jet_leptonDeltaR) ;
    NNvectorVar_.push_back(Jet_neHEF) ;
    NNvectorVar_.push_back(Jet_neEmEF) ;
    NNvectorVar_.push_back(Jet_vtxPt) ;
    NNvectorVar_.push_back(Jet_vtxMass) ;//10
    NNvectorVar_.push_back(Jet_vtx3dL) ;
    NNvectorVar_.push_back(Jet_vtxNtrk) ;
    NNvectorVar_.push_back(Jet_vtx3deL) ;
    NNvectorVar_.push_back(Jet_numDaughters_pt03) ;//this variable has changed order, in bdt it was last, check why
    NNvectorVar_.push_back(Jet_energyRing_dR0_em_Jet_e) ;//15
    NNvectorVar_.push_back(Jet_energyRing_dR1_em_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR2_em_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR3_em_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR4_em_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR0_neut_Jet_e) ;//20
    NNvectorVar_.push_back(Jet_energyRing_dR1_neut_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR2_neut_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR3_neut_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR4_neut_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR0_ch_Jet_e) ;//25
    NNvectorVar_.push_back(Jet_energyRing_dR1_ch_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR2_ch_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR3_ch_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR4_ch_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR0_mu_Jet_e) ;//30
    NNvectorVar_.push_back(Jet_energyRing_dR1_mu_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR2_mu_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR3_mu_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR4_mu_Jet_e) ;
    NNvectorVar_.push_back(Jet_chHEF);//35
    NNvectorVar_.push_back(Jet_chEmEF);
    NNvectorVar_.push_back(Jet_leptonPtRelInv);
    NNvectorVar_.push_back(isEle);
    NNvectorVar_.push_back(isMu);
    NNvectorVar_.push_back(isOther);//40
    NNvectorVar_.push_back(Jet_mass);
    NNvectorVar_.push_back(Jet_withPtd);

}
    
std::vector<float> bRegressionProducer::EvaluateNN(){
    tensorflow::Tensor input(tensorflow::DT_FLOAT, {1,43});
    for (unsigned int i = 0; i < NNvectorVar_.size(); i++){
        input.matrix<float>()(0,i) =  float(NNvectorVar_[i]);
    }
    std::vector<tensorflow::Tensor> outputs;
    tensorflow::run(session, { { "ffwd_inp:0",input } }, { "ffwd_out/BiasAdd:0" }, &outputs);

    std::vector<float> correction(3);//3 outputs, first value is mean and then other 2 quantiles
    for (unsigned int i = 0; i < 3; i++) {
        correction[i] = outputs[0].matrix<float>()(0, i);
    }
        
    return correction;
}//end EvaluateNN
    
DEFINE_FWK_MODULE( bRegressionProducer );
