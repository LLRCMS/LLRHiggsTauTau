import FWCore.ParameterSet.Config as cms

#TRIGGERLIST=cms.vstring()
TRIGGERLIST=[]
#list triggers and filter paths here!
# channel: kemu=0, ketau=1,kmutau=2,ktautau=3
HLTLIST = cms.VPSet(

### ===================== MC, Spring15 ============================
    ### e mu 
    cms.PSet (
        HLT = cms.string("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1"),
        path1 = cms.vstring ("hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"), #ele filters
        path2 = cms.vstring ("hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23"), # mu filters
        channel = cms.int32(0)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v1"),
        path1 = cms.vstring ("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"), # ele filters
        path2 = cms.vstring ("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8"), # muon filters
        channel = cms.int32(0)
        ),
    cms.PSet ( #NOT IN BASELINE
        HLT = cms.string("HLT_IsoMu24_eta2p1_v1"),
        path1 = cms.vstring (""), # ele filters
        path2 = cms.vstring ("hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09"), # muon filters
        channel = cms.int32(0)
        ),
    cms.PSet ( #NOT IN BASELINE
        HLT = cms.string("HLT_IsoMu27_v1"),
        path1 = cms.vstring (""), # ele filters
        path2 = cms.vstring ("hltL3crIsoL1sMu25L1f0L2f10QL3f27QL3trkIsoFiltered0p09"), # muon filters
        channel = cms.int32(0)
        ),

    ### e tauh
    cms.PSet (
        HLT = cms.string("HLT_Ele22_eta2p1_WP75_Gsf_LooseIsoPFTau20_v1"),
        path1 = cms.vstring ("hltEle22WP75L1IsoEG20erTau20erGsfTrackIsoFilter", "hltOverlapFilterIsoEle22WP75GsfLooseIsoPFTau20"), # e filters
        path2 = cms.vstring ("hltPFTau20TrackLooseIso", "hltOverlapFilterIsoEle22WP75GsfLooseIsoPFTau20"), # tauh filters
        channel = cms.int32(1)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Ele32_eta2p1_WP75_Gsf_v1"),
        path1 = cms.vstring ("hltEle32WP75GsfTrackIsoFilter"), # e filters
        path2 = cms.vstring (""), # tauh filters
        channel = cms.int32(1)
        ),

  ### mu tauh
    cms.PSet (
        HLT = cms.string("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v1"),
        path1 = cms.vstring ("hltL3crIsoL1sMu16erTauJet20erL1f0L2f10QL3f17QL3trkIsoFiltered0p09", "hltOverlapFilterIsoMu17LooseIsoPFTau20"), # mu filters
        path2 = cms.vstring ("hltPFTau20TrackLooseIsoAgainstMuon", "hltOverlapFilterIsoMu17LooseIsoPFTau20"), # tauh filters
        channel = cms.int32(2)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu24_eta2p1_v1"),
        path1 = cms.vstring ("hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09"), # mu filters
        path2 = cms.vstring (""), # tauh filters
        channel = cms.int32(2)
        ),
    cms.PSet ( #NOT IN BASELINE
        HLT = cms.string("HLT_IsoMu27_v1"),
        path1 = cms.vstring ("hltL3crIsoL1sMu25L1f0L2f10QL3f27QL3trkIsoFiltered0p09"), # mu filters
        path2 = cms.vstring (""), # tauh filters
        channel = cms.int32(2)
        ),

### tauh tauh
    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v1"),
        path1 = cms.vstring ("hltDoublePFTau40TrackPt1MediumIsolationDz02Reg"), # tauh filters
        path2 = cms.vstring ("hltDoublePFTau40TrackPt1MediumIsolationDz02Reg"), # tauh filters (replicated)
        channel = cms.int32(3)
        ),

### ===================== DATA, August 2015 ============================
    ### e mu 
    cms.PSet (
        HLT = cms.string("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2"),
        path1 = cms.vstring ("hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"), #ele filters
        path2 = cms.vstring ("hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23"), # mu filters
        channel = cms.int32(0)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v2"),
        path1 = cms.vstring ("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"), # ele filters
        path2 = cms.vstring ("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8"), # muon filters
        channel = cms.int32(0)
        ),

    ### e tauh
    cms.PSet (
        HLT = cms.string("HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v1"),
        path1 = cms.vstring ("hltEle22WP75L1IsoEG20erTau20erGsfTrackIsoFilter", "hltOverlapFilterIsoEle22WP75GsfLooseIsoPFTau20"), # e filters
        path2 = cms.vstring ("hltPFTau20TrackLooseIso", "hltOverlapFilterIsoEle22WP75GsfLooseIsoPFTau20"), # tauh filters
        channel = cms.int32(1)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Ele32_eta2p1_WPTight_Gsf_v1"),
        path1 = cms.vstring ("hltEle32WP75GsfTrackIsoFilter"), # e filters
        path2 = cms.vstring (""), # tauh filters
        channel = cms.int32(1)
        ),

    ### mu tauh
    cms.PSet (
        HLT = cms.string("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v2"),
        path1 = cms.vstring ("hltL3crIsoL1sMu16erTauJet20erL1f0L2f10QL3f17QL3trkIsoFiltered0p09", "hltOverlapFilterIsoMu17LooseIsoPFTau20"), # mu filters
        path2 = cms.vstring ("hltPFTau20TrackLooseIsoAgainstMuon", "hltOverlapFilterIsoMu17LooseIsoPFTau20"), # tauh filters
        channel = cms.int32(2)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu24_eta2p1_IterTrk02_v2"),
        path1 = cms.vstring ("hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09"), # mu filters
        path2 = cms.vstring (""), # tauh filters
        channel = cms.int32(2)
        ),
    
    ### tauh tauh
    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v2"),
        path1 = cms.vstring ("hltDoublePFTau40TrackPt1MediumIsolationDz02Reg"), # tauh filters
        path2 = cms.vstring ("hltDoublePFTau40TrackPt1MediumIsolationDz02Reg"), # tauh filters (replicated)
        channel = cms.int32(3)
        )
    )



#now I create the trigger list for HLTconfig
for i in range(len(HLTLIST)):
    tmpl =  str(HLTLIST[i].HLT).replace('cms.string(\'','') ##CMSSW Vaffanculo
    tmpl = tmpl.replace('\')','') ##CMSSW Vaffanculo x 2
    TRIGGERLIST.append(tmpl)
#print TRIGGERLIST