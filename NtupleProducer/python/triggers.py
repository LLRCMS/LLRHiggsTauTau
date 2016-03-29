import FWCore.ParameterSet.Config as cms

## TRIGGER VERSIONING: insert paths without the explicit version number, e.g. HLT_anything_v 
## do not remove any other part of path name or the check will give meaningless results!

## leg numbering: write the object type the two legs refer to
## 11: electrons, 13: muons, 15: taus, 999: others/unused (e.g. leg 2 in single muon)
##

TRIGGERLIST=[]
#list triggers and filter paths here!
# channel: kemu=0, ketau=1,kmutau=2,ktautau=3
HLTLIST = cms.VPSet(

### === Single muon triggers
    cms.PSet (
        HLT = cms.string("HLT_IsoMu24_eta2p1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu27_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu25L1f0L2f10QL3f27QL3trkIsoFiltered0p09"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu17_eta2p1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sSingleMu16erL1f0L2f10QL3f17QL3trkIsoFiltered0p09"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu18_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu22_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(999)
        ),
### === Single electron triggers
    cms.PSet (
        HLT = cms.string("HLT_Ele32_eta2p1_WP75_Gsf_v"),
        path1 = cms.vstring ("hltEle32WP75GsfTrackIsoFilter"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Ele22_eta2p1_WP75_Gsf_v"),
        path1 = cms.vstring ("hltSingleEle22WP75GsfTrackIsoFilter"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Ele32_eta2p1_WPTight_Gsf_v"),
        path1 = cms.vstring ("hltEle32WPTightGsfTrackIsoFilter"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Ele23_WPLoose_Gsf_v"),
        path1 = cms.vstring ("hltEle23WPLooseGsfTrackIsoFilter"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(999)
        ),
### === mu tauh triggers
    cms.PSet (
        HLT = cms.string("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu16erTauJet20erL1f0L2f10QL3f17QL3trkIsoFiltered0p09", "hltOverlapFilterIsoMu17LooseIsoPFTau20"),
        path2 = cms.vstring ("hltPFTau20TrackLooseIsoAgainstMuon", "hltOverlapFilterIsoMu17LooseIsoPFTau20"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
        ),
### === ele tauh triggers
    cms.PSet (
        HLT = cms.string("HLT_Ele22_eta2p1_WP75_Gsf_LooseIsoPFTau20_v"),
        path1 = cms.vstring ("hltEle22WP75L1IsoEG20erTau20erGsfTrackIsoFilter", "hltOverlapFilterIsoEle22WP75GsfLooseIsoPFTau20"),
        path2 = cms.vstring ("hltPFTau20TrackLooseIso", "hltOverlapFilterIsoEle22WP75GsfLooseIsoPFTau20"),
        leg1 = cms.int32(11),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v"),
        path1 = cms.vstring ("hltSingleEle22WPLooseGsfTrackIsoFilter", "hltOverlapFilterIsoEle22WPLooseGsfLooseIsoPFTau20"),
        path2 = cms.vstring ("hltPFTau20TrackLooseIso", "hltOverlapFilterIsoEle22WPLooseGsfLooseIsoPFTau20"),
        leg1 = cms.int32(11),
        leg2 = cms.int32(15)
        ),
### === tauh tauh triggers
    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau35TrackPt1MediumIsolationDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau35TrackPt1MediumIsolationDz02Reg"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau40TrackPt1MediumIsolationDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau40TrackPt1MediumIsolationDz02Reg"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),
### === ele mu triggers
    cms.PSet (
        HLT = cms.string("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v"),
        path1 = cms.vstring ("hltMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered17"),
        path2 = cms.vstring ("hltMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"),
        leg1 = cms.int32(11),
        leg2 = cms.int32(13)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v"),
        path1 = cms.vstring ("hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8"),
        path2 = cms.vstring ("hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"),
        leg1 = cms.int32(11),
        leg2 = cms.int32(13)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v"),
        path1 = cms.vstring ("hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23"),
        path2 = cms.vstring ("hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"),
        leg1 = cms.int32(11),
        leg2 = cms.int32(13)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v"),
        path1 = cms.vstring ("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8"),
        path2 = cms.vstring ("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"),
        leg1 = cms.int32(11),
        leg2 = cms.int32(13)
        ),
### === mu mu triggers 
    cms.PSet (
        HLT = cms.string("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(13)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(13)
        ),
    cms.PSet (
        HLT = cms.string("HLT_TripleMu_12_10_5_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(13)
        ),
### === ele ele triggers 
    cms.PSet (
        HLT = cms.string("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(13)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(13)
        )
    )


#now I create the trigger list for HLTconfig
for i in range(len(HLTLIST)):
    tmpl =  str(HLTLIST[i].HLT).replace('cms.string(\'','') ##CMSSW Vaffanculo
    tmpl = tmpl.replace('\')','') ##CMSSW Vaffanculo x 2
    TRIGGERLIST.append(tmpl)
#print TRIGGERLIST
