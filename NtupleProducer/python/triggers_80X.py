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
        HLT = cms.string("HLT_IsoMu20_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18L1f0L2f10QL3f20QL3trkIsoFiltered0p09"),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoTkMu20_v"),
        path1 = cms.vstring ("hltL3fL1sMu18L1f0Tkf20QL3trkIsoFiltered0p09"),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu22_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09"),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoTkMu22_v"),
        path1 = cms.vstring ("hltL3fL1sMu20L1f0Tkf22QL3trkIsoFiltered0p09"),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu22_eta2p1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sSingleMu20erL1f0L2f10QL3f22QL3trkIsoFiltered0p09"),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu24_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09"),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(999)
        ), ### ok
    cms.PSet (
        HLT = cms.string("HLT_IsoTkMu22_eta2p1_v"),
        path1 = cms.vstring ("hltL3fL1sMu20erL1f0Tkf22QL3trkIsoFiltered0p09"),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoTkMu24_v"),
        path1 = cms.vstring ("hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09"),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoTkMu24_eta2p1_v"),
        path1 = cms.vstring ("hltL3fL1sMu22erL1f0Tkf24QL3trkIsoFiltered0p09"),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoTkMu27_v"),
        path1 = cms.vstring ("hltL3fL1sMu22Or25L1f0Tkf27QL3trkIsoFiltered0p09"),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(999)
        ),
### === Single electron triggers
    cms.PSet (
        HLT = cms.string("HLT_Ele25_WPTight_Gsf_v"),
        path1 = cms.vstring ("hltEle25WPTightGsfTrackIsoFilter"), #FIXME: to check
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Ele25_eta2p1_WPLoose_Gsf_v"),
        path1 = cms.vstring ("hltEle25erWPLooseGsfTrackIsoFilter"), #FIXME: to check
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
         HLT = cms.string("HLT_Ele23_WPLoose_Gsf_v"),
         path1 = cms.vstring ("hltEle23WPLooseGsfTrackIsoFilter"),
         path2 = cms.vstring (""),
         path3 = cms.vstring (""),
         path4 = cms.vstring (""),
         leg1 = cms.int32(11),
         leg2 = cms.int32(999)
         ),
### === ele tauh triggers
    cms.PSet (
        HLT = cms.string("HLT_Ele25_eta2p1_WPTight_Gsf_v"),
        path1 = cms.vstring ("hltEle25erWPTightGsfTrackIsoFilter"), #FIXME: to check
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Ele27_WPTight_Gsf_v"),
        path1 = cms.vstring ("hltEle27WPTightGsfTrackIsoFilter"), #FIXME: to check
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(999)
        ), ### ok
    cms.PSet (
        HLT = cms.string("HLT_Ele27_eta2p1_WPLoose_Gsf_v"),
        path1 = cms.vstring ("hltEle27erWPLooseGsfTrackIsoFilter"), #FIXME: to check
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Ele27_eta2p1_WPTight_Gsf_v"),
        path1 = cms.vstring ("hltEle27erWPTightGsfTrackIsoFilter"), #FIXME: to check
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(999)
        ), ### ok
### === mu tauh triggers
    cms.PSet (
        HLT = cms.string("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu16erTauJet20erL1f0L2f10QL3f17QL3trkIsoFiltered0p09", "hltOverlapFilterIsoMu17LooseIsoPFTau20"),
        path2 = cms.vstring ("hltPFTau20TrackLooseIsoAgainstMuon", "hltOverlapFilterIsoMu17LooseIsoPFTau20"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sSingleMu16erL1f0L2f10QL3f17QL3trkIsoFiltered0p09", "hltOverlapFilterSingleIsoMu17LooseIsoPFTau20"),
        path2 = cms.vstring ("hltPFTau20TrackLooseIsoAgainstMuon", "hltOverlapFilterSingleIsoMu17LooseIsoPFTau20"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
        ), ### ok
    cms.PSet (
        HLT = cms.string("HLT_IsoMu19_eta2p1_LooseIsoPFTau20_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erTauJet20erL1f0L2f10QL3f19QL3trkIsoFiltered0p09", "hltOverlapFilterIsoMu19LooseIsoPFTau20"),
        path2 = cms.vstring ("hltPFTau20TrackLooseIsoAgainstMuon", "hltOverlapFilterIsoMu19LooseIsoPFTau20"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
        ), ### ok
    cms.PSet (
        HLT = cms.string("HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sSingleMu18erIorSingleMu20erL1f0L2f10QL3f19QL3trkIsoFiltered0p09", "hltOverlapFilterSingleIsoMu19LooseIsoPFTau20"),
        path2 = cms.vstring ("hltPFTau20TrackLooseIsoAgainstMuon", "hltOverlapFilterSingleIsoMu19LooseIsoPFTau20"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
        ), ### ok
    cms.PSet (
        HLT = cms.string("HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sSingleMu20erIorSingleMu22erL1f0L2f10QL3f21QL3trkIsoFiltered0p09", "hltOverlapFilterSingleIsoMu21LooseIsoPFTau20"),
        path2 = cms.vstring ("hltPFTau20TrackLooseIsoAgainstMuon", "hltOverlapFilterSingleIsoMu21LooseIsoPFTau20"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
        ), ### ok
### === ele tauh triggers
    cms.PSet (
        HLT = cms.string("HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v"),
        path1 = cms.vstring ("hltEle22WPLooseL1SingleIsoEG20erGsfTrackIsoFilter", "hltOverlapFilterSingleIsoEle22WPLooseGsfLooseIsoPFTau20"),
        path2 = cms.vstring ("hltPFTau20TrackLooseIso", "hltOverlapFilterSingleIsoEle22WPLooseGsfLooseIsoPFTau20"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(15)
        ), ### ok
    cms.PSet (
        HLT = cms.string("HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v"),
        path1 = cms.vstring ("hltEle24WPLooseL1IsoEG22erTau20erGsfTrackIsoFilter", "hltOverlapFilterIsoEle24WPLooseGsfLooseIsoPFTau20"),
        path2 = cms.vstring ("hltPFTau20TrackLooseIso", "hltOverlapFilterIsoEle24WPLooseGsfLooseIsoPFTau20"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v"),
        path1 = cms.vstring ("hltEle24WPLooseL1SingleIsoEG22erGsfTrackIsoFilter", "hltOverlapFilterSingleIsoEle24WPLooseGsfLooseIsoPFTau20"),
        path2 = cms.vstring ("hltPFTau20TrackLooseIso", "hltOverlapFilterSingleIsoEle24WPLooseGsfLooseIsoPFTau20"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(15)
        ), ### ok
    cms.PSet (
        HLT = cms.string("HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau30_v"),
        path1 = cms.vstring ("hltEle24WPLooseL1IsoEG22erIsoTau26erGsfTrackIsoFilter", "hltOverlapFilterIsoEle24WPLooseGsfLooseIsoPFTau30"),
        path2 = cms.vstring ("hltPFTau30TrackLooseIso", "hltOverlapFilterIsoEle24WPLooseGsfLooseIsoPFTau30"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(15)
        ), ### ok
### === tauh tauh triggers
    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumIsoPFTau32_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau32TrackPt1MediumIsolationDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau32TrackPt1MediumIsolationDz02Reg"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau35TrackPt1MediumIsolationDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau35TrackPt1MediumIsolationDz02Reg"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau40TrackPt1MediumIsolationDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau40TrackPt1MediumIsolationDz02Reg"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau35TrackPt1MediumCombinedIsolationDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau35TrackPt1MediumCombinedIsolationDz02Reg"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau40TrackPt1MediumCombinedIsolationDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau40TrackPt1MediumCombinedIsolationDz02Reg"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v"),
        path1 = cms.vstring ("hltDoublePFTau40TrackPt1MediumCombinedIsolationDz02"),
        path2 = cms.vstring ("hltDoublePFTau40TrackPt1MediumCombinedIsolationDz02"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),
# ### === ele mu triggers
     cms.PSet (
         HLT = cms.string("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v"),
         path1 = cms.vstring ("hltMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered17"),
         path2 = cms.vstring ("hltMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"),
         path3 = cms.vstring (""),
         path4 = cms.vstring (""),
         leg1 = cms.int32(13),
         leg2 = cms.int32(11)
         ),
     cms.PSet (
         HLT = cms.string("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v"),
         path1 = cms.vstring ("hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8"),
         path2 = cms.vstring ("hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"),
         path3 = cms.vstring (""),
         path4 = cms.vstring (""),
         leg1 = cms.int32(13),
         leg2 = cms.int32(11)
         ),
    cms.PSet (
         HLT = cms.string("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v"),
         path1 = cms.vstring ("hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23"),
         path2 = cms.vstring ("hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"),
         path3 = cms.vstring (""),
         path4 = cms.vstring (""),
         leg1 = cms.int32(13),
         leg2 = cms.int32(11)
         ),
    cms.PSet (
         HLT = cms.string("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v"),
         path1 = cms.vstring ("hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23"),
         path2 = cms.vstring ("hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"),
         path3 = cms.vstring (""),
         path4 = cms.vstring (""),
         leg1 = cms.int32(13),
         leg2 = cms.int32(11)
         ),
    cms.PSet (
         HLT = cms.string("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v"),
         path1 = cms.vstring ("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8"),
         path2 = cms.vstring ("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"),
         path3 = cms.vstring (""),
         path4 = cms.vstring (""),
         leg1 = cms.int32(13),
         leg2 = cms.int32(11)
         ),
    cms.PSet (
         HLT = cms.string("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"),
         path1 = cms.vstring ("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8"),
         path2 = cms.vstring ("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"),
         path3 = cms.vstring (""),
         path4 = cms.vstring (""),
         leg1 = cms.int32(13),
         leg2 = cms.int32(11)
         ),
### === mu mu triggers 
### FIXME!! MuMu and EleEle paths have not been checked in 80X and filter names are dummy
    cms.PSet (
        HLT = cms.string("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(13)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(13)
        ),
    cms.PSet (
        HLT = cms.string("HLT_TripleMu_12_10_5_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(13)
        ),
### === ele ele triggers 
    cms.PSet (
        HLT = cms.string("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(11)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(11)
        ),    
    cms.PSet (
        HLT = cms.string("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(13)
        ),
### === 3 lepton triggers 
    cms.PSet (
        HLT = cms.string("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(13)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(11)
        ),
### === single tau triggers
    cms.PSet (
        HLT = cms.string("HLT_VLooseIsoPFTau140_Trk50_eta2p1_v"),
        path1 = cms.vstring ("hltPFTau140TrackPt50LooseAbsOrRelVLooseIso"),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_VLooseIsoPFTau120_Trk50_eta2p1_v"),
        path1 = cms.vstring ("hltPFTau120TrackPt50LooseAbsOrRelVLooseIso"),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(999)
        ),

### === MET triggers
    cms.PSet (
        HLT = cms.string("HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(999),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(999),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(999),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(999),
        leg2 = cms.int32(999)
        ),

### === Double muon triggers
    cms.PSet (
        HLT = cms.string("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(13)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(13)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(13)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(13)
        ),
    cms.PSet (
        HLT = cms.string("HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(13)
        ),
    cms.PSet (
        HLT = cms.string("HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(13)
        ),

    )

#now I create the trigger list for HLTconfig
for i in range(len(HLTLIST)):
    tmpl =  str(HLTLIST[i].HLT).replace('cms.string(\'','') ##CMSSW Vaffanculo
    tmpl = tmpl.replace('\')','') ##CMSSW Vaffanculo x 2
    TRIGGERLIST.append(tmpl)
#print TRIGGERLIST
