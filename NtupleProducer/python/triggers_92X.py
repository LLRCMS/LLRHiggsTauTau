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

### === Single muon triggers - OK
    cms.PSet (
        HLT = cms.string("HLT_IsoMu27_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07"),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(999),
        pt1 = cms.double(29),
        pt2 = cms.double(-999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu24_v"), # ---------- it's prescaled: TO BE REMOVED
        path1 = cms.vstring ("hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07"),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(999),
        pt1 = cms.double(26),
        pt2 = cms.double(-999)
        ),
### === Single electron triggers - OK
    cms.PSet (
        HLT = cms.string("HLT_Ele32_WPTight_Gsf_v"),
        path1 = cms.vstring ("hltEle32WPTightGsfTrackIsoFilter"),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(999),
        pt1 = cms.double(35),
        pt2 = cms.double(-999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Ele35_WPTight_Gsf_v"),
        path1 = cms.vstring ("hltEle35noerWPTightGsfTrackIsoFilter"),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(999),
        pt1 = cms.double(38),
        pt2 = cms.double(-999)
        ),
### emulated HLT_Ele32_WPTight_Gsf_v for Data// HLT_Ele32_WPTight_Gsf_v was online only from Run2017C
    cms.PSet (
        HLT = cms.string("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v"),
        path1 = cms.vstring ("hltEle32L1DoubleEGWPTightGsfTrackIsoFilter","hltEGL1SingleEGOrFilter"),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(999),
        pt1 = cms.double(35),
        pt2 = cms.double(-999)
        ),

### === mu tauh triggers - OK
    cms.PSet (
        HLT = cms.string("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07","hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded"),
        path2 = cms.vstring ("hltSelectedPFTau27LooseChargedIsolationAgainstMuonL1HLTMatched","hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15),
        pt1 = cms.double(22),
        pt2 = cms.double(32)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07","hltOverlapFilterIsoMu20MediumChargedIsoPFTau27L1Seeded"),
        path2 = cms.vstring ("hltSelectedPFTau27MediumChargedIsolationAgainstMuonL1HLTMatched","hltOverlapFilterIsoMu20MediumChargedIsoPFTau27L1Seeded"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15),
        pt1 = cms.double(22),
        pt2 = cms.double(32)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07","hltOverlapFilterIsoMu20TightChargedIsoPFTau27L1Seeded"),
        path2 = cms.vstring ("hltSelectedPFTau27TightChargedIsolationAgainstMuonL1HLTMatched","hltOverlapFilterIsoMu20TightChargedIsoPFTau27L1Seeded"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15),
        pt1 = cms.double(22),
        pt2 = cms.double(32)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07","hltOverlapFilterIsoMu20LooseChargedIsoTightOOSCPhotonsPFTau27L1Seeded"),
        path2 = cms.vstring ("hltSelectedPFTau27LooseChargedIsolationTightOOSCPhotonsAgainstMuonL1HLTMatched","hltOverlapFilterIsoMu20LooseChargedIsoTightOOSCPhotonsPFTau27L1Seeded"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15),
        pt1 = cms.double(22),
        pt2 = cms.double(32)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_TightID_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07","hltOverlapFilterIsoMu20MediumChargedIsoTightOOSCPhotonsPFTau27L1Seeded"),
        path2 = cms.vstring("hltSelectedPFTau27MediumChargedIsolationTightOOSCPhotonsAgainstMuonL1HLTMatched","hltOverlapFilterIsoMu20MediumChargedIsoTightOOSCPhotonsPFTau27L1Seeded"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15),
        pt1 = cms.double(22),
        pt2 = cms.double(32)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_TightID_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07","hltOverlapFilterIsoMu20TightChargedIsoTightOOSCPhotonsPFTau27L1Seeded"),
        path2 = cms.vstring ("hltSelectedPFTau27TightChargedIsolationTightOOSCPhotonsAgainstMuonL1HLTMatched","hltOverlapFilterIsoMu20TightChargedIsoTightOOSCPhotonsPFTau27L1Seeded"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15),
        pt1 = cms.double(22),
        pt2 = cms.double(32)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_eta2p1_SingleL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sSingleMu22erL1f0L2f10QL3f24QL3trkIsoFiltered0p07","hltOverlapFilterIsoMu24LooseChargedIsoPFTau20"),
        path2 = cms.vstring ("hltPFTau20TrackLooseChargedIsoAgainstMuon","hltOverlapFilterIsoMu24LooseChargedIsoPFTau20"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15),
        pt1 = cms.double(26),
        pt2 = cms.double(25)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_eta2p1_SingleL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sSingleMu22erL1f0L2f10QL3f24QL3trkIsoFiltered0p07","hltOverlapFilterIsoMu24MediumChargedIsoPFTau20"),
        path2 = cms.vstring ("hltPFTau20TrackMediumChargedIsoAgainstMuon","hltOverlapFilterIsoMu24MediumChargedIsoPFTau20"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15),
        pt1 = cms.double(26),
        pt2 = cms.double(25)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_eta2p1_SingleL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sSingleMu22erL1f0L2f10QL3f24QL3trkIsoFiltered0p07","hltOverlapFilterIsoMu24TightChargedIsoPFTau20"),
        path2 = cms.vstring ("hltPFTau20TrackTightChargedIsoAgainstMuon","hltOverlapFilterIsoMu24TightChargedIsoPFTau20"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15),
        pt1 = cms.double(26),
        pt2 = cms.double(25)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_eta2p1_TightID_SingleL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sSingleMu22erL1f0L2f10QL3f24QL3trkIsoFiltered0p07","hltOverlapFilterIsoMu24LooseChargedIsoTightOOSCPhotonsPFTau20"),
        path2 = cms.vstring ("hltPFTau20TrackLooseChargedIsoAgainstMuon","hltOverlapFilterIsoMu24LooseChargedIsoTightOOSCPhotonsPFTau20"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15),
        pt1 = cms.double(26),
        pt2 = cms.double(25)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_eta2p1_TightID_SingleL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sSingleMu22erL1f0L2f10QL3f24QL3trkIsoFiltered0p07","hltOverlapFilterIsoMu24MediumChargedIsoTightOOSCPhotonsPFTau20"),
        path2 = cms.vstring ("hltPFTau20TrackMediumChargedIsoAgainstMuon","hltOverlapFilterIsoMu24MediumChargedIsoTightOOSCPhotonsPFTau20"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15),
        pt1 = cms.double(26),
        pt2 = cms.double(25)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_eta2p1_TightID_SingleL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sSingleMu22erL1f0L2f10QL3f24QL3trkIsoFiltered0p07","hltOverlapFilterIsoMu24TightChargedIsoTightOOSCPhotonsPFTau20"),
        path2 = cms.vstring ("hltPFTau20TrackTightChargedIsoAgainstMuon","hltOverlapFilterIsoMu24TightChargedIsoTightOOSCPhotonsPFTau20"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15),
        pt1 = cms.double(26),
        pt2 = cms.double(25)
        ),
### === ele tauh triggers - OK
    cms.PSet (
        HLT = cms.string("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v"),
        path1 = cms.vstring ("hltEle24erWPTightGsfTrackIsoFilterForTau","hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30"),
        path2 = cms.vstring ("hltSelectedPFTau30LooseChargedIsolationL1HLTMatched","hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(15),
        pt1 = cms.double(27),
        pt2 = cms.double(35)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1_v"),
        path1 = cms.vstring ("hltEle24erWPTightGsfTrackIsoFilterForTau","hltOverlapFilterIsoEle24WPTightGsfMediumChargedIsoPFTau30"),
        path2 = cms.vstring ("hltSelectedPFTau30MediumChargedIsolationL1HLTMatched","hltOverlapFilterIsoEle24WPTightGsfMediumChargedIsoPFTau30"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(15),
        pt1 = cms.double(27),
        pt2 = cms.double(35)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1_v"),
        path1 = cms.vstring ("hltEle24erWPTightGsfTrackIsoFilterForTau","hltOverlapFilterIsoEle24WPTightGsfTightChargedIsoPFTau30"),
        path2 = cms.vstring ("hltSelectedPFTau30TightChargedIsolationL1HLTMatched","hltOverlapFilterIsoEle24WPTightGsfTightChargedIsoPFTau30"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(15),
        pt1 = cms.double(27),
        pt2 = cms.double(35)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_TightID_CrossL1_v"),
        path1 = cms.vstring ("hltEle24erWPTightGsfTrackIsoFilterForTau","hltOverlapFilterIsoEle24WPTightGsfLooseIsoTightOOSCPhotonsPFTau30"),
        path2 = cms.vstring ("hltSelectedPFTau30LooseChargedIsolationTightOOSCPhotonsL1HLTMatched","hltOverlapFilterIsoEle24WPTightGsfLooseIsoTightOOSCPhotonsPFTau30"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(15),
        pt1 = cms.double(27),
        pt2 = cms.double(35)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_TightID_CrossL1_v"),
        path1 = cms.vstring ("hltEle24erWPTightGsfTrackIsoFilterForTau","hltOverlapFilterIsoEle24WPTightGsfMediumIsoTightOOSCPhotonsPFTau30"),
        path2 = cms.vstring ("hltSelectedPFTau30MediumChargedIsolationTightOOSCPhotonsL1HLTMatched","hltOverlapFilterIsoEle24WPTightGsfMediumIsoTightOOSCPhotonsPFTau30"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(15),
        pt1 = cms.double(27),
        pt2 = cms.double(35)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_TightID_CrossL1_v"),
        path1 = cms.vstring ("hltEle24erWPTightGsfTrackIsoFilterForTau","hltOverlapFilterIsoEle24WPTightGsfTightIsoTightOOSCPhotonsPFTau30"),
        path2 = cms.vstring ("hltSelectedPFTau30TightChargedIsolationTightOOSCPhotonsL1HLTMatched","hltOverlapFilterIsoEle24WPTightGsfTightIsoTightOOSCPhotonsPFTau30"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(15),
        pt1 = cms.double(27),
        pt2 = cms.double(35)
        ),
### === tauh tauh triggers - OK
    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg_v"), # ---------- it's prescaled: TO BE REMOVED
        path1 = cms.vstring ("hltDoublePFTau35TrackPt1MediumChargedIsolationDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau35TrackPt1MediumChargedIsolationDz02Reg"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15),
        pt1 = cms.double(40),
        pt2 = cms.double(40)
        ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15),
        pt1 = cms.double(40),
        pt2 = cms.double(40)
        ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15),
        pt1 = cms.double(40),
        pt2 = cms.double(40)
        ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleTightChargedIsoPFTau35_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau35TrackPt1TightChargedIsolationDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau35TrackPt1TightChargedIsolationDz02Reg"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15),
        pt1 = cms.double(40),
        pt2 = cms.double(40)
        ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau40TrackPt1TightChargedIsolationDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau40TrackPt1TightChargedIsolationDz02Reg"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15),
        pt1 = cms.double(45),
        pt2 = cms.double(45)
        ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleTightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau40TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau40TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15),
        pt1 = cms.double(45),
        pt2 = cms.double(45)
        ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau40TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau40TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15),
        pt1 = cms.double(45),
        pt2 = cms.double(45)
        ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumChargedIsoPFTau40_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau40TrackPt1MediumChargedIsolationDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau40TrackPt1MediumChargedIsolationDz02Reg"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15),
        pt1 = cms.double(45),
        pt2 = cms.double(45)
        ),
# ### === ele mu triggers -- TO BE DONE
#     cms.PSet (
#         HLT = cms.string("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v"),
#         path1 = cms.vstring ("hltMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered17"),
#         path2 = cms.vstring ("hltMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"),
#         leg1 = cms.int32(13),
#         leg2 = cms.int32(11)
#         ),
#     cms.PSet (
#         HLT = cms.string("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v"),
#         path1 = cms.vstring ("hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8"),
#         path2 = cms.vstring ("hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"),
#         leg1 = cms.int32(13),
#         leg2 = cms.int32(11)
#         ),
#    cms.PSet (
#         HLT = cms.string("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v"),
#         path1 = cms.vstring ("hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23"),
#         path2 = cms.vstring ("hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"),
#         leg1 = cms.int32(13),
#         leg2 = cms.int32(11)
#         ),
#    cms.PSet (
#         HLT = cms.string("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v"),
#         path1 = cms.vstring ("hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23"),
#         path2 = cms.vstring ("hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"),
#         leg1 = cms.int32(13),
#         leg2 = cms.int32(11)
#         ),
#    cms.PSet (
#         HLT = cms.string("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v"),
#         path1 = cms.vstring ("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8"),
#         path2 = cms.vstring ("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"),
#         leg1 = cms.int32(13),
#         leg2 = cms.int32(11)
#         ),
#    cms.PSet (
#         HLT = cms.string("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"),
#         path1 = cms.vstring ("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8"),
#         path2 = cms.vstring ("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"),
#         leg1 = cms.int32(13),
#         leg2 = cms.int32(11)
#         ),

### === mu mu triggers  -- TO BE DONE
### FIXME!! MuMu and EleEle paths have not been checked in 80X and filter names are dummy
#    cms.PSet (
#        HLT = cms.string("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v"),
#        path1 = cms.vstring (""),
#        path2 = cms.vstring (""),
#        leg1 = cms.int32(13),
#        leg2 = cms.int32(13)
#        ),
#    cms.PSet (
#        HLT = cms.string("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"),
#        path1 = cms.vstring (""),
#        path2 = cms.vstring (""),
#        leg1 = cms.int32(13),
#        leg2 = cms.int32(13)
#        ),
#    cms.PSet (
#        HLT = cms.string("HLT_TripleMu_12_10_5_v"),
#        path1 = cms.vstring (""),
#        path2 = cms.vstring (""),
#        leg1 = cms.int32(13),
#        leg2 = cms.int32(13)
#        ),

### === ele ele triggers  -- TO BE DONE
#    cms.PSet (
#        HLT = cms.string("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"),
#        path1 = cms.vstring (""),
#        path2 = cms.vstring (""),
#        leg1 = cms.int32(11),
#        leg2 = cms.int32(11)
#        ),
#    cms.PSet (
#        HLT = cms.string("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"),
#        path1 = cms.vstring (""),
#        path2 = cms.vstring (""),
#        leg1 = cms.int32(11),
#        leg2 = cms.int32(11)
#        ),
#    cms.PSet (
#        HLT = cms.string("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v"),
#        path1 = cms.vstring (""),
#        path2 = cms.vstring (""),
#        leg1 = cms.int32(13),
#        leg2 = cms.int32(13)
#        ),

### === 3 lepton triggers  -- TO BE DONE
#    cms.PSet (
#        HLT = cms.string("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v"),
#        path1 = cms.vstring (""),
#        path2 = cms.vstring (""),
#        leg1 = cms.int32(13),
#        leg2 = cms.int32(13)
#        ),
#    cms.PSet (
#        HLT = cms.string("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v"),
#        path1 = cms.vstring (""),
#        path2 = cms.vstring (""),
#        leg1 = cms.int32(11),
#        leg2 = cms.int32(11)
#        ),
### === single tau triggers - OK
    cms.PSet (
        HLT = cms.string("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v"),
        path1 = cms.vstring ("hltSelectedPFTau180MediumChargedIsolationL1HLTMatched","hltPFTau180TrackPt50LooseAbsOrRelMediumHighPtRelaxedIsoIso"),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(999),
        pt1 = cms.double(185),
        pt2 = cms.double(-999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr_v"),
        path1 = cms.vstring ("hltSelectedPFTau180MediumChargedIsolationL1HLTMatched1Prong","hltPFTau180TrackPt50LooseAbsOrRelMediumHighPtRelaxedIso1Prong"),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(999),
        pt1 = cms.double(185),
        pt2 = cms.double(-999)
        ),
## === VBF + double-tau triggers
    cms.PSet (
        HLT = cms.string("HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau20TrackPt1LooseChargedIsolationReg"),                    # hadronic tau
        path2 = cms.vstring ("hltDoublePFTau20TrackPt1LooseChargedIsolationReg"),                    # hadronic tau
        path3 = cms.vstring ("hltMatchedVBFTwoPFJets2CrossCleanedFromDoubleLooseChargedIsoPFTau20"), # 2 jets with pt>40
        path4 = cms.vstring ("hltMatchedVBFOnePFJet2CrossCleanedFromDoubleLooseChargedIsoPFTau20"),  # 1 jet with pt>115
        leg1 = cms.int32(15),
        leg2 = cms.int32(15),
        pt1 = cms.double(25),
        pt2 = cms.double(25)
        ),
    cms.PSet (
        HLT = cms.string("HLT_VBF_DoubleMediumChargedIsoPFTau20_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau20TrackPt1MediumChargedIsolationReg"),                    # hadronic tau
        path2 = cms.vstring ("hltDoublePFTau20TrackPt1MediumChargedIsolationReg"),                    # hadronic tau
        path3 = cms.vstring ("hltMatchedVBFTwoPFJets2CrossCleanedFromDoubleMediumChargedIsoPFTau20"), # 2 jets with pt>40
        path4 = cms.vstring ("hltMatchedVBFOnePFJet2CrossCleanedFromDoubleMediumChargedIsoPFTau20"),  # 1 jet with pt>115
        leg1 = cms.int32(15),
        leg2 = cms.int32(15),
        pt1 = cms.double(25),
        pt2 = cms.double(25)
        ),
    cms.PSet (
        HLT = cms.string("HLT_VBF_DoubleTightChargedIsoPFTau20_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau20TrackPt1TightChargedIsolationReg"),                    # hadronic tau
        path2 = cms.vstring ("hltDoublePFTau20TrackPt1TightChargedIsolationReg"),                    # hadronic tau
        path3 = cms.vstring ("hltMatchedVBFTwoPFJets2CrossCleanedFromDoubleTightChargedIsoPFTau20"), # 2 jets with pt>40
        path4 = cms.vstring ("hltMatchedVBFOnePFJet2CrossCleanedFromDoubleTightChargedIsoPFTau20"),  # 1 jet with pt>115
        leg1 = cms.int32(15),
        leg2 = cms.int32(15),
        pt1 = cms.double(25),
        pt2 = cms.double(25)
        ),
## === Tau + MET
    cms.PSet (
        HLT = cms.string("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_v"),
        path1 = cms.vstring ("hltPFTau50TrackPt30MediumAbsOrRelIso1Prong", "hltSelectedPFTau50MediumChargedIsolationL1HLTMatched"),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15),
        pt1 = cms.double(55),
        pt2 = cms.double(-999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100_v"),
        path1 = cms.vstring ("hltPFTau50TrackPt30MediumAbsOrRelIso1Prong", "hltSelectedPFTau50MediumChargedIsolationL1HLTMatched"),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15),
        pt1 = cms.double(55),
        pt2 = cms.double(-999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110_v"),
        path1 = cms.vstring ("hltPFTau50TrackPt30MediumAbsOrRelIso1Prong", "hltSelectedPFTau50MediumChargedIsolationL1HLTMatched"),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15),
        pt1 = cms.double(55),
        pt2 = cms.double(-999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130_v"),
        path1 = cms.vstring ("hltPFTau50TrackPt30MediumAbsOrRelIso1Prong", "hltSelectedPFTau50MediumChargedIsolationL1HLTMatched"),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15),
        pt1 = cms.double(55),
        pt2 = cms.double(-999)
        ),

# === MET triggers
    cms.PSet (
        HLT = cms.string("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v"),
        path1 = cms.vstring ("",""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(999),
        leg2 = cms.int32(999),
        ),
    cms.PSet (
        HLT = cms.string("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v"),
        path1 = cms.vstring ("",""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(999),
        leg2 = cms.int32(999),
        ),

### === Double muon triggers
    cms.PSet (
        HLT = cms.string("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"), # unprescaled, included in run >= 299368
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(13),
        ),
    cms.PSet (
        HLT = cms.string("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v"), # unprescaled
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(13),
        ),

## === 4 jets (for VBF) # FILTERS TO BE FIXED
    #cms.PSet (
    #    HLT = cms.string("HLT_QuadPFJet103_88_75_15_BTagCSV_p013_VBF2_v7"), #FIXME
    #    path1 = cms.vstring (""), # tau filters
    #    path2 = cms.vstring (""), # tau filters
    #    path3 = cms.vstring ("hltVBFPFJetCSVSortedMqq460Detaqq3p5"), # 4 jets
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(999),
    #    leg2 = cms.int32(999)
    #    ),
    #cms.PSet (
    #    HLT = cms.string("HLT_QuadPFJet103_88_75_15_DoubleBTagCSV_p013_p08_VBF1_v7"),  #FIXME
    #    path1 = cms.vstring (""), # tau filters
    #    path2 = cms.vstring (""), # tau filters
    #    path3 = cms.vstring ("hltVBFPFJetCSVSortedMqq200Detaqq1p5"), # 6 jets
    #    path4 = cms.vstring (""),  # 1 jet with pt>115
    #    leg1 = cms.int32(999),
    #    leg2 = cms.int32(999)
    #    ),
    )

#now I create the trigger list for HLTconfig
for i in range(len(HLTLIST)):
    tmpl =  str(HLTLIST[i].HLT).replace('cms.string(\'','') ##CMSSW Vaffanculo
    tmpl = tmpl.replace('\')','') ##CMSSW Vaffanculo x 2
    TRIGGERLIST.append(tmpl)
#print TRIGGERLIST
