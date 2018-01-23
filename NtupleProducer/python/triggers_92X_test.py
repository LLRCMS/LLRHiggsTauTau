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
        HLT = cms.string("HLT_IsoMu27_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07"),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(999)
        ),
### === Single electron triggers
    cms.PSet (
        HLT = cms.string("HLT_Ele35_WPTight_Gsf_v"),
        path1 = cms.vstring ("hltEle35noerWPTightGsfTrackIsoFilter"),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(999)
        ),
### === single tau triggers
    cms.PSet (
        HLT = cms.string("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v"),
        path1 = cms.vstring ("hltSelectedPFTau180MediumChargedIsolationL1HLTMatched"),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(999)
        ),
### === mu tauh triggers
    cms.PSet (
        HLT = cms.string("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07","hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded"),
        path2 = cms.vstring ("hltSelectedPFTau27LooseChargedIsolationAgainstMuonL1HLTMatched","hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
        ),
### === ele tauh triggers
    cms.PSet (
        HLT = cms.string("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v"),
        path1 = cms.vstring ("hltEle24erWPTightGsfTrackIsoFilterForTau","hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30"),
        path2 = cms.vstring ("hltSelectedPFTau30LooseChargedIsolationL1HLTMatched","hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(15)
        ),
### === tauh tauh triggers
    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau35TrackPt1MediumChargedIsolationDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau35TrackPt1MediumChargedIsolationDz02Reg"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),
## === VBF + double-tau triggers
    cms.PSet (
        HLT = cms.string("HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau20TrackPt1LooseChargedIsolationReg"),                    # hadronic tau
        path2 = cms.vstring ("hltDoublePFTau20TrackPt1LooseChargedIsolationReg"),                    # hadronic tau
        path3 = cms.vstring ("hltMatchedVBFTwoPFJets2CrossCleanedFromDoubleLooseChargedIsoPFTau20"), # 2 jets with pt>40
        path4 = cms.vstring ("hltMatchedVBFOnePFJet2CrossCleanedFromDoubleLooseChargedIsoPFTau20"),  # 1 jet with pt>115
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),
    #cms.PSet (
    #    HLT = cms.string("HLT_VBF_DoubleMediumChargedIsoPFTau20_Trk1_eta2p1_Reg_v"),
    #    path1 = cms.vstring ("hltDoublePFTau20TrackPt1MediumChargedIsolationReg"),                    # hadronic tau
    #    path2 = cms.vstring ("hltDoublePFTau20TrackPt1MediumChargedIsolationReg"),                    # hadronic tau
    #    path3 = cms.vstring ("hltMatchedVBFTwoPFJets2CrossCleanedFromDoubleMediumChargedIsoPFTau20"), # 2 jets with pt>40
    #    path4 = cms.vstring ("hltMatchedVBFOnePFJet2CrossCleanedFromDoubleMediumChargedIsoPFTau20"),  # 1 jet with pt>115
    #    leg1 = cms.int32(15),
    #    leg2 = cms.int32(15)
    #    ),
    #cms.PSet (
    #    HLT = cms.string("HLT_VBF_DoubleTightChargedIsoPFTau20_Trk1_eta2p1_Reg_v"),
    #    path1 = cms.vstring ("hltDoublePFTau20TrackPt1TightChargedIsolationReg"),                    # hadronic tau
    #    path2 = cms.vstring ("hltDoublePFTau20TrackPt1TightChargedIsolationReg"),                    # hadronic tau
    #    path3 = cms.vstring ("hltMatchedVBFTwoPFJets2CrossCleanedFromDoubleTightChargedIsoPFTau20"), # 2 jets with pt>40
    #    path4 = cms.vstring ("hltMatchedVBFOnePFJet2CrossCleanedFromDoubleTightChargedIsoPFTau20"),  # 1 jet with pt>115
    #    leg1 = cms.int32(15),
    #    leg2 = cms.int32(15)
    #    ),
## === 4 jets (for VBF)
    cms.PSet (
        HLT = cms.string("HLT_QuadPFJet103_88_75_15_BTagCSV_p013_VBF2_v7"), #FIXME
        path1 = cms.vstring (""), # tau filters
        path2 = cms.vstring (""), # tau filters
        path3 = cms.vstring ("hltVBFPFJetCSVSortedMqq460Detaqq3p5"), # 4 jets
        path4 = cms.vstring (""),
        leg1 = cms.int32(999),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_QuadPFJet103_88_75_15_DoubleBTagCSV_p013_p08_VBF1_v7"),  #FIXME
        path1 = cms.vstring (""), # tau filters
        path2 = cms.vstring (""), # tau filters
        path3 = cms.vstring ("hltVBFPFJetCSVSortedMqq200Detaqq1p5"), # 6 jets
        path4 = cms.vstring (""),  # 1 jet with pt>115
        leg1 = cms.int32(999),
        leg2 = cms.int32(999)
        ),
    )

#now I create the trigger list for HLTconfig
for i in range(len(HLTLIST)):
    tmpl =  str(HLTLIST[i].HLT).replace('cms.string(\'','') ##CMSSW Vaffanculo
    tmpl = tmpl.replace('\')','') ##CMSSW Vaffanculo x 2
    TRIGGERLIST.append(tmpl)
#print TRIGGERLIST
