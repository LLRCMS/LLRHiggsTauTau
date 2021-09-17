import FWCore.ParameterSet.Config as cms

## TRIGGER VERSIONING: insert paths without the explicit version number, e.g. HLT_anything_v
## do not remove any other part of path name or the check will give meaningless results!

## leg numbering: write the object type the two legs refer to
## 11: electrons, 13: muons, 15: taus, 999: others/unused (e.g. leg 2 in sigle muon)
##

TRIGGERLIST=[]
#list triggers and filter paths here!
# channel: kemu=0, ketau=1,kmutau=2,ktautau=3
HLTLIST = cms.VPSet(

### === Single muon triggers   -- UPDATED FOR 2018 DATA
    cms.PSet (
        HLT = cms.string("HLT_IsoMu24_v"),
        path1 = cms.vstring ("hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07"),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu27_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07"),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(999)
        ),
### === Single electron triggers   -- UPDATED FOR 2018 DATA
    cms.PSet (
        HLT = cms.string("HLT_Ele32_WPTight_Gsf_v"),
        path1 = cms.vstring ("hltEle32WPTightGsfTrackIsoFilter"),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Ele35_WPTight_Gsf_v"),
        path1 = cms.vstring ("hltEle35noerWPTightGsfTrackIsoFilter"),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(999)
        ),

### === Single tauh triggers 
    cms.PSet (
        HLT = cms.string("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v"),
        path1 = cms.vstring ("hltSelectedPFTau180MediumChargedIsolationL1HLTMatched"),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_MediumChargedIsoPFTau200HighPtRelaxedIso_Trk50_eta2p1_v"),
        path1 = cms.vstring ("hltSelectedPFTau200MediumChargedIsolationL1HLTMatched"),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_MediumChargedIsoPFTau220HighPtRelaxedIso_Trk50_eta2p1_v"),
        path1 = cms.vstring ("hltSelectedPFTau220MediumChargedIsolationL1HLTMatched"),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr_v"),
        path1 = cms.vstring ("hltSelectedPFTau180MediumChargedIsolationL1HLTMatched1Prong"),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(999)
        ),

### === mu tauh triggers -- CHECK first muon filter, I think it's wrong on the TauTrigger twiki -- TauTrigger twiki has been corrected on 2019-04-29
    cms.PSet (
        HLT = cms.string("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoBigORMu18erTauXXer2p1L1f0L2f10QL3f20QL3trkIsoFiltered0p07","hltHpsOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded"),
        path2 = cms.vstring ("hltHpsSelectedPFTau27LooseChargedIsolationAgainstMuonL1HLTMatched","hltHpsOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoBigORMu18erTauXXer2p1L1f0L2f10QL3f20QL3trkIsoFiltered0p07","hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded"),
        path2 = cms.vstring ("hltSelectedPFTau27LooseChargedIsolationAgainstMuonL1HLTMatched","hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
        ),
### === ele tauh triggers  -- UPDATED FOR 2018 DATA
    cms.PSet (
        HLT = cms.string("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1_v"),
        path1 = cms.vstring ("hltEle24erWPTightGsfTrackIsoFilterForTau","hltHpsOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30"),
        path2 = cms.vstring ("hltHpsSelectedPFTau30LooseChargedIsolationL1HLTMatched","hltHpsOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v"),
        path1 = cms.vstring ("hltEle24erWPTightGsfTrackIsoFilterForTau","hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30"),
        path2 = cms.vstring ("hltSelectedPFTau30LooseChargedIsolationL1HLTMatched","hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(15)
        ),
### === tauh tauh triggers  -- UPDATED FOR 2018 DATA
    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltHpsDoublePFTau35TrackPt1MediumChargedIsolationDz02Reg"),
        path2 = cms.vstring ("hltHpsDoublePFTau35TrackPt1MediumChargedIsolationDz02Reg"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau40TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau40TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltDoublePFTau40TrackPt1TightChargedIsolationDz02Reg"),
        path2 = cms.vstring ("hltDoublePFTau40TrackPt1TightChargedIsolationDz02Reg"),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),

## === VBF + double-tau triggers -- Davide: These triggerslooks ok in HLT GUI
    cms.PSet (
        HLT = cms.string("HLT_VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1_v"),
        path1 = cms.vstring ("hltHpsDoublePFTau20TrackLooseChargedIsoAgainstMuon"),                    # hadronic tau
        path2 = cms.vstring ("hltHpsDoublePFTau20TrackLooseChargedIsoAgainstMuon"),                    # hadronic tau
        path3 = cms.vstring ("hltMatchedVBFTwoPFJets2CrossCleanedFromDoubleLooseChargedIsoPFTauHPS20"), # 2 jets with pt>40
        path4 = cms.vstring ("hltMatchedVBFOnePFJet2CrossCleanedFromDoubleLooseChargedIsoPFTauHPS20"),  # 1 jet with pt>115
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_v"),
        path1 = cms.vstring ("hltDoublePFTau20TrackLooseChargedIsoAgainstMuon"),
        path2 = cms.vstring ("hltDoublePFTau20TrackLooseChargedIsoAgainstMuon"),
        path3 = cms.vstring ("hltMatchedVBFTwoPFJets2CrossCleanedFromDoubleLooseChargedIsoPFTau20"),
        path4 = cms.vstring ("hltMatchedVBFOnePFJet2CrossCleanedFromDoubleLooseChargedIsoPFTau20"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),

## === Tau + MET -- UPDATED FOR 2018 DATA
    cms.PSet (
        HLT = cms.string("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_v"), # prescaled
        path1 = cms.vstring ("hltSelectedPFTau50MediumChargedIsolationL1HLTMatched"),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100_v"),
        path1 = cms.vstring ("hltSelectedPFTau50MediumChargedIsolationL1HLTMatched"),
        path2 = cms.vstring (""), 
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110_v"),
        path1 = cms.vstring ("hltSelectedPFTau50MediumChargedIsolationL1HLTMatched"),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130_v"),
        path1 = cms.vstring ("hltSelectedPFTau50MediumChargedIsolationL1HLTMatched"),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),
	
## Monitoring triggers for 2018 data taking -- Filters should be checked
#    cms.PSet (
#        HLT = cms.string("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1_v"),
#        path1 = cms.vstring ("hltL3crIsoL1sBigOrMuXXerIsoTauYYerL1f0L2f10QL3f24QL3trkIsoFiltered0p07","hltHpsOverlapFilterIsoMu24MediumChargedIsoPFTau35MonitoringReg"),
#        path2 = cms.vstring ("hltHpsSelectedPFTau35TrackPt1MediumChargedIsolationL1HLTMatchedReg","hltHpsOverlapFilterIsoMu24MediumChargedIsoPFTau35MonitoringReg"),
#        path3 = cms.vstring (""),
#        path4 = cms.vstring (""),
#        leg1 = cms.int32(13),
#        leg2 = cms.int32(15)
#        ),
#    cms.PSet (
#        HLT = cms.string("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1_v"),
#        path1 = cms.vstring ("hltL3crIsoL1sBigOrMuXXerIsoTauYYerL1f0L2f10QL3f24QL3trkIsoFiltered0p07","hltHpsOverlapFilterIsoMu24MediumChargedIsoAndTightOOSCPhotonsPFTau35MonitoringReg"),
#        path2 = cms.vstring ("hltHpsSelectedPFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsL1HLTMatchedReg","hltHpsOverlapFilterIsoMu24MediumChargedIsoAndTightOOSCPhotonsPFTau35MonitoringReg"),
#        path3 = cms.vstring (""),
#        path4 = cms.vstring (""),
#        leg1 = cms.int32(13),
#        leg2 = cms.int32(15)
#        ),
#    cms.PSet (
#        HLT = cms.string("HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1_v"),
#        path1 = cms.vstring ("hltL3crIsoL1sBigOrMuXXerIsoTauYYerL1f0L2f10QL3f24QL3trkIsoFiltered0p07","hltHpsOverlapFilterIsoMu24TightChargedIsoPFTau35MonitoringReg"),
#        path2 = cms.vstring ("hltHpsSelectedPFTau35TrackPt1TightChargedIsolationL1HLTMatchedReg","hltHpsOverlapFilterIsoMu24TightChargedIsoPFTau35MonitoringReg"),
#        path3 = cms.vstring (""),
#        path4 = cms.vstring (""),
#        leg1 = cms.int32(13),
#        leg2 = cms.int32(15)
#        ),
#    cms.PSet (
#        HLT = cms.string("HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1_v"),
#        path1 = cms.vstring ("hltL3crIsoL1sBigOrMuXXerIsoTauYYerL1f0L2f10QL3f24QL3trkIsoFiltered0p07","hltHpsOverlapFilterIsoMu24TightChargedIsoAndTightOOSCPhotonsPFTau35MonitoringReg"),
#        path2 = cms.vstring ("hltHpsSelectedPFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsL1HLTMatchedReg","hltHpsOverlapFilterIsoMu24TightChargedIsoAndTightOOSCPhotonsPFTau35MonitoringReg"),
#       path3 = cms.vstring (""),
#        path4 = cms.vstring (""),
#        leg1 = cms.int32(13),
#        leg2 = cms.int32(15)
#        ),
#    cms.PSet (
#        HLT = cms.string("HLT_IsoMu27_LooseChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_v"),
#        path1 = cms.vstring ("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07","hltHpsOverlapFilterIsoMu27LooseChargedIsoPFTau20"),
#        path2 = cms.vstring ("hltHpsPFTau20TrackLooseChargedIsoAgainstMuon","hltHpsOverlapFilterIsoMu27LooseChargedIsoPFTau20"),
#        path3 = cms.vstring (""),
#        path4 = cms.vstring (""),
#        leg1 = cms.int32(13),
#        leg2 = cms.int32(15)
#        ),
#    cms.PSet (
#        HLT = cms.string("HLT_IsoMu27_MediumChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_v"),
#        path1 = cms.vstring ("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07","hltHpsOverlapFilterIsoMu27MediumChargedIsoPFTau20"),
#        path2 = cms.vstring ("hltHpsPFTau20TrackMediumChargedIsoAgainstMuon","hltHpsOverlapFilterIsoMu27MediumChargedIsoPFTau20"),
#        path3 = cms.vstring (""),
#        path4 = cms.vstring (""),
#        leg1 = cms.int32(13),
#        leg2 = cms.int32(15)
#        ),
#    cms.PSet (
#        HLT = cms.string("HLT_IsoMu27_TightChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_v"),
#        path1 = cms.vstring ("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07","hltHpsOverlapFilterIsoMu27TightChargedIsoPFTau20"),
#        path2 = cms.vstring ("hltHpsPFTau20TrackTightChargedIsoAgainstMuon","hltHpsOverlapFilterIsoMu27TightChargedIsoPFTau20"),
#        path3 = cms.vstring (""),
#        path4 = cms.vstring (""),
#        leg1 = cms.int32(13),
#        leg2 = cms.int32(15)
#        ),
#    cms.PSet (
#        HLT = cms.string("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_v"),
#        path1 = cms.vstring ("hltL3crIsoL1sMu22erIsoTau40erL1f0L2f10QL3f24QL3trkIsoFiltered0p07",""),
#        path2 = cms.vstring ("hltSelectedPFTau50MediumChargedIsolationL1HLTMatchedMu22IsoTau40",""),
#        path3 = cms.vstring (""),
#        path4 = cms.vstring (""),
#        leg1 = cms.int32(13),
#        leg2 = cms.int32(15)
#        ),
    ## Duplicate trigger names without "HPS" - See https://twiki.cern.ch/twiki/bin/view/CMS/TauTrigger#Trigger_table_for_2018 - HPS deployed in Menu v2.2
#    cms.PSet (
#        HLT = cms.string("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v"),
#        path1 = cms.vstring ("hltL3crIsoL1sBigOrMuXXerIsoTauYYerL1f0L2f10QL3f24QL3trkIsoFiltered0p07","hltOverlapFilterIsoMu24MediumChargedIsoPFTau35MonitoringReg"),
#        path2 = cms.vstring ("hltSelectedPFTau35TrackPt1MediumChargedIsolationL1HLTMatchedReg","hltOverlapFilterIsoMu24MediumChargedIsoPFTau35MonitoringReg"),
#        path3 = cms.vstring (""),
#        path4 = cms.vstring (""),
#        leg1 = cms.int32(13),
#        leg2 = cms.int32(15)
#        ),
#    cms.PSet (
#        HLT = cms.string("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1_v"),
#        path1 = cms.vstring ("hltL3crIsoL1sBigOrMuXXerIsoTauYYerL1f0L2f10QL3f24QL3trkIsoFiltered0p07","hltOverlapFilterIsoMu24MediumChargedIsoAndTightOOSCPhotonsPFTau35MonitoringReg"),
#        path2 = cms.vstring ("hltOverlapFilterIsoMu24MediumChargedIsoAndTightOOSCPhotonsPFTau35MonitoringReg","hltOverlapFilterIsoMu24MediumChargedIsoAndTightOOSCPhotonsPFTau35MonitoringReg"),
#        path3 = cms.vstring (""),
#        path4 = cms.vstring (""),
#        leg1 = cms.int32(13),
#        leg2 = cms.int32(15)
#        ),
#    cms.PSet (
#        HLT = cms.string("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v"),
#        path1 = cms.vstring ("hltL3crIsoL1sBigOrMuXXerIsoTauYYerL1f0L2f10QL3f24QL3trkIsoFiltered0p07","hltOverlapFilterIsoMu24TightChargedIsoPFTau35MonitoringReg"),
#        path2 = cms.vstring ("hltSelectedPFTau35TrackPt1TightChargedIsolationL1HLTMatchedReg","hltOverlapFilterIsoMu24TightChargedIsoPFTau35MonitoringReg"),
#        path3 = cms.vstring (""),
#        path4 = cms.vstring (""),
#        leg1 = cms.int32(13),
#        leg2 = cms.int32(15)
#        ),
#    cms.PSet (
#        HLT = cms.string("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1_v"),
#        path1 = cms.vstring ("hltL3crIsoL1sBigOrMuXXerIsoTauYYerL1f0L2f10QL3f24QL3trkIsoFiltered0p07","hltOverlapFilterIsoMu24TightChargedIsoAndTightOOSCPhotonsPFTau35MonitoringReg"),
#        path2 = cms.vstring ("hltSelectedPFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsL1HLTMatchedReg","hltOverlapFilterIsoMu24TightChargedIsoAndTightOOSCPhotonsPFTau35MonitoringReg"),
#        path3 = cms.vstring (""),
#        path4 = cms.vstring (""),
#        leg1 = cms.int32(13),
#        leg2 = cms.int32(15)
#        ),
#    cms.PSet (
#        HLT = cms.string("HLT_IsoMu27_LooseChargedIsoPFTau20_Trk1_eta2p1_SingleL1_v"),
#        path1 = cms.vstring ("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07","hltOverlapFilterIsoMu27LooseChargedIsoPFTau20"),
#        path2 = cms.vstring ("hltPFTau20TrackLooseChargedIsoAgainstMuon","hltOverlapFilterIsoMu27LooseChargedIsoPFTau20"),
#        path3 = cms.vstring (""),
#        path4 = cms.vstring (""),
#        leg1 = cms.int32(13),
#        leg2 = cms.int32(15)
#        ),
#    cms.PSet (
#        HLT = cms.string("HLT_IsoMu27_MediumChargedIsoPFTau20_Trk1_eta2p1_SingleL1_v"),
#        path1 = cms.vstring ("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07","hltOverlapFilterIsoMu27MediumChargedIsoPFTau20"),
#        path2 = cms.vstring ("hltPFTau20TrackMediumChargedIsoAgainstMuon","hltOverlapFilterIsoMu27MediumChargedIsoPFTau20"),
#        path3 = cms.vstring (""),
#        path4 = cms.vstring (""),
#        leg1 = cms.int32(13),
#        leg2 = cms.int32(15)
#        ),
#    cms.PSet (
#        HLT = cms.string("HLT_IsoMu27_TightChargedIsoPFTau20_Trk1_eta2p1_SingleL1_v"),
#        path1 = cms.vstring ("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07","hltOverlapFilterIsoMu27TightChargedIsoPFTau20"),
#        path2 = cms.vstring ("hltPFTau20TrackTightChargedIsoAgainstMuon","hltOverlapFilterIsoMu27TightChargedIsoPFTau20"),
#        path3 = cms.vstring (""),
#        path4 = cms.vstring (""),
#        leg1 = cms.int32(13),
#        leg2 = cms.int32(15)
#        )


## === MET/AK8 -- UPDATED FOR 2018 DATA -- LP: for Resonant analysis tests.
    cms.PSet (
        HLT = cms.string("HLT_PFHT500_PFMET100_PFMHT100_IDTight_v"),
        path1 = cms.vstring ("",""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(999),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_AK8PFJet400_TrimMass30_v"),
        path1 = cms.vstring ("",""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(999),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v"),
        path1 = cms.vstring ("",""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(999),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v"), 
        path1 = cms.vstring ("",""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(999),
        leg2 = cms.int32(999)
        ),
     )

## === To study -- paths from valeria's studies

#HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4
#HLT_PFMET120_PFMHT120_IDTight
#HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5
#HLT_RsqMR300_Rsq0p09_MR200
#HLT_Mu50
#HLT_PFHT450_SixPFJet36_PFBTagDeepCSV_1p59
#HLT_Photon35_TwoProngs35
#HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1

#now I create the trigger list for HLTconfig
for i in range(len(HLTLIST)):
    tmpl =  str(HLTLIST[i].HLT).replace('cms.string(\'','') ##CMSSW Vaffanculo
    tmpl = tmpl.replace('\')','') ##CMSSW Vaffanculo x 2
    TRIGGERLIST.append(tmpl)
#print TRIGGERLIST
