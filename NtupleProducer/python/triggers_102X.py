import FWCore.ParameterSet.Config as cms

## TRIGGER VERSIONING: insert paths without the explicit version number, e.g. HLT_anything_v 
## do not remove any other part of path name or the check will give meaningless results!

## leg numbering: write the object type the two legs refer to
## 11: electrons, 13: muons, 15: taus, 999: others/unused (e.g. leg 2 in single muon)

TRIGGERLIST=[]
#list triggers and filter paths here!
# channel: kemu=0, ketau=1,kmutau=2,ktautau=3
HLTLIST = cms.VPSet(

### === Single-mu ===

    #cms.PSet (
    #    HLT = cms.string("HLT_IsoMu20_v"),
    #    path1 = cms.vstring (""),
    #    path2 = cms.vstring (""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(13),
    #    leg2 = cms.int32(999)
    #    ),
    #cms.PSet (
    #    HLT = cms.string("HLT_IsoTkMu20_v"),
    #    path1 = cms.vstring (""),
    #    path2 = cms.vstring (""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(13),
    #    leg2 = cms.int32(999)
    #    ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu22_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoTkMu22_v"), 
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu22_eta2p1_v"), 
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoTkMu22_eta2p1_v"), 
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu24_v"), 
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoTkMu24_v"), 
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(999)
        ),
    #cms.PSet (
    #    HLT = cms.string("HLT_IsoMu24_eta2p1_v"), 
    #    path1 = cms.vstring (""),
    #    path2 = cms.vstring (""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(13),
    #    leg2 = cms.int32(999)
    #    ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu27_v"), 
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(999)
        ),


### === Single-ele ===

    cms.PSet (
        HLT = cms.string("HLT_Ele25_eta2p1_WPTight_Gsf_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Ele27_WPTight_Gsf_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Ele27_eta2p1_WPLoose_Gsf_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(999)
        ),
    #cms.PSet (
    #    HLT = cms.string("HLT_Ele27_eta2p1_WPTight_Gsf_v"),
    #    path1 = cms.vstring (""),
    #    path2 = cms.vstring (""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(11),
    #    leg2 = cms.int32(999)
    #    ),
    #cms.PSet (
    #    HLT = cms.string("HLT_Ele28_WPTight_Gsf_v"),
    #    path1 = cms.vstring (""),
    #    path2 = cms.vstring (""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(11),
    #    leg2 = cms.int32(999)
    #    ),
    #cms.PSet (
    #    HLT = cms.string("HLT_Ele30_WPTight_Gsf_v"),
    #    path1 = cms.vstring (""),
    #    path2 = cms.vstring (""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(11),
    #    leg2 = cms.int32(999)
    #    ),
    cms.PSet (
        HLT = cms.string("HLT_Ele32_WPTight_Gsf_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Ele35_WPTight_Gsf_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(999)
        ),
    #cms.PSet (
    #    HLT = cms.string("HLT_Ele38_WPTight_Gsf_v"),
    #    path1 = cms.vstring (""),
    #    path2 = cms.vstring (""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(11),
    #    leg2 = cms.int32(999)
    #    ),
    #cms.PSet (
    #    HLT = cms.string("HLT_Ele40_WPTight_Gsf_v"),
    #    path1 = cms.vstring (""),
    #    path2 = cms.vstring (""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(11),
    #    leg2 = cms.int32(999)
    #    ),


### === mu+tauh ===

    cms.PSet (
        HLT = cms.string("HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1_v"), 
        path1 = cms.vstring ("",""),
        path2 = cms.vstring ("",""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v"), 
        path1 = cms.vstring ("",""),
        path2 = cms.vstring ("",""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
        ),  
    cms.PSet (
        HLT = cms.string("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_v"), 
        path1 = cms.vstring ("",""),
        path2 = cms.vstring ("",""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
        ), 
    #cms.PSet (
    #    HLT = cms.string("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_v"), 
    #    path1 = cms.vstring ("",""),
    #    path2 = cms.vstring ("",""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(13),
    #    leg2 = cms.int32(15)
    #    ),  
    #cms.PSet (
    #    HLT = cms.string("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_CrossL1_v"), 
    #    path1 = cms.vstring ("",""),
    #    path2 = cms.vstring ("",""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(13),
    #    leg2 = cms.int32(15)
    #    ), 
    #cms.PSet (
    #    HLT = cms.string("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_v"), 
    #    path1 = cms.vstring ("",""),
    #    path2 = cms.vstring ("",""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(13),
    #    leg2 = cms.int32(15)
    #    ), 
    #cms.PSet (
    #    HLT = cms.string("HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_CrossL1_v"), 
    #    path1 = cms.vstring ("",""),
    #    path2 = cms.vstring ("",""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(13),
    #    leg2 = cms.int32(15)
    #    ), 
    #cms.PSet (
    #    HLT = cms.string("HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_v"), 
    #    path1 = cms.vstring ("",""),
    #    path2 = cms.vstring ("",""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(13),
    #    leg2 = cms.int32(15)
    #    ), 
    #cms.PSet (
    #    HLT = cms.string("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1_v"), 
    #    path1 = cms.vstring ("",""),
    #    path2 = cms.vstring ("",""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(13),
    #    leg2 = cms.int32(15)
    #    ), 
    #cms.PSet (
    #    HLT = cms.string("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1_v"), 
    #    path1 = cms.vstring ("",""),
    #    path2 = cms.vstring ("",""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(13),
    #    leg2 = cms.int32(15)
    #    ), 
    #cms.PSet (
    #    HLT = cms.string("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_TightID_CrossL1_v"), 
    #    path1 = cms.vstring ("",""),
    #    path2 = cms.vstring ("",""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(13),
    #    leg2 = cms.int32(15)
    #    ),
    #cms.PSet (
    #    HLT = cms.string("HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1_v"), 
    #    path1 = cms.vstring ("",""),
    #    path2 = cms.vstring ("",""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(13),
    #    leg2 = cms.int32(15)
    #    ), 
    #cms.PSet (
    #    HLT = cms.string("HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_TightID_CrossL1_v"), 
    #    path1 = cms.vstring ("",""),
    #    path2 = cms.vstring ("",""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(13),
    #    leg2 = cms.int32(15)
    #    ),


### === ele+tauh ===

    cms.PSet (
        HLT = cms.string("HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v"), 
        path1 = cms.vstring ("",""),
        path2 = cms.vstring ("",""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v"), 
        path1 = cms.vstring ("",""),
        path2 = cms.vstring ("",""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau30_v"), 
        path1 = cms.vstring ("",""),
        path2 = cms.vstring ("",""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v"), 
        path1 = cms.vstring ("",""),
        path2 = cms.vstring ("",""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1_v"), 
        path1 = cms.vstring ("",""),
        path2 = cms.vstring ("",""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(15)
        ),
    #cms.PSet (
    #    HLT = cms.string("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_v"), 
    #    path1 = cms.vstring ("",""),
    #    path2 = cms.vstring ("",""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(11),
    #    leg2 = cms.int32(15)
    #    ),
    #cms.PSet (
    #    HLT = cms.string("HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_CrossL1_v"), 
    #    path1 = cms.vstring ("",""),
    #    path2 = cms.vstring ("",""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(11),
    #    leg2 = cms.int32(15)
    #    ),
    #cms.PSet (
    #    HLT = cms.string("HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_v"), 
    #    path1 = cms.vstring ("",""),
    #    path2 = cms.vstring ("",""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(11),
    #    leg2 = cms.int32(15)
    #    ),
    #cms.PSet (
    #    HLT = cms.string("HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_CrossL1_v"), 
    #    path1 = cms.vstring ("",""),
    #    path2 = cms.vstring ("",""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(11),
    #    leg2 = cms.int32(15)
    #    ),
    #cms.PSet (
    #    HLT = cms.string("HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_v"), 
    #    path1 = cms.vstring ("",""),
    #    path2 = cms.vstring ("",""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(11),
    #    leg2 = cms.int32(15)
    #    ),
    #cms.PSet (
    #    HLT = cms.string("HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1_v"), 
    #    path1 = cms.vstring ("",""),
    #    path2 = cms.vstring ("",""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(11),
    #    leg2 = cms.int32(15)
    #    ),
    #cms.PSet (
    #    HLT = cms.string("HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1_v"), 
    #    path1 = cms.vstring ("",""),
    #    path2 = cms.vstring ("",""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(11),
    #    leg2 = cms.int32(15)
    #    ),


### === tauh+tauh ===

    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("",""),
        path2 = cms.vstring ("",""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("",""),
        path2 = cms.vstring ("",""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("",""),
        path2 = cms.vstring ("",""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v"),
        path1 = cms.vstring ("",""),
        path2 = cms.vstring ("",""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v"),
        path1 = cms.vstring ("",""),
        path2 = cms.vstring ("",""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ), 
    cms.PSet (
        HLT = cms.string("HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("",""),
        path2 = cms.vstring ("",""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ), 
    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("",""),
        path2 = cms.vstring ("",""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
        ), 
    #cms.PSet (
    #    HLT = cms.string("HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_v"),
    #    path1 = cms.vstring ("",""),
    #    path2 = cms.vstring ("",""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(15),
    #    leg2 = cms.int32(15)
    #    ),
    #cms.PSet (
    #    HLT = cms.string("HLT_DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v"),
    #    path1 = cms.vstring ("",""),
    #    path2 = cms.vstring ("",""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(15),
    #    leg2 = cms.int32(15)
    #    ), 
    #cms.PSet (
    #    HLT = cms.string("HLT_DoubleTightChargedIsoPFTau35_Trk1_eta2p1_Reg_"),
    #    path1 = cms.vstring ("",""),
    #    path2 = cms.vstring ("",""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(15),
    #    leg2 = cms.int32(15)
    #    ),
    #cms.PSet (
    #    HLT = cms.string("HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1_Reg_v"),
    #    path1 = cms.vstring ("",""),
    #    path2 = cms.vstring ("",""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(15),
    #    leg2 = cms.int32(15)
    #    ),
    #cms.PSet (
    #    HLT = cms.string("HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg_v"),
    #    path1 = cms.vstring ("",""),
    #    path2 = cms.vstring ("",""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(15),
    #    leg2 = cms.int32(15)
    #    ), 
    #cms.PSet (
    #    HLT = cms.string("HLT_DoubleMediumChargedIsoPFTau40_Trk1_eta2p1_Reg_v"),
    #    path1 = cms.vstring ("",""),
    #    path2 = cms.vstring ("",""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(15),
    #    leg2 = cms.int32(15)
    #    ),
    #cms.PSet (
    #    HLT = cms.string("HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v"),
    #    path1 = cms.vstring ("",""),
    #    path2 = cms.vstring ("",""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(15),
    #    leg2 = cms.int32(15)
    #    ), 
    #cms.PSet (
    #    HLT = cms.string("HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_v"),
    #    path1 = cms.vstring ("",""),
    #    path2 = cms.vstring ("",""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(15),
    #    leg2 = cms.int32(15)
    #    ),
    #cms.PSet (
    #    HLT = cms.string("HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_eta2p1_Reg_v"),
    #    path1 = cms.vstring ("",""),
    #    path2 = cms.vstring ("",""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(15),
    #    leg2 = cms.int32(15)
    #    ), 
    #cms.PSet (
    #    HLT = cms.string("HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg_v"),
    #    path1 = cms.vstring ("",""),
    #    path2 = cms.vstring ("",""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(15),
    #    leg2 = cms.int32(15)
    #    ), 
    #cms.PSet (
    #    HLT = cms.string("HLT_DoubleTightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v"),
    #    path1 = cms.vstring ("",""),
    #    path2 = cms.vstring ("",""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(15),
    #    leg2 = cms.int32(15)
    #    ),   


### === mu+ele ===

    cms.PSet (  
        HLT = cms.string("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(11)
        ),
    cms.PSet ( 
        HLT = cms.string("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(11)
        ),
    cms.PSet ( 
        HLT = cms.string("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(11)
        ),
    cms.PSet (  
        HLT = cms.string("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(11)
        ),
    cms.PSet (  
        HLT = cms.string("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(11)
        ),
    cms.PSet ( 
        HLT = cms.string("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(11)
        ),
    cms.PSet ( 
        HLT = cms.string("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(11)
        ),
    cms.PSet ( 
        HLT = cms.string("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(11)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(11)
        ),
	
### === mu+mu ===

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
    cms.PSet (  
        HLT = cms.string("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(13)
        ),  
    cms.PSet (  
        HLT = cms.string("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(13)
        ),   
    #cms.PSet (
    #    HLT = cms.string("HLT_Mu30_TkMu11_v"),
    #    path1 = cms.vstring (""),
    #    path2 = cms.vstring (""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(13),
    #    leg2 = cms.int32(13)
    #    ),

### === ele+ele ===

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
        HLT = cms.string("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"),
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
    #cms.PSet (
    #    HLT = cms.string("HLT_DoubleEle25_CaloIdL_MW_v"),
    #    path1 = cms.vstring (""),
    #    path2 = cms.vstring (""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(11),
    #    leg2 = cms.int32(11)
    #    ),
    #cms.PSet (
    #    HLT = cms.string("HLT_DoubleEle33_CaloIdL_MW_v"),
    #    path1 = cms.vstring (""),
    #    path2 = cms.vstring (""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(11),
    #    leg2 = cms.int32(11)
    #    ),
    #cms.PSet (
    #    HLT = cms.string("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v"),
    #    path1 = cms.vstring (""),
    #    path2 = cms.vstring (""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(11),
    #    leg2 = cms.int32(11)
    #    ),
        
### === tri-lepton ===

    cms.PSet (  
        HLT = cms.string("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(11),
        leg3 = cms.int32(11)
        ),
    cms.PSet (  
        HLT = cms.string("HLT_TripleMu_12_10_5_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(13),
        leg3 = cms.int32(13)
        ),
    #cms.PSet (
    #    HLT = cms.string("HLT_TripleMu_10_5_5_DZ_v"),
    #    path1 = cms.vstring (""),
    #    path2 = cms.vstring (""),
    #    path3 = cms.vstring (""),
    #    path4 = cms.vstring (""),
    #    leg1 = cms.int32(13),
    #    leg2 = cms.int32(13),
    #    leg3 = cms.int32(13)
    #    ),
    cms.PSet (  
        HLT = cms.string("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(11),
        leg3 = cms.int32(11)
        ),
    cms.PSet (  
        HLT = cms.string("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(11),
        leg3 = cms.int32(11)
        ),
    cms.PSet (  
        HLT = cms.string("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(13),
        leg3 = cms.int32(11)
        ),
    cms.PSet (  
        HLT = cms.string("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v"),
        path1 = cms.vstring (""),
        path2 = cms.vstring (""),
        path3 = cms.vstring (""),
        path4 = cms.vstring (""),
        leg1 = cms.int32(13),
        leg2 = cms.int32(13),
        leg3 = cms.int32(11)
        ),

    )

#create the trigger list for HLTconfig
for i in range(len(HLTLIST)):
    tmpl =  str(HLTLIST[i].HLT).replace('cms.string(\'','') 
    tmpl = tmpl.replace('\')','') 
    TRIGGERLIST.append(tmpl)
#print TRIGGERLIST
