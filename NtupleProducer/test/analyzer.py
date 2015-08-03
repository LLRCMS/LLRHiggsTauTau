#-----------------------------------------
#
#Producer controller
#
#-----------------------------------------
import os, re
PyFilePath = os.environ['CMSSW_BASE']+"/src/LLRHiggsTauTau/NtupleProducer/"

#samples list (it could be moved to a cfg file for better reading
#samples = [
#]
#apply corrections?
APPLYMUCORR=False
APPLYELECORR=True
APPLYFSR=False #this is by far the slowest module (not counting SVFit so far)
#Cuts on the Objects (add more cuts with &&)
#MUCUT="(isGlobalMuon || (isTrackerMuon && numberOfMatches>0)) && abs(eta)<2.4 && pt>8"
#ELECUT="abs(eta)<2.5 && gsfTrack.trackerExpectedHitsInner.numberOfHits<=1 && pt>10"
#TAUCUT="pt>15"
#JETCUT="pt>15"

USEPAIRMET=False
SVFITBYPASS=True # use SVFitBypass module, no SVfit computation, adds dummy userfloats for MET and SVfit mass
RUN_NTUPLIZER=True
BUILDONLYOS=False #If true don't create the collection of SS candidates (and thus don't run SV fit on them)

#relaxed sets for testing purposes
TAUDISCRIMINATOR="byIsolationMVA3oldDMwoLTraw"
PVERTEXCUT="!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2" #cut on good primary vertexes
MUCUT="(isGlobalMuon || (isTrackerMuon && numberOfMatches>0)) && pt>8"
ELECUT="userFloat('missingHit')<=1 && pt>10"#"gsfTrack.hitPattern().numberOfHits(HitPattern::MISSING_INNER_HITS)<=1 && pt>10"
TAUCUT="pt>18" #miniAOD tau from hpsPFTauProducer have pt>18 and decaymodefinding ID
JETCUT="pt>15"
LLCUT="mass>0"
BCUT="pt>5"

IsMC=True

# spring15 triggers
# NOTE: wildcards for the moment are only allowed for trigger filter but not to set trigger bits in ntuples
# NOTE2 (FIXME) : no trigger matching with reco objects for the moment!

#TRIGGERLIST = cms.vstring (
#      ### e mu triggers
#  "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1",
#  "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v1",
#  "HLT_IsoMu24_eta2p1_v1", # NOT IN BASELINE
#  "HLT_IsoMu27_v1", # NOT IN BASELINE
#  
#      ### mu tau triggers
#  "HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v1",
#  "HLT_IsoMu24_eta2p1_v1",
#  "HLT_IsoMu27_v1", # NOT IN BASELINE
#  
#      ### e tau triggers
#  "HLT_Ele22_eta2p1_WP75_Gsf_LooseIsoPFTau20_v1",
#  "HLT_Ele32_eta2p1_WP75_Gsf_v1",
#
#      ### tau tau triggers
#  "HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v1",
#)
#
##DATA samples
#if not  IsMC:
#    TRIGGERLIST = cms.vstring (
#        ### e mu triggers
#        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2",
#        "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v2",
#
#        ### mu tau
#        "HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v2",
#        "HLT_IsoMu24_eta2p1_IterTrk02_v2",
#
#        ### e tau
#        "HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v1",
#        "HLT_Ele32_eta2p1_WPTight_Gsf_v1",
#
#        ### tau tau
#        "HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v2"	
#        )

#word = cms.string("vaffanculo")
#print word.c_str()
#print word.data()
#print word.Data()
#print str(word)

#TRIGGERLIST=cms.vstring()
TRIGGERLIST=[]
#list triggers and filter paths here!
# channel: kemu=0, ketau=1,kmutau=2,ktautau=3
HLTLIST = cms.VPSet(

### ----------------- SPRING 15 - MC
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

### ----------------- DATA, August 2015
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


##
## Standard sequence
##
execfile(PyFilePath+"python/HiggsTauTauProducer.py")

### ----------------------------------------------------------------------
### Source, better to use sample to run on batch
### ----------------------------------------------------------------------
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      '/store/mc/RunIISpring15DR74/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v2/00000/06B5178E-F008-E511-A2CF-00261894390B.root'
    #  'file:/data_CMS/cms/salerno/Lambdam4_74x_step2/miniAOD_lambdam4_2_300000Events_0Skipped_1435089918.88/output_65.root',
    #	 'file:/data_CMS/cms/salerno/Lambda20_74x_step2/miniAOD_lambda20_2_300000Events_0Skipped_1434450687.86/output_0.root',
    #	 '/store/mc/RunIISpring15DR74/TTTo2L2Nu_13TeV-powheg/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/0457C8B3-BA09-E511-80D0-008CFA56D58C.root',
    #  'file:/data_CMS/cms/ortona/Lambda20_step3/miniAOD_lambda20_3_300000Events_0Skipped_1425913946.36/output_0.root',
    #  '/store/mc/Phys14DR/GluGluToHToTauTau_M-125_13TeV-powheg-pythia6/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/1813E94A-D36E-E411-8EDC-3417EBE34D08.root'
    )
)

#Limited nEv for testing purposes. -1 to run all events
process.maxEvents.input = 1000

##
## Output file
##

if RUN_NTUPLIZER:
    process.TFileService=cms.Service('TFileService',fileName=cms.string('HTauTauAnalysis.root'))
else:
    process.out = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string('Enriched_miniAOD.root'),
        #outputCommands = cms.untracked.vstring(['keep *']),
        fastCloning     = cms.untracked.bool(False)
    )
    process.end = cms.EndPath(process.out)

#process.options = cms.PSet(skipEvent =  cms.untracked.vstring('ProductNotFound')),
#process.p = cms.EndPath(process.HTauTauTree)
process.p = cms.Path(process.Candidates)

# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1
#process.MessageLogger.categories.append('onlyError')
#process.MessageLogger.cerr.onlyError=cms.untracked.PSet(threshold  = cms.untracked.string('ERROR'))
#process.MessageLogger.cerr.threshold='ERROR'
#process.MessageLogger = cms.Service("MessageLogger",
#	destinations = cms.untracked.vstring('log.txt')
#)
#process.MessageLogger.threshold = cms.untracked.string('ERROR')

