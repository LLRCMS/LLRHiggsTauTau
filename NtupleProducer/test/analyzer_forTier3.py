#-----------------------------------------
#
#Producer controller
#
#-----------------------------------------
import os
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

#Samples:
IsMC=XXX_ISMC_XXX

##
## Standard sequence
##

execfile(PyFilePath+"python/HiggsTauTauProducer.py")
#execfile(PyFilePath+"python/DoubleHiggsProducer.py") #doesn't have mu and e --> it is just hadronic tau final states...

### ----------------------------------------------------------------------
### Source, better to use sample to run on batch
### ----------------------------------------------------------------------

#extract the file list from an external file (defines the FILELIST variable)
execfile(PyFilePath+"test/XXX_SAMPLEFILENAME_XXX")

process.source = cms.Source("PoolSource",
    fileNames = FILELIST
    )

#Limited nEv for testing purposes. -1 to run all events
process.maxEvents.input = cms.untracked.int32 (XXX_MAXEVENTS_XXX)
process.source.skipEvents = cms.untracked.uint32 (XXX_SKIPEVENTS_XXX)

# JSON mask for data --> to upload each time with the proper list using the macro in tools/
#from JSON file
if not IsMC:
  process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange( *(
      '251244:85-251244:86',
      '251244:88-251244:93',
      '251244:96-251244:121',
      '251244:123-251244:156',
      '251244:158-251244:428',
      '251244:430-251244:442',
      '251251:1-251251:31',
      '251251:33-251251:97',
      '251251:99-251251:167',
      '251252:1-251252:283',
      '251252:285-251252:505',
      '251252:507-251252:554',
      '251561:1-251561:94',
      '251562:1-251562:439',
      '251562:443-251562:691',
      '251643:1-251643:216',
      '251643:222-251643:606',
      '251721:21-251721:36',
      '251721:123-251721:244',
      '251883:56-251883:56',
      '251883:58-251883:60',
      '251883:62-251883:144',
      '251883:156-251883:437',
  ))

##
## Output file
##
if RUN_NTUPLIZER:
    process.TFileService=cms.Service('TFileService',fileName=cms.string("XXX_OUTPUTFILE_XXX"))
else:
    process.out = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string("XXX_OUTPUTFILE_XXX"),
        #outputCommands = cms.untracked.vstring(['keep *']),
        fastCloning     = cms.untracked.bool(False)
    )
    process.end = cms.EndPath(process.out)

#process.options = cms.PSet(skipEvent =  cms.untracked.vstring('ProductNotFound')),
#process.p = cms.EndPath(process.HTauTauTree)
process.p = cms.Path(process.Candidates)

# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
#process.MessageLogger.categories.append('onlyError')
#process.MessageLogger.cerr.onlyError=cms.untracked.PSet(threshold  = cms.untracked.string('ERROR'))
#process.MessageLogger.cerr.threshold='ERROR'
#process.MessageLogger = cms.Service("MessageLogger",
#	destinations = cms.untracked.vstring('log.txt')
#)
#process.MessageLogger.threshold = cms.untracked.string('ERROR')

