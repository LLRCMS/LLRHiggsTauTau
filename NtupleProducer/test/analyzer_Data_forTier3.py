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

TRIGGERLIST = [#"HLT_*", #["HLT_Mu17_Mu8_v*", "HLT_Mu17_TkMu8_v*"] # to run on DATA/MC 2012 # "HLT_*" is a empty path
  "HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v1",
  "HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v1",
  "HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v1",
  "HLT_IsoMu24_eta2p1_IterTrk02_v1",
  "HLT_Ele22_eta2p1_WP85_Gsf_LooseIsoPFTau20_v1",
  "HLT_Ele27_eta2p1_WP85_Gsf_v1",
  "HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v1"
#    "HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v1",
#    "HLT_IsoMu17_eta2p1_v1",
#    "HLT_IsoMu17_eta2p1_MediumIsoPFTau40_Trk1_eta2p1_Reg_v1",
#    "HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1_v1",
#    "HLT_IsoMu24_eta2p1_IterTrk01_v1",
#    "HLT_IsoMu24_eta2p1_IterTrk02_v1",
#    "HLT_IsoMu24_eta2p1_IterTrk02_LooseIsoPFTau20_v1",
#    "HLT_Ele22_eta2p1_WP85_Gsf_LooseIsoPFTau20_v1",
#    "HLT_Ele32_eta2p1_WP85_Gsf_v1",
#    "HLT_Ele32_eta2p1_WP85_Gsf_LooseIsoPFTau20_v1",
#    "HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v1",
#    "HLT_IsoMu16_eta2p1_CaloMET30_LooseIsoPFTau50_Trk30_eta2p1_v1",
#    "HLT_IsoMu16_eta2p1_CaloMET30_v1",
#    "HLT_Mu16_eta2p1_CaloMET30_v1",
#    "HLT_LooseIsoPFTau50_Trk30_eta2p1_v1",
#    "HLT_DoubleIsoMu17_eta2p1_v1",
#    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v1",
#    "HLT_Ele27_eta2p1_WP85_Gsf_LooseIsoPFTau20_v1",
#    "HLT_Ele27_eta2p1_WP85_Gsf_v1"
]

#Samples:
IsMC=False

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

#from JSON file
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

