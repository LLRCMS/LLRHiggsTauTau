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

#Samples:
IsMC=XXX_ISMC_XXX
Is25ns=True

#relaxed sets for testing purposes
TAUDISCRIMINATOR="byIsolationMVA3oldDMwoLTraw"
PVERTEXCUT="!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2" #cut on good primary vertexes
MUCUT="isLooseMuon && pt>8"
ELECUT="userFloat('missingHit')<=1 && pt>10"#"gsfTrack.hitPattern().numberOfHits(HitPattern::MISSING_INNER_HITS)<=1 && pt>10"
TAUCUT="tauID('byCombinedIsolationDeltaBetaCorrRaw3Hits') < 1000.0 && pt>18" #miniAOD tau from hpsPFTauProducer have pt>18 and decaymodefinding ID
JETCUT="pt>15"
LLCUT="mass>0"
BCUT="pt>5"


##
## Standard sequence
##
execfile(PyFilePath+"python/triggers.py") # contains the list of triggers and filters
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

# JSON mask for data --> defined in the lumiMask file
# from JSON file
if not IsMC:
  execfile(PyFilePath+"python/lumiMask.py")
  process.source.lumisToProcess = LUMIMASK

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

