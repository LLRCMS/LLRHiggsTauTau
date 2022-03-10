# coding: utf-8

import os
import re

from FWCore.ParameterSet.VarParsing import VarParsing


# setup options
options = VarParsing("python")

# set defaults of common options
# options.setDefault("inputFiles", "")  # currently, a list is used to toggle the input for testing
options.setDefault("outputFile", "HTauTauAnalysis.root")
options.setDefault("maxEvents", -1)

# register new options
options.register("year", 2018, VarParsing.multiplicity.singleton, VarParsing.varType.int,
    "year of data-taking or MC campaign")
options.register("mc", True, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
    "whether the input dataset is from MC simulation or not")
options.register("reportEvery", 5000, VarParsing.multiplicity.singleton, VarParsing.varType.int,
    "number of events after which a report message is written")
options.register("wantSummary", True, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
    "print a summary at the end")

# get arguments of sys.argv
options.parseArguments()


#-----------------------------------------
#
#Producer controller
#
#-----------------------------------------

PyFilePath = os.environ['CMSSW_BASE']+"/src/LLRHiggsTauTau/NtupleProducer/"

# Year/Period
YEAR   = options.year
PERIOD = '' # use 'postVFP' for 2016 VFP samples, can be left empty if running on 2017 and 2018

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

USEPAIRMET=False # input to SVfit: true: MVA pair MET; false: PFmet (HF inclusion set using USE_NOHFMET)
APPLYMETCORR=False # flag to enable (True) and disable (False) Z-recoil corrections for MVA MET response and resolution
USE_NOHFMET = False # True to exclude HF and run on silver json


SVFITBYPASS=True # use SVFitBypass module, no SVfit computation, adds dummy userfloats for MET and SVfit mass
#USECLASSICSVFIT=True # if True use the ClassicSVfit package, if False use the SVFitStandAlone package

BUILDONLYOS=False #If true don't create the collection of SS candidates (and thus don't run SV fit on them)
APPLYTESCORRECTION=True # shift the central value of the tau energy scale before computing up/down variations
COMPUTEUPDOWNSVFIT=False # compute SVfit for up/down TES variation
COMPUTEMETUPDOWNSVFIT=False # compute SVfit for up/down MET JES variation
doCPVariables=False # compute CP variables and PV refit
COMPUTEQGVAR = False # compute QG Tagger for jets
IsMC = options.mc
Is25ns=True
HLTProcessName='HLT' #Different names possible, check e.g. at https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD.
if not IsMC:
    HLTProcessName='HLT' #It always 'HLT' for real data
print "HLTProcessName: ",HLTProcessName

#relaxed sets for testing purposes
TAUDISCRIMINATOR="byIsolationMVA3oldDMwoLTraw"
PVERTEXCUT="!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2" #cut on good primary vertexes
MUCUT="isLooseMuon && pt>10"#"isLooseMuon && pt>5"
ELECUT="pt>10"#"pt>7"#"gsfTracsk.hitPattern().numberOfHits(HitPattern::MISSING_INNER_HITS)<=1 && pt>10"
TAUCUT="pt>15"#"tauID('byCombinedIsolationDeltaBetaCorrRaw3Hits') < 1000.0 && pt>18" #miniAOD tau from hpsPFTauProducer have pt>18 and decaymodefinding ID
JETCUT="pt>0" # was 10, is now 0 to save all the jets and be able to copute JEC MET in KLUB
LLCUT="mass>-99999"
BCUT="pt>5"

# ------------------------
DO_ENRICHED=False # do True by default, both ntuples and enriched outputs are saved!
STORE_ENRICHEMENT_ONLY=True # When True and DO_ENRICHED=True only collection additional to MiniAOD standard are stored. They can be used to reproduce ntuples when used together with oryginal MiniAOD with two-file-solution
# ------------------------

is106X = True if 'CMSSW_10_6' in os.environ['CMSSW_VERSION'] else False
print "is106X:" , is106X
is102X = True if 'CMSSW_10_2' in os.environ['CMSSW_VERSION'] else False
print "is102X:" , is102X
is94X = True if 'CMSSW_9' in os.environ['CMSSW_VERSION'] else False# True to run in 92X (2017), False to run in 80X (2016) or 76X (2015)
print "is94X: " , is94X
is80X = True if 'CMSSW_8' in os.environ['CMSSW_VERSION'] else False# True to run in 80X (2016), False to run in 76X (2015)
print "is80X: " , is80X

##
## Standard sequence
##

if is106X:
    execfile(PyFilePath+"python/HiggsTauTauProducer_106X.py")
elif is102X:
    execfile(PyFilePath+"python/HiggsTauTauProducer_102X.py")
elif is94X:
    execfile(PyFilePath+"python/HiggsTauTauProducer_94X.py")
elif is80X:
    execfile(PyFilePath+"python/HiggsTauTauProducer_80X.py")
else :
    execfile(PyFilePath+"python/HiggsTauTauProducer.py")

### ----------------------------------------------------------------------
### Source, better to use sample to run on batch
### ----------------------------------------------------------------------

# example files mapped to (year, mc)
example_files = {
    # sync files, legacy 2016
    # (2016, False): "/store/data/Run2016H/SingleElectron/MINIAOD/17Jul2018-v1/00000/0026BF69-A58A-E811-BBA8-1866DA87AB31.root",
    # (2016, False): "/store/data/Run2016F/SingleMuon/MINIAOD/17Jul2018-v1/00000/002F631B-E98D-E811-A8FF-1CB72C1B64E6.root",
    # (2016, False): "/store/data/Run2016C/Tau/MINIAOD/17Jul2018-v1/40000/FC60B37F-468A-E811-8C7B-0CC47A7C3424.root",
    # (2016, True): "/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2/270000/EC774BCF-16E9-E811-B699-AC1F6B0DE2F4.root",
    # (2016, True): "/store/mc/RunIISummer16MiniAODv3/VBFHHTo2B2Tau_CV_1_C2V_1_C3_1_13TeV-madgraph/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/70000/E8B2E181-FF3E-E911-B72D-0025904C7FC2.root",

    # sync files, legacy 2017
    # (2017, False): "/store/data/Run2017D/SingleElectron/MINIAOD/31Mar2018-v1/90000/FEFBCFEA-5939-E811-99DA-0025905B85D6.root",
    # (2017, False): "/store/data/Run2017F/SingleMuon/MINIAOD/31Mar2018-v1/80004/EEC206A7-5737-E811-8E0A-0CC47A2B0388.root",
    # (2017, False): "/store/data/Run2017C/Tau/MINIAOD/31Mar2018-v1/90000/FE83BD44-CD37-E811-8731-842B2B180922.root",
    # (2017, True): "/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/20000/0056F0F4-4F44-E811-B415-FA163E5B5253.root",
    # (2017, True): "/store/mc/RunIIFall17MiniAODv2/VBFHHTo2B2Tau_CV_1_C2V_1_C3_1_13TeV-madgraph/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/90000/A6C340CF-9588-E811-9E52-782BCB46E733.root",
    # (2017, True): "/store/mc/RunIIFall17MiniAODv2/VBFHHTo2B2Tau_CV_1_C2V_1_C3_1_13TeV-madgraph/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/90000/941B9AB4-A189-E811-BA1C-0242AC130002.root",

    # sync files, legacy 2018
    # (2018, False): "/store/data/Run2018D/EGamma/MINIAOD/22Jan2019-v2/110000/2FD9167C-A197-E749-908D-3809F7B83A23.root",
    # (2018, False): "/store/data/Run2018D/SingleMuon/MINIAOD/22Jan2019-v2/110000/0011A0F3-23A3-D444-AD4D-2A2FFAE23796.root",
    # (2018, False): "/store/data/Run2018D/Tau/MINIAOD/PromptReco-v2/000/320/497/00000/584D46DF-5F95-E811-9646-FA163E293146.root",
    # (2018, True): "/store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/B3F93EA2-04C6-E04E-96AF-CB8FAF67E6BA.root",
    # (2018, True): "/store/mc/RunIIAutumn18MiniAOD/VBFHHTo2B2Tau_CV_1_C2V_1_C3_1_TuneCP5_PSWeights_13TeV-madgraph-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/250000/F22484D3-E820-F040-8E63-0864695D025B.root",

    # ultra-legacy test files
    (2016, True): "/store/mc/RunIISummer20UL16MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_v17-v1/2540000/127921F5-CAC9-AB4B-8C99-62981DB57E45.root",
    (2017, True): "/store/mc/RunIISummer20UL17MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_mc2017_realistic_v9-v2/120000/05668F8F-6DE3-3649-B228-9D7620F7C279.root",
    (2018, True): "/store/mc/RunIISummer20UL18MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/120000/001C8DDF-599C-5E45-BF2C-76F887C9ADE9.root",
}

process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(example_files.get((options.year, options.mc), "")))

# process.source.skipEvents = cms.untracked.uint32(968)
#process.source.eventsToProcess = cms.untracked.VEventRange("1:2347130-1:2347130") # run only on event=2347130 (syntax= from run:evt - to run:evt)

#Limited nEv for testing purposes. -1 to run all events
process.maxEvents.input = options.maxEvents

# overwrite some process options
process.options.wantSummary = cms.untracked.bool(options.wantSummary)

# JSON mask for data --> defined in the lumiMask file
# from JSON file
if not IsMC:
  if YEAR == 2016:
    execfile(PyFilePath+"python/lumiMask_2016.py")
  if YEAR == 2017:
    execfile(PyFilePath+"python/lumiMask_2017.py")
  if YEAR == 2018:
    execfile(PyFilePath+"python/lumiMask_2018.py")
  process.source.lumisToProcess = LUMIMASK

##
## Output file
##
process.TFileService=cms.Service('TFileService',fileName=cms.string('HTauTauAnalysis.root'))
# process.TFileService=cms.Service('TFileService', fileName=cms.string(options.__getattr__("outputFile", noTags=True)))
#process.TFileService=cms.Service('TFileService',fileName=cms.string('refFiles/Mu16_sync.root'))

if DO_ENRICHED:
    process.out = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string('Enriched_miniAOD.root'),
        outputCommands = cms.untracked.vstring('keep *'),
        fastCloning     = cms.untracked.bool(False),
        #Compression settings from MiniAOD allowing to save about 10% of disc space compared to defults ->
        compressionAlgorithm = cms.untracked.string('LZMA'),
        compressionLevel = cms.untracked.int32(4),
        dropMetaData = cms.untracked.string('ALL'),
        eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
        overrideInputFileSplitLevels = cms.untracked.bool(True)
        # <-
    )
    if STORE_ENRICHEMENT_ONLY:
        # Store only additional collections compared to MiniAOD necessary to reproduce ntuples (basically MVAMET, lepton pairs with SVFit and corrected jets)
        # Size of about 10% of full EnrichedMiniAOD
        process.out.outputCommands.append('drop *')
        process.out.outputCommands.append('keep *_SVllCand_*_*')
        process.out.outputCommands.append('keep *_SVbypass_*_*')
        process.out.outputCommands.append('keep *_barellCand_*_*')
        process.out.outputCommands.append('keep *_corrMVAMET_*_*')
        process.out.outputCommands.append('keep *_MVAMET_*_*')
        process.out.outputCommands.append('keep *_jets_*_*')
        process.out.outputCommands.append('keep *_patJetsReapplyJEC_*_*')
        process.out.outputCommands.append('keep *_softLeptons_*_*')
        process.out.outputCommands.append('keep *_genInfo_*_*')
        #process.out.fileName = 'EnrichementForMiniAOD.root' #FIXME: change name of output file?
    process.end = cms.EndPath(process.out)


process.p = cms.Path(process.Candidates)

# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
#process.MessageLogger.categories.append('onlyError')
#process.MessageLogger.cerr.onlyError=cms.untracked.PSet(threshold  = cms.untracked.string('ERROR'))
#process.MessageLogger.cerr.threshold='ERROR'
#process.MessageLogger = cms.Service("MessageLogger",
#	destinations = cms.untracked.vstring('log.txt')
#)
#process.MessageLogger.threshold = cms.untracked.string('ERROR')

#processDumpFile = open('process.dump' , 'w')
#print >> processDumpFile, process.dumpPython()
