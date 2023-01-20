#-----------------------------------------
#
#Producer controller
#
#-----------------------------------------
import os, re
PyFilePath = os.environ['CMSSW_BASE']+"/src/LLRHiggsTauTau/NtupleProducer/"

# Year/Period
YEAR   = 2018
#YEAR   = 2016
PERIOD = ''
#PERIOD = 'postVFP' # use 'postVFP' for 2016 VFP samples, can be left empty if running on 2017 and 2018
IsMC=False

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
process.source = cms.Source("PoolSource",
    # Files listed here used for synch with Pisa
    fileNames = cms.untracked.vstring(
        #EGamma dataset
#        '/store/data/Run2018A/EGamma/MINIAOD/UL2018_MiniAODv2_GT36-v1/2820000/015BEACB-338C-894D-8EB5-B5AB2A7B8E81.root'
        #SingleMuon dataset
#        '/store/data/Run2018B/SingleMuon/MINIAOD/UL2018_MiniAODv2_GT36-v1/2520000/0139F1B3-834C-334E-ABA7-08828EDA119B.root'
        #Tau dataset
#        '/store/data/Run2018C/Tau/MINIAOD/UL2018_MiniAODv2_GT36-v1/2830000/03143451-FB26-F24A-B542-4165F66D5416.root'
        #MET dataset#
#        '/store/data/Run2018D/MET/MINIAOD/UL2018_MiniAODv2_GT36-v1/2820000/10089943-A30E-6446-962D-93332C37EE71.root'
        #DY NLO
        '/store/mc/RunIISummer20UL18MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/230000/0191D778-A59B-5149-9EA2-1FF39D787429.root'
        #Signal (ggf spin 0, mX=850)
#        '/store/mc/RunIISummer20UL18MiniAODv2/GluGluToRadionToHHTo2B2Tau_M-850_TuneCP5_PSWeights_narrow_13TeV-madgraph-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/30000/5CAFEC30-CAB2-4741-9E25-6931416D3F59.root'
    )
)

# process.source.skipEvents = cms.untracked.uint32(968)
#process.source.eventsToProcess = cms.untracked.VEventRange("1:2347130-1:2347130") # run only on event=2347130 (syntax= from run:evt - to run:evt)

#Limited nEv for testing purposes. -1 to run all events
process.maxEvents.input = -1 #10000

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
process.MessageLogger.cerr.FwkReport.reportEvery = 5000
#process.MessageLogger.categories.append('onlyError')
#process.MessageLogger.cerr.onlyError=cms.untracked.PSet(threshold  = cms.untracked.string('ERROR'))
#process.MessageLogger.cerr.threshold='ERROR'
#process.MessageLogger = cms.Service("MessageLogger",
#	destinations = cms.untracked.vstring('log.txt')
#)
#process.MessageLogger.threshold = cms.untracked.string('ERROR')

#processDumpFile = open('process.dump' , 'w')
#print >> processDumpFile, process.dumpPython()
