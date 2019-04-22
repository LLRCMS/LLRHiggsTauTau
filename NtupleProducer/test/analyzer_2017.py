#-----------------------------------------
#
#Producer controller
#
#-----------------------------------------
import os, re
PyFilePath = os.environ['CMSSW_BASE']+"/src/LLRHiggsTauTau/NtupleProducer/"

#YEAR = 2016
YEAR = 2017
#YEAR = 2018
print "YEAR: ",YEAR

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
USECLASSICSVFIT=True # if True use the ClassicSVfit package, if False use the SVFitStandAlone package

BUILDONLYOS=False #If true don't create the collection of SS candidates (and thus don't run SV fit on them)
APPLYTESCORRECTION=False # shift the central value of the tau energy scale before computing up/down variations
COMPUTEUPDOWNSVFIT=False # compute SVfit for up/down TES variation
COMPUTEMETUPDOWNSVFIT=False # compute SVfit for up/down MET JES variation
doCPVariables=False # compute CP variables and PV refit
COMPUTEQGVAR = True # compute QG Tagger for jets

IsMC=True
#IsMC=False
print "IsMC: ", IsMC
Is25ns=True

HLTProcessName='HLT' #Different names possible, check e.g. at https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD.
if not IsMC:
    HLTProcessName='HLT' #It always 'HLT' for real data
#print "HLTProcessName: ",HLTProcessName

#relaxed sets for testing purposes
TAUDISCRIMINATOR="byIsolationMVA3oldDMwoLTraw"
PVERTEXCUT="!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2" #cut on good primary vertexes
MUCUT="isLooseMuon && pt>5"
ELECUT="pt>5"#"pt>7"#"gsfTracsk.hitPattern().numberOfHits(HitPattern::MISSING_INNER_HITS)<=1 && pt>10"
TAUCUT="tauID('byCombinedIsolationDeltaBetaCorrRaw3Hits') < 1000.0 && pt>18"
JETCUT="pt>0" # was 10, is now 0 to save all the jets and be able to copute JEC MET in KLUB
LLCUT="mass>-99999"
BCUT="pt>5"

# ------------------------
DO_ENRICHED=False # do True by default, both ntuples and enriched outputs are saved!
STORE_ENRICHEMENT_ONLY=True # When True and DO_ENRICHED=True only collection additional to MiniAOD standard are stored. They can be used to reproduce ntuples when used together with oryginal MiniAOD with two-file-solution
# ------------------------

Is102X = True if 'CMSSW_10' in os.environ['CMSSW_VERSION'] else False
if Is102X: print "RELEASE: 102X"

Is92X = True if 'CMSSW_9' in os.environ['CMSSW_VERSION'] else False
if Is92X: print "RELEASE: 94X"

Is80X = True if 'CMSSW_8' in os.environ['CMSSW_VERSION'] else False
if Is80X: print "RELEASE: 80X"

##
## Standard sequence
##
if Is102X:
   execfile(PyFilePath+"python/HiggsTauTauProducer_102X.py")
elif Is92X:
   execfile(PyFilePath+"python/HiggsTauTauProducer_94X.py")
elif Is80X:
   execfile(PyFilePath+"python/HiggsTauTauProducer_80X.py")

### ----------------------------------------------------------------------
### Source, better to use sample to run on batch
### ----------------------------------------------------------------------
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    
    #2016 ttH
    #'/store/mc/RunIISummer16MiniAODv3/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/120000/F24F2D5E-DDEC-E811-AF50-90B11C08AD7D.root',
 
    #2017 ttH
    '/store/mc/RunIIFall17MiniAODv2/ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/100000/7C60AC2B-E76F-E811-9D60-0025905B860C.root',

    #2018 ttH 
    #'/store/mc/RunIIAutumn18MiniAOD/ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/270000/A912BBFA-D1A1-8544-A430-8C98C0767737.root',

    )
)

# process.source.skipEvents = cms.untracked.uint32(968)
#process.source.eventsToProcess = cms.untracked.VEventRange("1:2347130-1:2347130") # run only on event=2347130 (syntax= from run:evt - to run:evt)

#Limited nEv for testing purposes. -1 to run all events
process.maxEvents.input = 1000
#process.maxEvents.input = -1

# JSON mask for data --> defined in the lumiMask file
# from JSON file
#if not IsMC:
#  execfile(PyFilePath+"python/lumiMask.py")
#  process.source.lumisToProcess = LUMIMASK

##
## Output file
##
process.TFileService=cms.Service('TFileService',fileName=cms.string('HTauTauAnalysis_2017.root'))

# L1 trigger objects (as suggested on: https://twiki.cern.ch/twiki/bin/view/CMS/HiggsToTauTauWorking2017#Trigger_Information )
#  ----> TO BE FIXED <----
#process.out.outputCommands.append('keep *_caloStage2Digis_*_*') #FRA
#process.out.outputCommands.append('keep *_gmtStage2Digis_*_*')  #FRA 


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

#process.options = cms.PSet(skipEvent =  cms.untracked.vstring('ProductNotFound'))
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

#processDumpFile = open('process.dump' , 'w')
#print >> processDumpFile, process.dumpPython()
