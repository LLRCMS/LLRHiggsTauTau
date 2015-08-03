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

##
## Standard sequence
##
execfile(PyFilePath+"python/triggers.py") # contains the list of triggers and filters
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

