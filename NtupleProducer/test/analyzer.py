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

USEPAIRMET=False # input to SVfit: true: MVA pair MET; false: PFmet (HF inclusion set using USE_NOHFMET)
USE_NOHFMET = False # True to exclude HF and run on silver json

SVFITBYPASS=True # use SVFitBypass module, no SVfit computation, adds dummy userfloats for MET and SVfit mass
BUILDONLYOS=False #If true don't create the collection of SS candidates (and thus don't run SV fit on them)

IsMC=True
Is25ns=True

#relaxed sets for testing purposes
TAUDISCRIMINATOR="byIsolationMVA3oldDMwoLTraw"
PVERTEXCUT="!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2" #cut on good primary vertexes
MUCUT="isLooseMuon && pt>8"
ELECUT="pt>10"#"gsfTrack.hitPattern().numberOfHits(HitPattern::MISSING_INNER_HITS)<=1 && pt>10"
TAUCUT="tauID('byCombinedIsolationDeltaBetaCorrRaw3Hits') < 1000.0 && pt>18" #miniAOD tau from hpsPFTauProducer have pt>18 and decaymodefinding ID
JETCUT="pt>15"
LLCUT="mass>0"
BCUT="pt>5"


# ------------------------
DO_ENRICHED=False # do True by default, both ntuples and enriched outputs are saved!
# ------------------------

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
    '/store/mc/RunIIFall15MiniAODv2/SUSYGluGluToHToTauTau_M-160_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/50000/12184969-3DB8-E511-879B-001E67504A65.root',
    #'/store/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/70000/002ABFCA-A0B9-E511-B9BA-0CC47A57CD6A.root',
    #'/store/mc/RunIISpring15MiniAODv2/SUSYGluGluToBBHToTauTau_M-1000_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/60000/D4177994-FB78-E511-9D05-C4346BC76CD8.root',
    #'/store/mc/RunIISpring15MiniAODv2/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/003F1529-D36D-E511-9E33-001E6724816F.root',
    #'/store/mc/RunIISpring15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2_ext3-v1/10000/003964D7-D06E-E511-A8DA-001517F7F524.root', # miniAOD v2
    #'/store/data/Run2015D/Tau/MINIAOD/PromptReco-v4/000/259/861/00000/A2D00E3D-457E-E511-BDEF-02163E012584.root', # DATA 2015D
    #'/store/data/Run2015D/Tau/MINIAOD/PromptReco-v4/000/258/159/00000/06064263-DD6B-E511-99BA-02163E013861.root',
    #  '/store/data/Run2015D/Tau/MINIAOD/05Oct2015-v1/30000/00756739-616F-E511-BC4F-003048FFCBB0.root',
    #  '/store/data/Run2015D/SingleMuon/MINIAOD/PromptReco-v4/000/258/159/00000/6CA1C627-246C-E511-8A6A-02163E014147.root',
    #  '/store/data/Run2015D/SingleMuon/MINIAOD/05Oct2015-v1/10000/021FD3F0-876F-E511-99D2-0025905A6060.root',
    #  '/store/data/Run2015D/Tau/MINIAOD/05Oct2015-v1/30000/00756739-616F-E511-BC4F-003048FFCBB0.root',   
    #  '/store/data/Run2015D/SingleElectron/MINIAOD/PromptReco-v4/000/258/159/00000/D0E6617B-186C-E511-8FEE-02163E01452F.root',
    #  '/store/data/Run2015C/SingleMuon/MINIAOD/PromptReco-v1/000/253/944/00000/BA286EB7-5141-E511-8829-02163E01437D.root',
    #  '/store/data/Run2015C/SingleMuon/MINIAOD/PromptReco-v1/000/253/888/00000/DA473838-0941-E511-9644-02163E014405.root',
    #  '/store/mc/RunIISpring15DR74/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/048FB1EE-33FD-E411-A2BA-0025905A6094.root',
    #  'file:0E885FA4-FF01-E511-B2DB-002590A3C95E.root',
    #  '/store/mc/RunIISpring15DR74/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/048FB1EE-33FD-E411-A2BA-0025905A6094.root',
    #  '/store/mc/RunIISpring15DR74/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v2/00000/06B5178E-F008-E511-A2CF-00261894390B.root'
    #  'file:/data_CMS/cms/salerno/Lambdam4_74x_step2/miniAOD_lambdam4_2_300000Events_0Skipped_1435089918.88/output_65.root',
    #  'file:/data_CMS/cms/salerno/Lambda20_74x_step2/miniAOD_lambda20_2_300000Events_0Skipped_1434450687.86/output_0.root',
    #  '/store/mc/RunIISpring15DR74/TTTo2L2Nu_13TeV-powheg/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/0457C8B3-BA09-E511-80D0-008CFA56D58C.root',
    #  'file:/data_CMS/cms/ortona/Lambda20_step3/miniAOD_lambda20_3_300000Events_0Skipped_1425913946.36/output_0.root',
    #  '/store/mc/Phys14DR/GluGluToHToTauTau_M-125_13TeV-powheg-pythia6/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/1813E94A-D36E-E411-8EDC-3417EBE34D08.root'
    )
)

#Limited nEv for testing purposes. -1 to run all events
process.maxEvents.input = 100

# JSON mask for data --> defined in the lumiMask file
# from JSON file
#if not IsMC:
#  execfile(PyFilePath+"python/lumiMask.py")
#  process.source.lumisToProcess = LUMIMASK

##
## Output file
##

process.TFileService=cms.Service('TFileService',fileName=cms.string('HTauTauAnalysis.root'))

if DO_ENRICHED:
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
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
#process.MessageLogger.categories.append('onlyError')
#process.MessageLogger.cerr.onlyError=cms.untracked.PSet(threshold  = cms.untracked.string('ERROR'))
#process.MessageLogger.cerr.threshold='ERROR'
#process.MessageLogger = cms.Service("MessageLogger",
#	destinations = cms.untracked.vstring('log.txt')
#)
#process.MessageLogger.threshold = cms.untracked.string('ERROR')

