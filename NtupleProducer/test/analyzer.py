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
IsMC=True

##
## Standard sequence
##
execfile(PyFilePath+"python/HiggsTauTauProducer.py")

### ----------------------------------------------------------------------
### Source, better to use sample to run on batch
### ----------------------------------------------------------------------
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    'file:00C90EFC-3074-E411-A845-002590DB9262.root'
#         'file:/data_CMS/cms/ortona/Lambda20_step3/miniAOD_lambda20_3_300000Events_0Skipped_1425913946.36/output_0.root',
#        '/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00C90EFC-3074-E411-A845-002590DB9262.root',
#        'root://cms-xrd-global.cern.ch//store/mc/Phys14DR/GluGluToHToTauTau_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/2405749F-8B6F-E411-88EE-848F69FD2910.root',
#        'root://cms-xrd-global.cern.ch//store/mc/Phys14DR/GluGluToHToTauTau_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/3C05111C-8C6F-E411-A93E-7845C4F91450.root',
#        'root://cms-xrd-global.cern.ch//store/mc/Phys14DR/GluGluToHToTauTau_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/4851965D-EF6F-E411-9C94-3417EBE338FA.root',
#        'root://cms-xrd-global.cern.ch//store/mc/Phys14DR/GluGluToHToTauTau_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/4CC32E9F-8B6F-E411-9E56-00266CF9B4C0.root',
#        'root://cms-xrd-global.cern.ch//store/mc/Phys14DR/GluGluToHToTauTau_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/7E046B54-EF6F-E411-81DC-7845C4FC35D2.root',
#        'root://cms-xrd-global.cern.ch//store/mc/Phys14DR/GluGluToHToTauTau_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/186C2230-BA6F-E411-9FBD-7845C4FC365C.root',
#        'root://cms-xrd-global.cern.ch//store/mc/Phys14DR/GluGluToHToTauTau_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/3E56199D-1D70-E411-B425-7845C4FC36CB.root',
#        'root://cms-xrd-global.cern.ch//store/mc/Phys14DR/GluGluToHToTauTau_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/4A64F4C5-B96F-E411-B470-00266CF25490.root',
#        'root://cms-xrd-global.cern.ch//store/mc/Phys14DR/GluGluToHToTauTau_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/646AB792-1D70-E411-A935-00266CF9C0F0.root',
#        'root://cms-xrd-global.cern.ch//store/mc/Phys14DR/GluGluToHToTauTau_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/86CFA7C5-B96F-E411-B077-00266CF25490.root',
#        'root://cms-xrd-global.cern.ch//store/mc/Phys14DR/GluGluToHToTauTau_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/B6ACCFEC-C66F-E411-A6DA-7845C4FC378E.root',
#        'root://cms-xrd-global.cern.ch//store/mc/Phys14DR/GluGluToHToTauTau_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/D293773F-BA6F-E411-950E-00A0D1EE25D0.root'
        #"file:"+PyFilePath+'/test/2405749F-8B6F-E411-88EE-848F69FD2910.root' #local copy
        #'/store/mc/Phys14DR/GluGluToHToTauTau_M-125_13TeV-powheg-pythia6/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/1813E94A-D36E-E411-8EDC-3417EBE34D08.root'
        #"/store/mc/Spring14miniaod//GluGluToHToTauTau_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/1062D5FF-2D09-E411-943C-0025900EB52A.root",
        #"/store/mc/Spring14miniaod//GluGluToHToTauTau_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/2640D54D-2D09-E411-9FAA-003048D47670.root",
        #"/store/mc/Spring14miniaod//GluGluToHToTauTau_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/304D2104-2D09-E411-9BBC-0025900EB52A.root",
        #"/store/mc/Spring14miniaod//GluGluToHToTauTau_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/540AB9B2-2D09-E411-B413-001517357DDE.root",
        #"/store/mc/Spring14miniaod//GluGluToHToTauTau_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/6A959CA1-2D09-E411-8B84-0025900EB1A0.root",
        #"/store/mc/Spring14miniaod//GluGluToHToTauTau_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/78F2C486-2D09-E411-A993-0025903451A8.root",
        #"/store/mc/Spring14miniaod//GluGluToHToTauTau_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/A871CFAE-2D09-E411-B079-0025900EB232.root",
        #"/store/mc/Spring14miniaod//GluGluToHToTauTau_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/C214279B-2C09-E411-B39F-0025903451A8.root",
        #"/store/mc/Spring14miniaod//GluGluToHToTauTau_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/DA6898EF-2C09-E411-BE79-003048D410ED.root",
        #"/store/mc/Spring14miniaod//GluGluToHToTauTau_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/E42DB840-2E09-E411-9B3A-003048D410ED.root"
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
process.MessageLogger.cerr.FwkReport.reportEvery = 50
#process.MessageLogger.categories.append('onlyError')
#process.MessageLogger.cerr.onlyError=cms.untracked.PSet(threshold  = cms.untracked.string('ERROR'))
#process.MessageLogger.cerr.threshold='ERROR'
#process.MessageLogger = cms.Service("MessageLogger",
#	destinations = cms.untracked.vstring('log.txt')
#)
#process.MessageLogger.threshold = cms.untracked.string('ERROR')

