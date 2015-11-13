SVFITBYPASS = False
IsMC = XXX_ISMC_XXX
APPLYFSR = False

import os
import FWCore.ParameterSet.Config as cms
process = cms.Process("HTauTauNtuples")

PyFilePath = os.environ['CMSSW_BASE']+"/src/LLRHiggsTauTau/NtupleProducer/"
execfile(PyFilePath+"python/triggers.py") # contains the list of triggers and filters

### ----------------------------------------------------------------------
### Set the GT
### ----------------------------------------------------------------------
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
if IsMC:
    #process.GlobalTag.globaltag = 'PLS170_V6AN1::All'#'GR_70_V2_AN1::All'   #MC in 70X, cf https://twiki.cern.ch/twiki/bin/view/CMS/MiniAOD
    #process.GlobalTag.globaltag = 'PHYS14_25_V1::All' #MC in PHYS14
    #process.GlobalTag.globaltag = 'PHYS14_ST_V1::All'
    process.GlobalTag.globaltag = 'MCRUN2_74_V9::All' #MC in Spring15 miniAOD
else :
    process.GlobalTag.globaltag = 'GR_70_V2_AN1::All'   # data in 70X, cf https://twiki.cern.ch/twiki/bin/view/CMS/MiniAOD
print process.GlobalTag.globaltag

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

# use local file
#process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring(
#         '/store/user/lcadamur/HTauTauTrees/test_prod_mem_1apr/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/HHTauTauTrees_noNtupl_SVbypass/150402_090719/0000/Enriched_miniAOD_1.root',
#         )
#    )

#source external file with LFN of enriched miniAOD (VBF)
execfile(PyFilePath+"test/XXX_SAMPLEFILENAME_XXX")
process.source = cms.Source("PoolSource",
    fileNames = FILELIST
    )
process.load("FWCore.MessageService.MessageLogger_cfi")

#Limited nEv for testing purposes. -1 to run all events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(XXX_MAXEVENTS_XXX)
)
process.source.skipEvents = cms.untracked.uint32 (XXX_SKIPEVENTS_XXX)

# JSON mask for data --> to upload each time with the proper list using the macro in tools/
#from JSON file
# not sure it is needed here...
# CRAB has already filtered by the desired lumi
# if not IsMC:
#   process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange( *(
#       '251244:85-251244:86',
#       '251244:88-251244:93',
#       '251244:96-251244:121',
#       '251244:123-251244:156',
#       '251244:158-251244:428',
#       '251244:430-251244:442',
#       '251251:1-251251:31',
#       '251251:33-251251:97',
#       '251251:99-251251:167',
#       '251252:1-251252:283',
#       '251252:285-251252:505',
#       '251252:507-251252:554',
#       '251561:1-251561:94',
#       '251562:1-251562:439',
#       '251562:443-251562:691',
#       '251643:1-251643:216',
#       '251643:222-251643:606',
#       '251721:21-251721:36',
#       '251721:123-251721:244',
#       '251883:56-251883:56',
#       '251883:58-251883:60',
#       '251883:62-251883:144',
#       '251883:156-251883:437',
#   ))


#Global configuration
process.HTauTauTree = cms.EDAnalyzer("HTauTauNtuplizer",
                      fileName = cms.untracked.string ("CosaACaso"),
                      skipEmptyEvents = cms.bool(True),
                      applyFSR = cms.bool(APPLYFSR),
                      IsMC = cms.bool(IsMC),
                      triggerResultsLabel = cms.InputTag("TriggerResults", "", "HLT"),
                      triggerSet = cms.InputTag("selectedPatTrigger"),
                      triggerList = HLTLIST
                      )

if SVFITBYPASS:
    process.HTauTauTree.CandCollection = cms.untracked.string("SVbypass")
else:
    process.HTauTauTree.CandCollection = cms.untracked.string("SVllCand")

process.TFileService=cms.Service('TFileService',fileName=cms.string("XXX_OUTPUTFILENTUPLE_XXX"))
process.p = cms.Path(process.HTauTauTree)
