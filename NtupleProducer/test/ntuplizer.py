SVFITBYPASS = False
IsMC = True
APPLYFSR = False

import FWCore.ParameterSet.Config as cms
process = cms.Process("HTauTauNtuples")

### ----------------------------------------------------------------------
### Set the GT
### ----------------------------------------------------------------------
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
if IsMC:
    #process.GlobalTag.globaltag = 'PLS170_V6AN1::All'#'GR_70_V2_AN1::All'   #MC in 70X, cf https://twiki.cern.ch/twiki/bin/view/CMS/MiniAOD
    process.GlobalTag.globaltag = 'PHYS14_25_V1::All' #MC in PHYS14
    #process.GlobalTag.globaltag = 'PHYS14_ST_V1::All'
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
import os
PyFilePath = os.environ['CMSSW_BASE']+"/src/LLRHiggsTauTau/NtupleProducer/"
execfile(PyFilePath+"test/VBF_EnrichMiniAOD_11Mag2015.py")
process.source = cms.Source("PoolSource",
    fileNames = FILELIST
    )
#process.load("FWCore.MessageService.MessageLogger_cfi")

#Global configuration
process.Ntuplizer = cms.EDAnalyzer("HTauTauNtuplizer",
                      fileName = cms.untracked.string ("CosaACaso"),
                      skipEmptyEvents = cms.bool(True),
                      applyFSR = cms.bool(APPLYFSR),
                      IsMC = cms.bool(IsMC),
                      triggerResultsLabel = cms.InputTag("TriggerResults", "", "HLT"),
                      triggerSet = cms.InputTag("selectedPatTrigger")
                      )

if SVFITBYPASS:
    process.Ntuplizer.CandCollection = cms.untracked.string("SVbypass")
else:
    process.Ntuplizer.CandCollection = cms.untracked.string("SVllCand")

process.TFileService=cms.Service('TFileService',fileName=cms.string('HTauTauAnalysis_ntuple.root'))
process.p = cms.Path(process.Ntuplizer)
