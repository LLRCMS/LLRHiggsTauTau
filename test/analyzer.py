#-----------------------------------------
#
#Producer controller
#
#-----------------------------------------
import os
PyFilePath = os.environ['CMSSW_BASE']+"/src/LLRHiggsTauTau/"

#samples list (it could be moved to a cfg file for better reading
#samples = [
#]
APPLYMUCORR=False
##
## Standard sequence
##
execfile(PyFilePath+"python/HiggsTauTauProducer.py")

### ----------------------------------------------------------------------
### Source, better to use sample to run on batch
### ----------------------------------------------------------------------
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'/store/cmst3/user/cmgtools/CMG/GluGluToHToZZTo4L_M-130_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/V5/PAT_CMG_V5_2_0/patTuple_1.root'
        "/store/mc/Spring14miniaod//GluGluToHToTauTau_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/1062D5FF-2D09-E411-943C-0025900EB52A.root",
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
process.TFileService=cms.Service('TFileService',fileName=cms.string('HTauTauAnalysis.root'))

#Global configuration
TreeSetup = cms.EDAnalyzer("HTauTauNtupleMaker",
                           )

process.HTauTauTree = TreeSetup.clone()

process.p = cms.EndPath(process.HTauTauTree)


# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
