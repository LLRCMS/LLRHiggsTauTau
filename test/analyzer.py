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

##
## Standard sequence
##
execfile(PyFilePath+"python/HiggsTauTauProducer.py")

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
