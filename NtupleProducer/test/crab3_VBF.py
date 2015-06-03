from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'HTauTau_VBF'
config.General.workArea = 'crab3'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'analyzer.py'

config.section_("Data")
config.Data.inputDataset = '/VBF_HToTauTau_M-125_13TeV-powheg-pythia6/Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v2/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 4000 #number of files per jobs
config.Data.totalUnits = -1 #number of event
config.Data.outLFNDirBase = '/store/user/lcadamur/EnrichedMiniAOD/VBF_HToTauTau_pfMET_prod_03_06_2015/'
config.Data.publication = True
config.Data.publishDataName = 'VBF_HTauTau'

config.section_("Site")
config.Site.storageSite = 'T2_FR_GRIF_LLR'
