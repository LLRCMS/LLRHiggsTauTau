from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'HTauTau_TTbar_NoSVFit_12Giu2015'
config.General.workArea = 'crab3'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'analyzer.py'

config.section_("Data")
config.Data.inputDataset = '/TTTo2L2Nu_13TeV-powheg/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 4000 #number of files per jobs
config.Data.totalUnits = -1 #number of event
config.Data.outLFNDirBase = '/store/user/lcadamur/HHTauTauFrameworkNtuples/HTauTau_TTbar_NoSVFit_12Giu2015/'
config.Data.publication = False
config.Data.publishDataName = 'HTauTau_TTbar_NoSVFit_12Giu2015'

config.section_("Site")
config.Site.storageSite = 'T2_FR_GRIF_LLR'
