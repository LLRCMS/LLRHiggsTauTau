# TEMPLATE used for automatic script submission of multiple datasets

from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'DefaultReqName'
config.General.workArea = 'DefaultCrab3Area'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'analyzer.py'

config.section_("Data")
config.Data.inputDataset = '/my/precious/dataset'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 4000 #number of events per jobs
config.Data.totalUnits = -1 #number of event
config.Data.outLFNDirBase = '/store/user/lcadamur/HHNtuples/DefaultOutLFNDirBase'
config.Data.publication = True
config.Data.publishDataName = 'DefaultPublishName'

config.section_("Site")
config.Site.storageSite = 'T2_FR_GRIF_LLR'
