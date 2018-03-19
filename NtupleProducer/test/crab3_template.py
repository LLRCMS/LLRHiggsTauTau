# TEMPLATE used for automatic script submission of multiple datasets

from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'DefaultReqName'
config.General.workArea = 'DefaultCrab3Area'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'analyzer.py' # to produce LLR ntuples or EnrichedMiniAOD according to the RunNtuplizer bool
config.JobType.sendExternalFolder = True #Needed until the PR including the Spring16 ele MVA ID is integrated in CMSSW/cms-data.

config.section_("Data")
config.Data.inputDataset = '/my/precious/dataset'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 20000 #number of events per jobs # 18K FOR SINGLE ELE, 10k for others
config.Data.totalUnits = -1 #number of event
config.Data.outLFNDirBase = '/store/user/lcadamur/HHNtuples/DefaultOutLFNDirBase'
config.Data.publication = True
config.Data.outputDatasetTag = 'DefaultPublishName'

config.section_("Site")
# PARIGI
#config.Site.storageSite = 'T2_FR_GRIF_LLR'
# MILANO
config.Site.storageSite = 'T3_IT_MIB'
