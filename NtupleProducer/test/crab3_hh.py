from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'HHTauTau_TTJets'
config.General.workArea = 'crab3'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'analyzer.py'

config.section_("Data")
config.Data.inputDataset = '/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
#config.Data.inputDataset = '/TTbarH_HToTauTau_M-125_13TeV_amcatnlo-pythia8-tauola/Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased' #eventBased for events
config.Data.unitsPerJob = 1 #number of files per jobs
config.Data.totalUnits = -1 #number of event
config.Data.outLFN = '/store/user/gortona/HTauTauTrees/prod_2015_02_03/'
config.Data.publication = False
config.Data.publishDataName = 'HHTauTauTrees'

config.section_("Site")
config.Site.storageSite = 'T2_FR_GRIF_LLR'
