from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()
config.section_('General')
config.General.requestName = 'SUSYggH160_v2'
config.General.workArea = 'crab_projects'
config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'analyzer.py'
config.section_("Data")
config.Data.inputDataset = '/SUSYGluGluToHToTauTau_M-160_TuneCUETP8M1_13TeV-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 20000 #number of files per jobs
config.Data.totalUnits = -1 #number of event
config.Data.outLFNDirBase = '/store/user/ccaillol/Ntuples_76X_v1/SUSYggH160_27_01_16_v2/'
config.Data.publication = False
config.Data.outputDatasetTag = 'SUSYggH160_27_01_16_v2'
config.section_("Site")
config.Site.storageSite = 'T2_BE_IIHE'
config.Site.blacklist = ['T1_US_FNAL']
