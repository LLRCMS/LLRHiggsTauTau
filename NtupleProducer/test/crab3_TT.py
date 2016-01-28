from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()
config.section_('General')
config.General.requestName = 'TT_v1'
config.General.workArea = 'crab_projects'
config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'analyzer.py'
config.section_("Data")
config.Data.inputDataset = '/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 20000 #number of files per jobs
config.Data.totalUnits = -1 #number of event
config.Data.outLFNDirBase = '/store/user/ccaillol/Ntuples_76X_v1/TT_28_01_16/'
config.Data.publication = False
config.Data.outputDatasetTag = 'TT_28_01_16'
config.section_("Site")
config.Site.storageSite = 'T2_BE_IIHE'
config.Site.blacklist = ['T1_US_FNAL']
