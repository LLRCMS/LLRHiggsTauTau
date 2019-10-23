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
config.JobType.inputFiles=['JECUncertaintySources'] # FRA: adding to the sandobx the directory with JEC files (https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3FAQ#How_are_the_inputFiles_handled_i)
config.JobType.allowUndistributedCMSSW = True #Davide June 2019: This line is necessary to run with CMSSW_10_2_14 on slc7_amd64_gcc700
#config.JobType.maxJobRuntimeMin = 2200 #Davide June 2019: This is necesessary for the analysis of FileBased splitted dataset (2018 DY)
#config.JobType.maxMemoryMB = 2500 
#config.JobType.maxMemoryMB = 8000 #Davide June 2019: This is necesessary for the analysis of FileBased splitted dataset (2018 DY)

config.section_("Data")
config.Data.inputDataset = '/my/precious/dataset'
config.Data.inputDBS = 'global'
#config.Data.splitting = 'FileBased' # Use this split algorithm for huge datasets
#config.Data.splitting = 'EventAwareLumiBased'
config.Data.splitting = 'Automatic'
#config.Data.unitsPerJob = 18000 #number of events per jobs when splitting mode is EventAwareLumiBased # Number of files per job when splitting mode is FileBased  # 18K FOR SOME TT BKG, 10k for others
config.Data.totalUnits = -1 #number of event
config.Data.outLFNDirBase = '/store/user/lcadamur/HHNtuples/DefaultOutLFNDirBase'
config.Data.publication = True
config.Data.outputDatasetTag = 'DefaultPublishName'
config.Data.allowNonValidInputDataset = True

# to run on dedicated 2k cores
#config.section_("Debug")
#config.Debug.extraJDL = [ '+DESIRED_Sites="T3_IT_Opportunistic_hnsci"','+JOB_CMSSite="T3_IT_Opportunistic_hnsci"','+AccountingGroup="highprio.spiga"' ]


config.section_("Site")
# PARIGI
#config.Site.storageSite = 'T2_FR_GRIF_LLR'
# MILANO
config.Site.storageSite = 'T3_IT_MIB'
