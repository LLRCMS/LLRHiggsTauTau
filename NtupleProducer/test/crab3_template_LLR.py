# TEMPLATE used for automatic script submission of multiple datasets

from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'DefaultReqName'
config.General.workArea = 'DefaultCrab3Area'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'analyzer_LLR.py' # to produce LLR ntuples or EnrichedMiniAOD according to the RunNtuplizer bool
config.JobType.sendExternalFolder = True #Needed until the PR including the Spring16 ele MVA ID is integrated in CMSSW/cms-data.
config.JobType.inputFiles=['JECUncertaintySources'] # FRA: adding to the sandobx the directory with JEC files (https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3FAQ#How_are_the_inputFiles_handled_i)
#config.JobType.maxMemoryMB=4000

config.section_("Data")
config.Data.inputDataset = '/my/precious/dataset'
config.Data.inputDBS = 'global'
config.Data.splitting = 'Automatic'
#config.Data.unitsPerJob = 28000 #number of events per jobs # 18K FOR SINGLE ELE, 10k for others, 28K for muons
config.Data.totalUnits = -1 #number of event
config.Data.publication = True
config.Data.outputDatasetTag = 'DefaultPublishName'

#private 2K cores
#config.section_("Debug")
#config.Debug.extraJDL = [ '+DESIRED_Sites="T3_IT_Opportunistic_hnsci"','+JOB_CMSSite="T3_IT_Opportunistic_hnsci"','+AccountingGroup="highprio.spiga"' ]
config.section_("Site")
# PARIGI
config.Site.storageSite = 'T2_FR_GRIF_LLR'
# MILANO
#config.Site.storageSite = 'T3_IT_MIB'
