# coding: utf-8

from WMCore.Configuration import Configuration


config = Configuration()

config.section_("General")
config.General.requestName = "OVERWRITTEN"
config.General.workArea = "OVERWRITTEN"

config.section_("JobType")
config.JobType.pluginName = "Analysis"
config.JobType.psetName = "analyzer.py"
config.JobType.pyCfgParams = ["year=2017", "mc=True", "maxEvents=-1"]
config.JobType.sendExternalFolder = True
config.JobType.inputFiles = ["JECUncertaintySources"]
config.JobType.allowUndistributedCMSSW = True
# config.JobType.maxJobRuntimeMin = 2200 #Davide June 2019: This is necesessary for the analysis of FileBased splitted dataset (2018 DY)
# config.JobType.maxMemoryMB = 2500
# config.JobType.maxMemoryMB = 8000 #Davide June 2019: This is necesessary for the analysis of FileBased splitted dataset (2018 DY)

config.section_("Data")
config.Data.inputDataset = "OVERWRITTEN"
config.Data.inputDBS = "global"
config.Data.ignoreLocality = False
config.Data.splitting = "FileBased"
# config.Data.splitting = "EventAwareLumiBased"
# config.Data.splitting = "Automatic"
config.Data.unitsPerJob = 1
config.Data.totalUnits = -1
config.Data.outLFNDirBase = "OVERWRITTEN"
config.Data.publication = False
config.Data.outputDatasetTag = "DefaultPublishName"
config.Data.allowNonValidInputDataset = True

config.section_("Site")
config.Site.storageSite = "T2_DE_DESY"

config.section_("Debug")
config.Debug.extraJDL = ["+CMS_ALLOW_OVERFLOW=False"]
