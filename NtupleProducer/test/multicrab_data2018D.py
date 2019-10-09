if __name__ == '__main__':
 #####
 ##   Multicrab configuration
 #####
 import sys
 from CRABClient.UserUtilities import config, getUsernameFromSiteDB
 config = config()
 from CRABAPI.RawCommand import crabCommand
 from CRABClient.ClientExceptions import ClientException
 from httplib import HTTPException
 config.General.workArea = 'Crab_projects'

 def submit(config):
  try:
   crabCommand('submit', config = config)
  except HTTPException as hte:
   print "Failed submitting task: %s" % (hte.headers)
  except ClientException as cle:
   print "Failed submitting task: %s" % (cle)
 
 #####
 ##   Crab configuration
 #####
 
 datasetnames  = [

 # EGamma dataset (SingleElectron + SinglePhoton + DoubleEG)
 'Oct19v1_Data_2018_Oct19_EGamma_BlockD',
 # SingleMuon dataset 
 'Oct19v1_Data_2018_Oct19_SingleMu_BlockD',
 # DoubleMuon
 'Oct19v1_Data_2018_Oct19_DoubleMu_BlockD',
 # MuonEG
 'Oct19v1_Data_2018_Oct19_MuonEG_BlockD',
 # Tau
 'Oct19v1_Data_2018_Oct19_Tau_BlockD',

   ]
 
datasetinputs = [

 # EGamma dataset (SingleElectron + SinglePhoton + DoubleEG)
 '/EGamma/Run2018D-22Jan2019-v2/MINIAOD',
 # SingleMuon dataset : AT LEAST 1 high-energy muon in the event.
 '/SingleMuon/Run2018D-22Jan2019-v2/MINIAOD',
 # DoubleMuon dataset : AT LEAST 2 high-energy muon in the event.
 '/DoubleMuon/Run2018D-PromptReco-v2/MINIAOD',
 # MuonEG dataset : AT LEAST 1 high-energy electron and 1 high-energy muon in the event.
 '/MuonEG/Run2018D-PromptReco-v2/MINIAOD',
 # Tau dataset : AT LEAST 1 high-energy tau
 '/Tau/Run2018D-PromptReco-v2/MINIAOD'
                
    ]

for d in range(0,len(datasetnames)):

    print 'multicrab.py: Running datasetname: ', datasetnames[d]

    config.section_('General')
    config.General.requestName  = datasetnames[d]
    config.General.workArea     = 'crab3'
    config.General.transferLogs = True

    config.section_('JobType')
    config.JobType.pluginName  = 'Analysis'
    config.JobType.psetName    = 'analyzer_data2018D.py'
    config.JobType.sendExternalFolder = True
    config.JobType.maxMemoryMB = 2000 # Default == 2Gb : maximum guaranteed to run on all sites
    config.JobType.inputFiles  = ['JECUncertaintySources']

    config.section_('Data')
    config.Data.allowNonValidInputDataset = True
    config.Data.inputDataset   = datasetinputs[d]
    config.Data.inputDBS       = 'global'
    config.Data.splitting      = 'LumiBased'
    config.Data.unitsPerJob    = 10
    config.Data.lumiMask       = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'
    config.Data.outLFNDirBase  = '/store/user/cmartinp/ttH_Legacy/Data_2018_Oct19/'
    config.Data.publication    = True
    config.Data.outputDatasetTag = datasetnames[d] 

    print 'multicrab.py: outLFNDirBase = /store/user/cmartinp/ttH_Legacy/Data_2018_Oct19/'

    config.section_('Site')
    config.Site.storageSite    = 'T2_FR_GRIF_IRFU' #'T2_FR_GRIF_LLR'
    print 'multicrab.py: Submitting Jobs'
    submit(config)
