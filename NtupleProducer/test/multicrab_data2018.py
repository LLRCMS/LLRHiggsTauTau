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
'Data_2018_v1_EGamma_BlockA',
'Data_2018_v1_EGamma_BlockB',
'Data_2018_v1_EGamma_BlockC',
'Data_2018_v1_EGamma_BlockD',
# SingleMuon dataset 
'Data_2018_v1_SingleMu_BlockA',
'Data_2018_v1_SingleMu_BlockB',
'Data_2018_v1_SingleMu_BlockC',
'Data_2018_v1_SingleMu_BlockD',
# DoubleMuon
'Data_2018_v1_DoubleMu_BlockA',
'Data_2018_v1_DoubleMu_BlockB',
'Data_2018_v1_DoubleMu_BlockC',
'Data_2018_v1_DoubleMu_BlockD',
# MuonEG
'Data_2018_v1_MuonEG_BlockA',
'Data_2018_v1_MuonEG_BlockB',
'Data_2018_v1_MuonEG_BlockC',
'Data_2018_v1_MuonEG_BlockD',
# Tau
'Data_2018_v1_Tau_BlockA',
'Data_2018_v1_Tau_BlockB',
'Data_2018_v1_Tau_BlockC',
'Data_2018_v1_Tau_BlockD',

   ]
 
datasetinputs = [

# EGamma dataset (SingleElectron + SinglePhoton + DoubleEG)
'/EGamma/Run2018A-17Sep2018-v2/MINIAOD',
'/EGamma/Run2018B-17Sep2018-v1/MINIAOD',
'/EGamma/Run2018C-17Sep2018-v1/MINIAOD',
'/EGamma/Run2018D-PromptReco-v2/MINIAOD',
# SingleMuon dataset : AT LEAST 1 high-energy muon in the event.
'/SingleMuon/Run2018A-17Sep2018-v2/MINIAOD',
'/SingleMuon/Run2018B-17Sep2018-v1/MINIAOD',
'/SingleMuon/Run2018C-17Sep2018-v1/MINIAOD',
'/SingleMuon/Run2018D-PromptReco-v2/MINIAOD',
# DoubleMuon dataset : AT LEAST 2 high-energy muon in the event.
'/DoubleMuon/Run2018A-17Sep2018-v2/MINIAOD',
'/DoubleMuon/Run2018B-17Sep2018-v1/MINIAOD',
'/DoubleMuon/Run2018C-17Sep2018-v1/MINIAOD',
'/DoubleMuon/Run2018D-PromptReco-v2/MINIAOD',
# MuonEG dataset : AT LEAST 1 high-energy electron and 1 high-energy muon in the event.
'/MuonEG/Run2018A-17Sep2018-v1/MINIAOD',
'/MuonEG/Run2018B-17Sep2018-v1/MINIAOD',
'/MuonEG/Run2018C-17Sep2018-v1/MINIAOD',
'/MuonEG/Run2018D-PromptReco-v2/MINIAOD',
# Tau dataset : AT LEAST 1 high-energy tau
'/Tau/Run2018A-17Sep2018-v1/MINIAOD',
'/Tau/Run2018B-17Sep2018-v1/MINIAOD',
'/Tau/Run2018C-17Sep2018-v1/MINIAOD',
'/Tau/Run2018D-PromptReco-v2/MINIAOD'
                
    ]

for d in range(0,len(datasetnames)):

    print 'multicrab.py: Running datasetname: ', datasetnames[d]

    config.section_('General')
    config.General.requestName = datasetnames[d]
    config.General.workArea    = 'crab3'
    config.General.transferLogs = True

    config.section_('JobType')
    config.JobType.pluginName  = 'Analysis'
    config.JobType.psetName    = 'analyzer_2018.py'
    config.JobType.sendExternalFolder = True
    config.JobType.inputFiles  = ['JECUncertaintySources']

    config.section_('Data')
    config.Data.inputDataset   = datasetinputs[d]
    config.Data.inputDBS       = 'global'
    config.Data.splitting      = 'LumiBased'
    config.Data.unitsPerJob    = 30
    config.Data.lumiMask       = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'
    config.Data.outLFNDirBase  = '/store/user/cmartinp/ttH_Legacy/Data_2018_v1/'
    config.Data.publication    = True
    config.Data.outputDatasetTag = datasetnames[d] 

    print 'multicrab.py: outLFNDirBase = /store/user/cmartinp/'

    config.section_('Site')
    config.Site.storageSite    = 'T2_FR_GRIF_LLR'
    print 'multicrab.py: Submitting Jobs'
    submit(config)
