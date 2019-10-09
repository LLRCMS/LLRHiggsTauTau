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

# SingleElectron dataset :
'Oct19v1_Data_2017_Oct19_SingleEle_BlockB',
'Oct19v1_Data_2017_Oct19_SingleEle_BlockC',
'Oct19v1_Data_2017_Oct19_SingleEle_BlockD',
'Oct19v1_Data_2017_Oct19_SingleEle_BlockE',
'Oct19v1_Data_2017_Oct19_SingleEle_BlockF',
# SingleMuon dataset 
'Oct19v1_Data_2017_Oct19_SingleMu_BlockB',
'Oct19v1_Data_2017_Oct19_SingleMu_BlockC',
'Oct19v1_Data_2017_Oct19_SingleMu_BlockD',
'Oct19v1_Data_2017_Oct19_SingleMu_BlockE',
'Oct19v1_Data_2017_Oct19_SingleMu_BlockF',
# DoubleEG dataset
'Oct19v1_Data_2017_Oct19_DoubleEG_BlockB',
'Oct19v1_Data_2017_Oct19_DoubleEG_BlockC',
'Oct19v1_Data_2017_Oct19_DoubleEG_BlockD',
'Oct19v1_Data_2017_Oct19_DoubleEG_BlockE',
'Oct19v1_Data_2017_Oct19_DoubleEG_BlockF',
# DoubleMuon
'Oct19v1_Data_2017_Oct19_DoubleMu_BlockB',
'Oct19v1_Data_2017_Oct19_DoubleMu_BlockC',
'Oct19v1_Data_2017_Oct19_DoubleMu_BlockD',
'Oct19v1_Data_2017_Oct19_DoubleMu_BlockE',
'Oct19v1_Data_2017_Oct19_DoubleMu_BlockF',
# MuonEG
'Oct19v1_Data_2017_Oct19_MuonEG_BlockB',
'Oct19v1_Data_2017_Oct19_MuonEG_BlockC',
'Oct19v1_Data_2017_Oct19_MuonEG_BlockD',
'Oct19v1_Data_2017_Oct19_MuonEG_BlockE',
'Oct19v1_Data_2017_Oct19_MuonEG_BlockF',
# Tau
'Oct19v1_Data_2017_Oct19_Tau_BlockB',
'Oct19v1_Data_2017_Oct19_Tau_BlockC',
'Oct19v1_Data_2017_Oct19_Tau_BlockD',
'Oct19v1_Data_2017_Oct19_Tau_BlockE',
'Oct19v1_Data_2017_Oct19_Tau_BlockF'

   ]
 
datasetinputs = [

# SingleElectron dataset : AT LEAST 1 high-energy electron in the event.
'/SingleElectron/Run2017B-31Mar2018-v1/MINIAOD',
'/SingleElectron/Run2017C-31Mar2018-v1/MINIAOD',
'/SingleElectron/Run2017D-31Mar2018-v1/MINIAOD',
'/SingleElectron/Run2017E-31Mar2018-v1/MINIAOD',
'/SingleElectron/Run2017F-31Mar2018-v1/MINIAOD',
# SingleMuon dataset : AT LEAST 1 high-energy muon in the event.
'/SingleMuon/Run2017B-31Mar2018-v1/MINIAOD',
'/SingleMuon/Run2017C-31Mar2018-v1/MINIAOD',
'/SingleMuon/Run2017D-31Mar2018-v1/MINIAOD',
'/SingleMuon/Run2017E-31Mar2018-v1/MINIAOD',
'/SingleMuon/Run2017F-31Mar2018-v1/MINIAOD',
# DoubleEG dataset : AT LEAST 2 high-energy electron in the event.
'/DoubleEG/Run2017B-31Mar2018-v1/MINIAOD',
'/DoubleEG/Run2017C-31Mar2018-v1/MINIAOD',
'/DoubleEG/Run2017D-31Mar2018-v1/MINIAOD',
'/DoubleEG/Run2017E-31Mar2018-v1/MINIAOD',
'/DoubleEG/Run2017F-31Mar2018-v1/MINIAOD',
# DoubleMuon dataset : AT LEAST 2 high-energy muon in the event.
'/DoubleMuon/Run2017B-31Mar2018-v1/MINIAOD',
'/DoubleMuon/Run2017C-31Mar2018-v1/MINIAOD',
'/DoubleMuon/Run2017D-31Mar2018-v1/MINIAOD',
'/DoubleMuon/Run2017E-31Mar2018-v1/MINIAOD',
'/DoubleMuon/Run2017F-31Mar2018-v1/MINIAOD',
# MuonEG dataset : AT LEAST 1 high-energy electron and 1 high-energy muon in the event.
'/MuonEG/Run2017B-31Mar2018-v1/MINIAOD',
'/MuonEG/Run2017C-31Mar2018-v1/MINIAOD',
'/MuonEG/Run2017D-31Mar2018-v1/MINIAOD',
'/MuonEG/Run2017E-31Mar2018-v1/MINIAOD',
'/MuonEG/Run2017F-31Mar2018-v1/MINIAOD',
# Tau dataset : AT LEAST 1 high-energy tau
'/Tau/Run2017B-31Mar2018-v1/MINIAOD',
'/Tau/Run2017C-31Mar2018-v1/MINIAOD',
'/Tau/Run2017D-31Mar2018-v1/MINIAOD',
'/Tau/Run2017E-31Mar2018-v1/MINIAOD',
'/Tau/Run2017F-31Mar2018-v1/MINIAOD'
                
    ]

for d in range(0,len(datasetnames)):

    print 'multicrab.py: Running datasetname: ', datasetnames[d]

    config.section_('General')
    config.General.requestName  = datasetnames[d]
    config.General.workArea     = 'crab3'
    config.General.transferLogs = True

    config.section_('JobType')
    config.JobType.pluginName  = 'Analysis'
    config.JobType.psetName    = 'analyzer_data2017.py'
    config.JobType.sendExternalFolder = True
    config.JobType.maxMemoryMB = 2000 # Default == 2Gb : maximum guaranteed to run on all sites
    config.JobType.inputFiles  = ['JECUncertaintySources']

    config.section_('Data')
    config.Data.allowNonValidInputDataset = True
    config.Data.inputDataset   = datasetinputs[d]
    config.Data.inputDBS       = 'global'
    config.Data.splitting      = 'LumiBased'
    config.Data.unitsPerJob    = 30    
    config.Data.lumiMask       = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt'
    config.Data.outLFNDirBase  = '/store/user/cmartinp/ttH_Legacy/Data_2016_Oct19/'
    config.Data.publication    = True
    config.Data.outputDatasetTag = datasetnames[d] 

    print 'multicrab.py: outLFNDirBase = /store/user/cmartinp/ttH_Legacy/Data_2016_Oct19/'

    config.section_('Site')
    config.Site.storageSite    = 'T2_FR_GRIF_IRFU' #'T2_FR_GRIF_LLR'
    print 'multicrab.py: Submitting Jobs'
    submit(config)
