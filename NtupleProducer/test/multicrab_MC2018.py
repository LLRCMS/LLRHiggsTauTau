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

# signal
'MC_2018_Oct19_ttHToNonbb', #0
'MC_2018_Oct19_ttHJetToNonbb', #1
'MC_2018_Oct19_TTH_ctcvcp', #2
# TH
'MC_2018_Oct19_THQ_ctcvcp', #3
'MC_2018_Oct19_THW_ctcvcp', #4
# TTZ
'MC_2018_Oct19_TTZ_M10', #5
'MC_2018_Oct19_TTZ_M1to10', #6
# TTW
'MC_2018_Oct19_TTWJets', #7
# TTWW
'MC_2018_Oct19_TTWW', #8
# ST
'MC_2018_Oct19_ST_sCh_lepDecay', #9
'MC_2018_Oct19_ST_tCh_top', #10
'MC_2018_Oct19_ST_tCh_antitop', #11
'MC_2018_Oct19_ST_tW_top', #12
'MC_2018_Oct19_ST_tW_antitop', #13
'MC_2018_Oct19_tWll', #14
# TT
'MC_2018_Oct19_TTTo2L', #15
'MC_2018_Oct19_TTToSemiLep', #16
'MC_2018_Oct19_TTToHad', #17
'MC_2018_Oct19_TTJets_DiLep', #18
'MC_2018_Oct19_TTJets_TToSingleLep', #19
'MC_2018_Oct19_TTJets_TbarToSingleLep', #20
# Triboson
'MC_2018_Oct19_WWW', #21
'MC_2018_Oct19_WWZ', #22
'MC_2018_Oct19_WZZ', #23
'MC_2018_Oct19_ZZZ', #24
# Conversions
'MC_2018_Oct19_WZG', #25
'MC_2018_Oct19_WGToLNuG', #26
'MC_2018_Oct19_ZGToLLG',  #27
'MC_2018_Oct19_TGJetsLep', #28
'MC_2018_Oct19_TTGJets', #29
# ggF
'MC_2018_Oct19_GGHToTauTau', #30
# VBF
'MC_2018_Oct19_VBFHToTauTau', #31
# tZq
'MC_2018_Oct19_tZq', #32
# WpWp
'MC_2018_Oct19_WpWpJJ', #33
# WW DS
'MC_2018_Oct19_WW_DS', #34
# TTTT
'MC_2018_Oct19_TTTT', #35
# VH
'MC_2018_Oct19_VHToNonbb', #36
'MC_2018_Oct19_ZHTobb', #37
'MC_2018_Oct19_ZHTobb_ext', #38
'MC_2018_Oct19_ZHToTauTau', #39
# DY
'MC_2018_Oct19_DYJets_M10to50', #40
'MC_2018_Oct19_DYJets_M50', #41
'MC_2018_Oct19_DYJets_M50_ext', #42
#EWK
'MC_2018_Oct19_WJets', #43
'MC_2018_Oct19_WWTo2LNu', #44
'MC_2018_Oct19_WZTo3LNu_ext', #45
'MC_2018_Oct19_ZZTo4L_ext1', #46
'MC_2018_Oct19_ZZTo4L_ext2', #47
# Training
'MC_2018_Oct19_ttZ', #48
'MC_2018_Oct19_ttW', #49

    ]

 datasetinputs = [
 
# signal
'/ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM', #0
'/ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #1
'/TTH_4f_ctcvcp_TuneCP5_13TeV_madgraph_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #2
# TH
'/THQ_ctcvcp_4f_Hincl_13TeV_madgraph_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #3
'/THW_ctcvcp_5f_Hincl_13TeV_madgraph_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #4
# TTZ
'/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM', #5
'/TTZToLL_M-1to10_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #6
# TTW
'/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM', #7
# TTWW
'/TTWW_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v1/MINIAODSIM', #8
# ST
'/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v4/MINIAODSIM', #9
'/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #10
'/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #11
'/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM', #12
'/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM', #13
'/ST_tWll_5f_LO_TuneCP5_PSweights_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM', #14
# TT
'/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #15
'/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #16
'/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #17
'/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #18
'/TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #19
'/TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #20
# Triboson
'/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM', #21
'/WWZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM', #22
'/WZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM', #23
'/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM', #24
# Conversions
'/WZG_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #25
'/WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #26
'/ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM', #27
'/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM', #28
'/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #29
# ggF
'/GluGluHToTauTau_M125_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM', #30
# VBF
'/VBFHToTauTau_M125_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM', #31
# tZq
'/tZq_ll_4f_ckm_NLO_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM', #32
# WpWp
'/WpWpJJ_EWK-QCD_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #33
# WW DS
'/WWTo2L2Nu_DoubleScattering_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #34
# TTTT
'/TTTT_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM', #35
# VH
'/VHToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM', #36
'/ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM', #37
'/ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM', #38
'/ZHToTauTau_M125_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM', #39
#DY
'/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM', #40
'/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #41
'/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v1/MINIAODSIM', #42
#EWK
'/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM', #43
'/WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #44
'/WZTo3LNu_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM', #45
'/ZZTo4L_TuneCP5_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM', #46
'/ZZTo4L_TuneCP5_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v2/MINIAODSIM', #47
# Training
'/ttZJets_TuneCP5_13TeV_madgraphMLM_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM', #48
'/ttWJets_TuneCP5_13TeV_madgraphMLM_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM', #49

    ]


for d in range(0,len(datasetnames)):
    print 'multicrab.py: Running datasetname: ', datasetnames[d]

    config.section_('General')
    config.General.requestName = datasetnames[d]
    config.General.workArea    = 'crab3_Oct19'
    config.General.transferLogs = True

    config.section_('JobType')
    config.JobType.pluginName  = 'Analysis'
    config.JobType.psetName    = 'analyzer_MC2018.py'
    config.JobType.sendExternalFolder = True
    config.JobType.inputFiles=['JECUncertaintySources']

    config.section_("Data")
    config.Data.inputDataset   = datasetinputs[d]
    config.Data.inputDBS = 'global'
    config.JobType.maxMemoryMB = 2000 # Default == 2Gb : maximum guaranteed to run on all sites
    config.Data.inputDBS       = 'global'
    #config.Data.splitting      = 'FileBased'
    #config.Data.totalUnits     = 40000 #With 'FileBased' splitting tells how many files to analyse
    config.Data.splitting      = 'Automatic'
    config.Data.unitsPerJob    = 180
    config.Data.unitsPerJob    = 1
    config.Data.outLFNDirBase  = '/store/user/cmartinp/ttH_Legacy/MC_2018_Oct19/'
    config.Data.publication = True
    config.Data.outputDatasetTag = datasetnames[d]    

    print 'multicrab.py: outLFNDirBase = /store/user/cmartinp/ttH_Legacy/MC_2018_Oct19/'
    #config.Data.publication = True

    config.section_('Site')
    config.Site.storageSite    = 'T2_FR_GRIF_LLR'
    print 'multicrab.py: Submitting Jobs'
    submit(config)
