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
'MC_2016_Oct19_TTHnobb', #0
'MC_2016_Oct19_ttHnobb', #1
# TH
'MC_2016_Oct19_THQ_TuneCP5_ctcvcp', #2
'MC_2016_Oct19_THQ_Tune8M1_ctcvcp', #3
'MC_2016_Oct19_THQ', #4
'MC_2016_Oct19_THW_TuneCP5_ctcvcp', #5
'MC_2016_Oct19_THW_Tune8M1_ctcvcp', #6
'MC_2016_Oct19_THW', #7
# TTZ
'MC_2016_Oct19_TTZ_M1to10', #8
'MC_2016_Oct19_TTZ_M10_ext1', #9
'MC_2016_Oct19_TTZ_M10_ext2', #10
# TTW
'MC_2016_Oct19_TTWJets', #11
# TTWW
'MC_2016_Oct19_TTWW', #12
# ST
'MC_2016_Oct19_ST_sCh_lepDecay', #13
'MC_2016_Oct19_ST_sCh_lepDecay_PS', #14
'MC_2016_Oct19_ST_tCh_top', #15
'MC_2016_Oct19_ST_tCh_antitop', #16
'MC_2016_Oct19_ST_tCh_antitop_PS', #17
'MC_2016_Oct19_ST_tW_top', #18
'MC_2016_Oct19_ST_tW_antitop', #19
'MC_2016_Oct19_tWll', #20 
# TT
'MC_2016_Oct19_TTTo2L_PS', #21
'MC_2016_Oct19_TTToSemiLep_PS', #22
'MC_2016_Oct19_TTToHad_PS', #23
'MC_2016_Oct19_TTJets_DiLep_v1', #24
'MC_2016_Oct19_TTJets_DiLep_ext', #25
'MC_2016_Oct19_TTJets_TToSingleLep_v1', #26
'MC_2016_Oct19_TTJets_TToSingleLep_ext', #27
'MC_2016_Oct19_TTJets_TbarToSingleLep_v1', #28
'MC_2016_Oct19_TTJets_TbarToSingleLep_ext', #29
# Triboson
'MC_2016_Oct19_WWW', #30 
'MC_2016_Oct19_WWZ', #31
'MC_2016_Oct19_WZZ', #32
'MC_2016_Oct19_ZZZ', #33
# Conversions
'MC_2016_Oct19_WZG', #34
'MC_2016_Oct19_WGToLNuG_ext1', #35
'MC_2016_Oct19_WGToLNuG_ext2', #36
'MC_2016_Oct19_WGToLNuG_ext3', #37
'MC_2016_Oct19_ZGToLLG',  #38
'MC_2016_Oct19_TGJetsLep', #39
'MC_2016_Oct19_TTGJets_v1', #40
'MC_2016_Oct19_TTGJets_ext', #41
# ggF
'MC_2016_Oct19_GGHToTauTau', #42
# VBF
'MC_2016_Oct19_VBFHToTauTau', #43
# tZq
'MC_2016_Oct19_tZq_ext', #44
'MC_2016_Oct19_tZq_PS', #45
# WpWp
'MC_2016_Oct19_WpWpJJ', #46
# WW DS
'MC_2016_Oct19_WW_DS', #47
# TTTT
'MC_2016_Oct19_TTTT', #48
# VH
'MC_2016_Oct19_VHToNonbb', #49
'MC_2016_Oct19_ZHTobb', #50
'MC_2016_Oct19_ZHToTauTau', #51
# DY
'MC_2016_Oct19_DYJets_M10to50', #52
'MC_2016_Oct19_DYJets_M50', #53
#EWK
'MC_2016_Oct19_WJets_v1', #54
'MC_2016_Oct19_WJets_ext', #55
'MC_2016_Oct19_WWTo2LNu', #56
'MC_2016_Oct19_WZTo3LNu', #57
'MC_2016_Oct19_ZZTo4L', #58
# Training
'MC_2016_Oct19_ttZ', #59
'MC_2016_Oct19_ttW', #60

    ]

 datasetinputs = [
 
# signal
'/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #0
'/ttHJetToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8_mWCutfix/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM', #1
# TH
'/THQ_ctcvcp_HIncl_M125_TuneCP5_13TeV-madgraph-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', #2
'/THQ_ctcvcp_Hincl_13TeV-madgraph-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #3
'/THQ_Hincl_13TeV-madgraph-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #4
'/THW_ctcvcp_HIncl_M125_TuneCP5_13TeV-madgraph-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', #5
'/THW_ctcvcp_Hincl_13TeV-madgraph-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #6
'/THW_Hincl_13TeV-madgraph-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #7
# TTZ
'/TTZToLL_M-1to10_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', #8
'/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/MINIAODSIM', #9
'/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext3-v1/MINIAODSIM', #10
# TTW
'/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/MINIAODSIM', #11
# TTWW
'/TTWW_TuneCUETP8M2T4_13TeV-madgraph-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/MINIAODSIM', #12
# ST
'/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', #13
'/ST_s-channel_4f_leptonDecays_13TeV_PSweights-amcatnlo-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', #14
'/ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', #15
'/ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', #16
'/ST_t-channel_antitop_4f_inclusiveDecays_13TeV_PSweights-powhegV2-madspin/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', #17
'/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/MINIAODSIM', #18
'/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/MINIAODSIM', #19
'/ST_tWll_5f_LO_13TeV-MadGraph-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', #20
# TT
'/TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', #21
'/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', #22
'/TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', #23
'/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #24
'/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM', #25
'/TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #26
'/TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM', #27
'/TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #28
'/TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM', #29
# Triboson
'/WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #30
'/WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #31
'/WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #32
'/ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #33
# Conversions
'/WZG_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', #34
'/WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/MINIAODSIM', #35
'/WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/MINIAODSIM', #36
'/WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext3-v1/MINIAODSIM', #37
'/ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/MINIAODSIM', #38
'/TGJets_leptonDecays_13TeV_amcatnlo_madspin_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #39
'/TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #40
'/TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #41
# ggF
'/GluGluHToTauTau_M125_13TeV_powheg_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #42
# VBF
'/VBFHToTauTau_M125_13TeV_powheg_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #43
# tZq
'/tZq_ll_4f_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/MINIAODSIM', #44
'/tZq_ll_4f_PSweights_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', #45
# WpWp
'/WpWpJJ_EWK-QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', #46
# WW DS
'/WWTo2L2Nu_DoubleScattering_13TeV-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #47
# TTTT
'/TTTT_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', # 48
# VH
'/VHToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #49
'/ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/MINIAODSIM', #50
'/ZHToTauTau_M125_13TeV_powheg_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #51
#DY
'/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #52
'/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/MINIAODSIM', #53
#EWK
'/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #54
'/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2/MINIAODSIM', #55
'/WWTo2L2Nu_13TeV-powheg/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #56
'/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/MINIAODSIM', #57
'/ZZTo4L_13TeV_powheg_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', #58
# Training
'/ttZJets_13TeV_madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', #59
'/ttWJets_13TeV_madgraphMLM/RunIISummer16MiniAODv3-94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', #60
    ]


for d in range(0,len(datasetnames)):
    print 'multicrab.py: Running datasetname: ', datasetnames[d]

    config.section_('General')
    config.General.requestName = datasetnames[d]
    config.General.workArea    = 'crab3_Oct19'
    config.General.transferLogs = True

    config.section_('JobType')
    config.JobType.pluginName  = 'Analysis'
    config.JobType.psetName    = 'analyzer_MC2016.py'
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
    config.Data.outLFNDirBase  = '/store/user/cmartinp/ttH_Legacy/MC_2016_Oct19/'
    config.Data.publication = True
    config.Data.outputDatasetTag = datasetnames[d]    

    print 'multicrab.py: outLFNDirBase = /store/user/cmartinp/ttH_Legacy/MC_2016_Oct19/'
    #config.Data.publication = True

    config.section_('Site')
    config.Site.storageSite    = 'T2_FR_GRIF_LLR'
    print 'multicrab.py: Submitting Jobs'
    submit(config)
