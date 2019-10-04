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

 # Signal
 'Oct19v1_MC_2018_ttHJetToNonbb', #0
 'Oct19v1_MC_2018_THQ_ctcvcp', #1
 'Oct19v1_MC_2018_THW_ctcvcp', #2
 'Oct19v1_MC_2018_TTH_ctcvcp', #3

 # TTZ
 'Oct19v1_MC_2018_TTZToLLNuNu_M-10', #4
 'Oct19v1_MC_2018_TTZToLL_M-1to10', #5

 # TTW
 'Oct19v1_MC_2018_TTWJetsToLNu', #6

 # TTWW
 'Oct19v1_MC_2018_TTWW', #7

 # TT
 'Oct19v1_MC_2018_ST_s-channel', #8
 'Oct19v1_MC_2018_ST_t-channel_top', #9
 'Oct19v1_MC_2018_ST_t-channel_antitop', #10
 'Oct19v1_MC_2018_ST_tW_top', #11
 'Oct19v1_MC_2018_ST_tW_antitop', #12
 'Oct19v1_MC_2018_ST_tWll', #13
 'Oct19v1_MC_2018_TTJets_DiLept', #14
 'Oct19v1_MC_2018_TTJets_SingleLeptFromT', #15
 'Oct19v1_MC_2018_TTJets_SingleLeptFromTbar', #16

 # ggH
 'Oct19v1_MC_2018_GluGluHToTauTau', #17
 'Oct19v1_MC_2018_GluGluHToZZTo4L', #18
 'Oct19v1_MC_2018_GluGluHToZZTo2L2Q', #19
 'Oct19v1_MC_2018_GluGluHToWWToLNuQQ', #20
 'Oct19v1_MC_2018_GluGluHToWWTo2L2Nu', #21
 'Oct19v1_MC_2018_GluGluHToMuMu', #22
 'Oct19v1_MC_2018_GluGluHToMuMu_ext1',  #23
 'Oct19v1_MC_2018_GluGluHToBB', #24
 'Oct19v1_MC_2018_GluGluHToGG', #25

 # qqH
 'Oct19v1_MC_2018_VBFHToTauTau', #26
 'Oct19v1_MC_2018_VBF_HToZZTo4L', #27
 'Oct19v1_MC_2018_VBFHToWWToLNuQQ', #28
 'Oct19v1_MC_2018_VBFHToWWTo2L2Nu', #29
 'Oct19v1_MC_2018_VBFHToMuMu', #30
 'Oct19v1_MC_2018_VBFHToBB', #31
 'Oct19v1_MC_2018_VBFHToGG', #32

 # Rares
 'Oct19v1_MC_2018_WWW', #33
 'Oct19v1_MC_2018_WWZ', #34
 'Oct19v1_MC_2018_WZZ', #35
 'Oct19v1_MC_2018_ZZZ', #36
 'Oct19v1_MC_2018_WZG', #37
 'Oct19v1_MC_2018_WGToLNuG', #38
 'Oct19v1_MC_2018_ZGToLLG', #39
 'Oct19v1_MC_2018_TGJets', #40
 'Oct19v1_MC_2018_TTGJets', #41
 'Oct19v1_MC_2018_tZq_ll', #42
 'Oct19v1_MC_2018_WpWpJJ', #43
 'Oct19v1_MC_2018_WWTo2L2Nu_DoubleScattering', #44
 'Oct19v1_MC_2018_TTTT', #45

 # VH
 'Oct19v1_MC_2018_VHToNonbb', #46
 'Oct19v1_MC_2018_ZH_HToBB_ZToLL', #47
 'Oct19v1_MC_2018_ZH_HToBB_ZToLL_ext1', #48
 'Oct19v1_MC_2018_ZHToTauTau', #49

 # EWK
 'Oct19v1_MC_2018_DYJetsToLL_M-10to50', #50
 'Oct19v1_MC_2018_DYJetsToLL_M-50', #51
 'Oct19v1_MC_2018_DYJetsToLL_M-50_ext2', #52
 'Oct19v1_MC_2018_WJetsToLNu', #53
 'Oct19v1_MC_2018_WWTo2L2Nu', #54
 'Oct19v1_MC_2018_WZTo3LNu', #55
 'Oct19v1_MC_2018_ZZTo4L', #56
 'Oct19v1_MC_2018_ZZTo4L_ext2', #57
 
 # TTVH
 'Oct19v1_MC_2018_TTWH', #58
 'Oct19v1_MC_2018_TTZH', #59

 # HH
 'GluGluToHHTo2B2Tau_node_SM', #60
 'GluGluToHHTo2B2Tau_node_2', #61
 'GluGluToHHTo2B2Tau_node_3', #62
 'GluGluToHHTo2B2Tau_node_4', #63
 'GluGluToHHTo2B2Tau_node_5', #64
 'GluGluToHHTo2B2Tau_node_6', #65
 'GluGluToHHTo2B2Tau_node_7', #66
 'GluGluToHHTo2B2Tau_node_8', #67
 'GluGluToHHTo2B2Tau_node_9', #68
 'GluGluToHHTo2B2Tau_node_10', #69
 'GluGluToHHTo2B2Tau_node_11', #70
 'GluGluToHHTo2B2Tau_node_12', #71


    ]

 datasetinputs = [

 # Signal
 '/ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #0
 '/THQ_ctcvcp_4f_Hincl_13TeV_madgraph_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #1
 '/THW_ctcvcp_5f_Hincl_13TeV_madgraph_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #2
 '/TTH_4f_ctcvcp_TuneCP5_13TeV_madgraph_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #3

 # TTZ
 '/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM', #4
 '/TTZToLL_M-1to10_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #5

 # TTW
 '/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM', #6

 # TTWW
 '/TTWW_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v1/MINIAODSIM', #7

 # TT
 '/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v4/MINIAODSIM', #8
 '/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #9
 '/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #10
 '/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM', #11
 '/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM', #12
 '/ST_tWll_5f_LO_TuneCP5_PSweights_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM', #13
 '/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #14
 '/TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #15
 '/TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #16

 # ggH
 '/GluGluHToTauTau_M125_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM', #17
 '/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM', #18
 '/GluGluHToZZTo2L2Q_M125_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #19
 '/GluGluHToWWToLNuQQ_M125_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #20
 '/GluGluHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #21
 '/GluGluHToMuMu_M-125_TuneCP5_PSweights_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM', #22
 '/GluGluHToMuMu_M-125_TuneCP5_PSweights_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM', #23
 '/GluGluHToBB_M125_13TeV_amcatnloFXFX_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #24
 '/GluGluHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #25

 # qqH
 '/VBFHToTauTau_M125_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM', #26
 '/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM', #27
 '/VBFHToWWToLNuQQ_M125_13TeV_powheg_JHUGen_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #28
 '/VBFHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #29
 '/VBFHToMuMu_M-125_TuneCP5_PSweights_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #30
 '/VBFHToBB_M-125_13TeV_powheg_pythia8_weightfix/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #31
 '/VBFHToGG_M125_13TeV_amcatnlo_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #32

 # Rares
 '/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM', #33
 '/WWZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM', #34
 '/WZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM', #35
 '/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM', #36
 '/WZG_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #37
 '/WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #38
 '/ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM', #39
 '/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM', #40
 '/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #41
 '/tZq_ll_4f_ckm_NLO_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM', #42
 '/WpWpJJ_EWK-QCD_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #43
 '/WWTo2L2Nu_DoubleScattering_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #44
 '/TTTT_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM', #45

 # VH
 '/VHToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM', #46
 '/ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM', #47
 '/ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM', #48
 '/ZHToTauTau_M125_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM', #49

 # EWK
 '/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM', #50
 '/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #51
 '/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v1/MINIAODSIM', #52
 '/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM', #53
 '/WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #54
 '/WZTo3LNu_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM', #55
 '/ZZTo4L_TuneCP5_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM', #56
 '/ZZTo4L_TuneCP5_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v2/MINIAODSIM', #57

 # TTVH
 '/TTWH_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM', #58
 '/TTZH_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM', #59

 # HH
 '/GluGluToHHTo2B2Tau_node_SM_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #60
 '/GluGluToHHTo2B2Tau_node_2_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #61
 '/GluGluToHHTo2B2Tau_node_3_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #62
 '/GluGluToHHTo2B2Tau_node_4_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #63
 '/GluGluToHHTo2B2Tau_node_5_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #64
 '/GluGluToHHTo2B2Tau_node_6_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #65
 '/GluGluToHHTo2B2Tau_node_7_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #66
 '/GluGluToHHTo2B2Tau_node_8_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #67
 '/GluGluToHHTo2B2Tau_node_9_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #68
 '/GluGluToHHTo2B2Tau_node_10_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #69
 '/GluGluToHHTo2B2Tau_node_11_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #70
 '/GluGluToHHTo2B2Tau_node_12_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', #71

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
    config.JobType.sendExternalFolder = True  # Default == 2Gb : maximum guaranteed to run on all sites
    config.JobType.maxMemoryMB = 2000
    config.JobType.inputFiles=['JECUncertaintySources']

    config.section_("Data")
    config.Data.inputDataset   = datasetinputs[d]
    config.Data.inputDBS = 'global'
    config.Data.inputDBS       = 'global'
    #config.Data.splitting      = 'FileBased'
    config.Data.splitting      = 'Automatic'
    #config.Data.totalUnits     = 40000 #With 'FileBased' splitting tells how many files to analyse
    config.Data.unitsPerJob    = 180
    config.Data.outLFNDirBase  = '/store/user/cmartinp/ttH_Legacy/MC_2018_Oct19v1/'
    config.Data.publication    = True
    config.Data.outputDatasetTag = datasetnames[d]    

    print 'multicrab.py: outLFNDirBase = /store/user/cmartinp/ttH_Legacy/MC_2018_Oct19v1/'
    #config.Data.publication = True

    config.section_('Site')
    config.Site.storageSite    = 'T2_FR_GRIF_LLR' #T2_FR_GRIF_IRFU
    print 'multicrab.py: Submitting Jobs'
    submit(config)
