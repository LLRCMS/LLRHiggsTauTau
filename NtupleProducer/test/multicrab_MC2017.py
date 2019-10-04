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
 'Oct19v1_MC_2017_ttHJetToNonbb', #0
 'Oct19v1_MC_2017_THQ_ctcvcp', #1
 'Oct19v1_MC_2017_THW_ctcvcp', #2
 'Oct19v1_MC_2017_TTH_ctcvcp', #3

 # TTZ
 'Oct19v1_MC_2017_TTZToLLNuNu_M-10', #4
 'Oct19v1_MC_2017_TTZToLLNuNu_M-10_PS', #5
 'Oct19v1_MC_2017_TTZToLL_M-1to10', #6

 # TTW
 'Oct19v1_MC_2017_TTWJetsToLNu', #7
 'Oct19v1_MC_2017_TTWJetsToLNu_PS', #8

 # TTWW
 'Oct19v1_MC_2017_TTWW', #9

 # TT
 'Oct19v1_MC_2017_ST_s-channel', #10
 'Oct19v1_MC_2017_ST_s-channel_PS', #11
 'Oct19v1_MC_2017_ST_t-channel_top', #12
 'Oct19v1_MC_2017_ST_t-channel_top_PS', #13
 'Oct19v1_MC_2017_ST_t-channel_antitop', #14
 'Oct19v1_MC_2017_ST_t-channel_antitop_PS', #15
 'Oct19v1_MC_2017_ST_tW_top', #16
 'Oct19v1_MC_2017_ST_tW_top_PS', #17
 'Oct19v1_MC_2017_ST_tW_antitop', #18
 'Oct19v1_MC_2017_ST_tW_antitop_PS', #19
 'Oct19v1_MC_2017_ST_tWll', #20
 'Oct19v1_MC_2017_TTJets_DiLept', #21
 'Oct19v1_MC_2017_TTJets_SingleLeptFromT', #22
 'Oct19v1_MC_2017_TTJets_SingleLeptFromTbar', #23

 # ggH
 'Oct19v1_MC_2017_GluGluHToTauTau', #24
 'Oct19v1_MC_2017_GluGluHToTauTau_ext1', #25
 'Oct19v1_MC_2017_GluGluHToZZTo4L_ext1', #26
 'Oct19v1_MC_2017_GluGluHToZZTo4L_ext3', #27
 'Oct19v1_MC_2017_GluGluHToZZTo4L_ext4', #28
 'Oct19v1_MC_2017_GluGluHToZZTo2L2Q', #29
 'Oct19v1_MC_2017_GluGluHToWWToLNuQQ', #30
 'Oct19v1_MC_2017_GluGluHToWWTo2L2Nu', #31
 'Oct19v1_MC_2017_GluGluHToMuMu', #32
 'Oct19v1_MC_2017_GluGluHToMuMu_ext1', #33
 'Oct19v1_MC_2017_GluGluHToBB', #34
 'Oct19v1_MC_2017_GluGluHToGG', #35

 # qqH
 'Oct19v1_MC_2017_VBFHToTauTau', #36
 'Oct19v1_MC_2017_VBF_HToZZTo4L_ext2', #37
 'Oct19v1_MC_2017_VBF_HToZZTo4L_ext1', #38
 'Oct19v1_MC_2017_VBF_HToZZTo4L', #39
 'Oct19v1_MC_2017_VBFHToWWToLNuQQ', #40
 'Oct19v1_MC_2017_VBFHToWWTo2L2Nu', #41
 'Oct19v1_MC_2017_VBFHToMuMu', #42
 'Oct19v1_MC_2017_VBFHToBB', #43
 'Oct19v1_MC_2017_VBFHToGG', #44
 'Oct19v1_MC_2017_VBFHToGG_PS', #45

 # Rares
 'Oct19v1_MC_2017_WWW', #46
 'Oct19v1_MC_2017_WWZ', #47
 'Oct19v1_MC_2017_WZZ', #48
 'Oct19v1_MC_2017_ZZZ', #49
 'Oct19v1_MC_2017_WZG', #50
 'Oct19v1_MC_2017_WGToLNuG', #51
 'Oct19v1_MC_2017_ZGToLLG', #52
 'Oct19v1_MC_2017_TGJets', #53
 'Oct19v1_MC_2017_TTGJets', #54
 'Oct19v1_MC_2017_TTGJets_ext1', #55
 'Oct19v1_MC_2017_tZq_ll', #56
 'Oct19v1_MC_2017_WpWpJJ', #57
 'Oct19v1_MC_2017_WWTo2L2Nu_DoubleScattering', #58
 'Oct19v1_MC_2017_TTTT', #59
 'Oct19v1_MC_2017_TTTT_PS', #60

 # VH
 'Oct19v1_MC_2017_VHToNonbb', #61
 'Oct19v1_MC_2017_ZH_HToBB_ZToLL', #62
 'Oct19v1_MC_2017_ZHToTauTau', #63

 # EWK
 'Oct19v1_MC_2017_DYJetsToLL_M-10to50', #64
 'Oct19v1_MC_2017_DYJetsToLL_M-10to50_ext1', #65
 'Oct19v1_MC_2017_DYJetsToLL_M-50', #66
 'Oct19v1_MC_2017_DYJetsToLL_M-50_ext1', #67
 'Oct19v1_MC_2017_WJetsToLNu', #68
 'Oct19v1_MC_2017_WJetsToLNu_ext1', #69
 'Oct19v1_MC_2017_WWTo2L2Nu', #70
 'Oct19v1_MC_2017_WWTo2L2Nu_ext1', #71
 'Oct19v1_MC_2017_WZTo3LNu', #72
 'Oct19v1_MC_2017_ZZTo4L', #73
 'Oct19v1_MC_2017_ZZTo4L_ext1', #74
 'Oct19v1_MC_2017_ZZTo4L_ext2', #75

 # TTVH
 'Oct19v1_MC_2017_TTWH', #76
 'Oct19v1_MC_2017_TTZH', #77

 # HH
 'Oct19v1_MC_2017_GluGluToHHTo2B2VTo2L2Nu_node_SM', #78
 'Oct19v1_MC_2017_GluGluToHHTo2B2VTo2L2Nu_node_2', #79
 'Oct19v1_MC_2017_GluGluToHHTo2B2VTo2L2Nu_node_3', #80
 'Oct19v1_MC_2017_GluGluToHHTo2B2VTo2L2Nu_node_7', #81
 'Oct19v1_MC_2017_GluGluToHHTo2B2VTo2L2Nu_node_9', #82
 'Oct19v1_MC_2017_GluGluToHHTo2B2VTo2L2Nu_node_12', #83
 'Oct19v1_MC_2017_GluGluToHHTo2B2Tau_node_SM', #84
 'Oct19v1_MC_2017_GluGluToHHTo2B2Tau_node_2', #85
 'Oct19v1_MC_2017_GluGluToHHTo2B2Tau_node_3', #86
 'Oct19v1_MC_2017_GluGluToHHTo2B2Tau_node_4', #87
 'Oct19v1_MC_2017_GluGluToHHTo2B2Tau_node_7', #88
 'Oct19v1_MC_2017_GluGluToHHTo2B2Tau_node_9', #89
 'Oct19v1_MC_2017_GluGluToHHTo2B2Tau_node_12', #90
 'Oct19v1_MC_2017_GluGluToHHTo4Tau_node_SM_13', #91
 'Oct19v1_MC_2017_GluGluToHHTo4Tau_node_2', #92
 'Oct19v1_MC_2017_GluGluToHHTo4Tau_node_3', #93
 'Oct19v1_MC_2017_GluGluToHHTo4Tau_node_7', #94
 'Oct19v1_MC_2017_GluGluToHHTo4Tau_node_9', #95
 'Oct19v1_MC_2017_GluGluToHHTo4Tau_node_12', #96
 'Oct19v1_MC_2017_GluGluToHHTo2V2Tau_node_SM', #97
 'Oct19v1_MC_2017_GluGluToHHTo2V2Tau_node_2', #98
 'Oct19v1_MC_2017_GluGluToHHTo2V2Tau_node_3', #99
 'Oct19v1_MC_2017_GluGluToHHTo2V2Tau_node_4', #100
 'Oct19v1_MC_2017_GluGluToHHTo2V2Tau_node_5', #101
 'Oct19v1_MC_2017_GluGluToHHTo2V2Tau_node_6', #102
 'Oct19v1_MC_2017_GluGluToHHTo2V2Tau_node_7', #103
 'Oct19v1_MC_2017_GluGluToHHTo2V2Tau_node_8', #104
 'Oct19v1_MC_2017_GluGluToHHTo2V2Tau_node_9', #105
 'Oct19v1_MC_2017_GluGluToHHTo2V2Tau_node_10', #106
 'Oct19v1_MC_2017_GluGluToHHTo2V2Tau_node_11', #107
 'Oct19v1_MC_2017_GluGluToHHTo2V2Tau_node_12', #108
 'Oct19v1_MC_2017_GluGluToHHTo4V_node_SM', #109
 'Oct19v1_MC_2017_GluGluToHHTo4V_node_2', #110
 'Oct19v1_MC_2017_GluGluToHHTo4V_node_3', #111
 'Oct19v1_MC_2017_GluGluToHHTo4V_node_4', #112
 'Oct19v1_MC_2017_GluGluToHHTo4V_node_5', #113
 'Oct19v1_MC_2017_GluGluToHHTo4V_node_6', #114
 'Oct19v1_MC_2017_GluGluToHHTo4V_node_7', #115
 'Oct19v1_MC_2017_GluGluToHHTo4V_node_8', #116
 'Oct19v1_MC_2017_GluGluToHHTo4V_node_9', #117
 'Oct19v1_MC_2017_GluGluToHHTo4V_node_10', #118
 'Oct19v1_MC_2017_GluGluToHHTo4V_node_11', #119
 'Oct19v1_MC_2017_GluGluToHHTo4V_node_12' #120
 
    ]

 datasetinputs = [

 # Signal
 '/ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM', #0
 '/THQ_ctcvcp_4f_Hincl_13TeV_madgraph_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #1
 '/THW_ctcvcp_5f_Hincl_13TeV_madgraph_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #2
 '/TTH_4f_ctcvcp_TuneCP5_13TeV_madgraph_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #3

 # TTZ
 '/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM', #4
 '/TTZToLLNuNu_M-10_TuneCP5_PSweights_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #5
 '/TTZToLL_M-1to10_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM', #6

 # TTW
 '/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM', #7
 '/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM', #8

 # TTWW
 '/TTWW_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM', #9

 # TT
 '/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM', #10
 '/ST_s-channel_4f_leptonDecays_TuneCP5_PSweights_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM', #11
 '/ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM', #12
 '/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #13
 '/ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM', #14
 '/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #15
 '/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM', #16
 '/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM', #17
 '/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM', #18
 '/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM', #19
 '/ST_tWll_5f_LO_TuneCP5_PSweights_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM', #20
 '/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #21
 '/TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #22
 '/TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #23

 # ggH
 '/GluGluHToTauTau_M125_13TeV_powheg_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v2/MINIAODSIM', #24
 '/GluGluHToTauTau_M125_13TeV_powheg_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM', #25
 '/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM', #26
 '/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext3-v1/MINIAODSIM', #27
 '/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext4-v1/MINIAODSIM', #28
 '/GluGluHToZZTo2L2Q_M125_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM', #29
 '/GluGluHToWWToLNuQQ_M125_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen710_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM', #30
 '/GluGluHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #31
 '/GluGluHToMuMu_M-125_TuneCP5_PSweights_13TeV_powheg_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #32
 '/GluGluHToMuMu_M-125_TuneCP5_PSweights_13TeV_powheg_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM', #33
 '/GluGluHToBB_M125_13TeV_amcatnloFXFX_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM', #34
 '/GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM', #35

 # qqH
 '/VBFHToTauTau_M125_13TeV_powheg_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM', #36
 '/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext2-v2/MINIAODSIM', #37
 '/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM', #38
 '/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #39
 '/VBFHToWWToLNuQQ_M125_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen710_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #40
 '/VBFHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #41
 '/VBFHToMuMu_M-125_TuneCP5_PSweights_13TeV_powheg_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #42
 '/VBFHToBB_M-125_13TeV_powheg_pythia8_weightfix/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #43
 '/VBFHToGG_M125_13TeV_amcatnlo_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #44
 '/VBFHToGG_M125_13TeV_amcatnlo_pythia8_PSWeights/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM', #45

 # Rares
 '/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM', #46
 '/WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM', #47
 '/WZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #48
 '/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM', #49
 '/WZG_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #50
 '/WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #51
 '/ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/MINIAODSIM', #52
 '/TGJets_leptonDecays_TuneCP5_PSweights_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #53
 '/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM', #54
 '/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2/MINIAODSIM', #55
 '/tZq_ll_4f_ckm_NLO_TuneCP5_PSweights_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v2/MINIAODSIM', #56
 '/WpWpJJ_EWK-QCD_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #57
 '/WWTo2L2Nu_DoubleScattering_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM', #58
 '/TTTT_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM', #59
 '/TTTT_TuneCP5_PSweights_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #60

 # VH
 '/VHToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM', #61
 '/ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #62
 '/ZHToTauTau_M125_13TeV_powheg_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #63

 # EWK
 '/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #64
 '/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2/MINIAODSIM', #65
 '/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM', #66
 '/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM', #67
 '/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/MINIAODSIM', #68
 '/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2/MINIAODSIM', #69
 '/WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #70
 '/WWTo2L2Nu_NNPDF31_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM', #71
 '/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM', #72
 '/ZZTo4L_13TeV_powheg_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM', #73
 '/ZZTo4L_13TeV_powheg_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM', #74
 '/ZZTo4L_13TeV_powheg_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext2-v1/MINIAODSIM', #75

 # TTVH
 '/TTWH_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM', #76
 '/TTZH_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM', #77

 # HH
 '/GluGluToHHTo2B2VTo2L2Nu_node_SM_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #78
 '/GluGluToHHTo2B2VTo2L2Nu_node_2_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #79
 '/GluGluToHHTo2B2VTo2L2Nu_node_3_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #80
 '/GluGluToHHTo2B2VTo2L2Nu_node_7_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #81
 '/GluGluToHHTo2B2VTo2L2Nu_node_9_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #82
 '/GluGluToHHTo2B2VTo2L2Nu_node_12_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #83
 '/GluGluToHHTo2B2Tau_node_SM_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #84
 '/GluGluToHHTo2B2Tau_node_2_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #85
 '/GluGluToHHTo2B2Tau_node_3_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #86
 '/GluGluToHHTo2B2Tau_node_4_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #87
 '/GluGluToHHTo2B2Tau_node_7_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #88
 '/GluGluToHHTo2B2Tau_node_9_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSI', #89
 '/GluGluToHHTo2B2Tau_node_12_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #90
 '/GluGluToHHTo4Tau_node_SM_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #91
 '/GluGluToHHTo4Tau_node_2_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #92
 '/GluGluToHHTo4Tau_node_3_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #93
 '/GluGluToHHTo4Tau_node_7_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #94
 '/GluGluToHHTo4Tau_node_9_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #95
 '/GluGluToHHTo4Tau_node_12_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #96
 '/GluGluToHHTo2V2Tau_node_SM_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #97
 '/GluGluToHHTo2V2Tau_node_2_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #98
 '/GluGluToHHTo2V2Tau_node_3_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #99
 '/GluGluToHHTo2V2Tau_node_4_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #100
 '/GluGluToHHTo2V2Tau_node_5_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #101
 '/GluGluToHHTo2V2Tau_node_6_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #102
 '/GluGluToHHTo2V2Tau_node_7_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #103
 '/GluGluToHHTo2V2Tau_node_8_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #104
 '/GluGluToHHTo2V2Tau_node_9_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #105
 '/GluGluToHHTo2V2Tau_node_10_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #106
 '/GluGluToHHTo2V2Tau_node_11_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #107
 '/GluGluToHHTo2V2Tau_node_12_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #108
 '/GluGluToHHTo4V_node_SM_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #109
 '/GluGluToHHTo4V_node_2_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #110
 '/GluGluToHHTo4V_node_3_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #111
 '/GluGluToHHTo4V_node_4_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #112
 '/GluGluToHHTo4V_node_5_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #113
 '/GluGluToHHTo4V_node_6_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #114
 '/GluGluToHHTo4V_node_7_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #115
 '/GluGluToHHTo4V_node_8_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #116
 '/GluGluToHHTo4V_node_9_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #117
 '/GluGluToHHTo4V_node_10_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #118
 '/GluGluToHHTo4V_node_11_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM', #119
 '/GluGluToHHTo4V_node_12_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM' #120

    ]


for d in range(0,len(datasetnames)):

    print 'multicrab.py: Running datasetname: ', datasetnames[d]

    config.section_('General')
    config.General.requestName = datasetnames[d]
    config.General.workArea    = 'crab3_Oct19'
    config.General.transferLogs = True

    config.section_('JobType')
    config.JobType.pluginName  = 'Analysis'
    config.JobType.psetName    = 'analyzer_MC2017.py'
    config.JobType.sendExternalFolder = True
    config.JobType.maxMemoryMB = 2000 # Default == 2Gb : maximum guaranteed to run on all sites
    config.JobType.inputFiles=['JECUncertaintySources']

    config.section_("Data")
    config.Data.inputDataset   = datasetinputs[d]
    config.Data.inputDBS = 'global'
    config.Data.inputDBS       = 'global'
    #config.Data.splitting      = 'FileBased'
    config.Data.splitting      = 'Automatic'
    #config.Data.totalUnits     = 40000 #With 'FileBased' splitting tells how many files to analyse
    config.Data.unitsPerJob    = 180
    config.Data.outLFNDirBase  = '/store/user/cmartinp/ttH_Legacy/MC_2017_Oct19v1/'
    config.Data.publication    = True
    config.Data.outputDatasetTag = datasetnames[d]    

    print 'multicrab.py: outLFNDirBase = /store/user/cmartinp/ttH_Legacy/MC_2017_Oct19v1/'
    #config.Data.publication = True

    config.section_('Site')
    config.Site.storageSite    = 'T2_FR_GRIF_LLR' #T2_FR_GRIF_IRFU
    print 'multicrab.py: Submitting Jobs'
    submit(config)
