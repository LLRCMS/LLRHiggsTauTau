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
 'Oct19v1_MC_2016_ttHJetToNonbb', #0
 'Oct19v1_MC_2016_THQ_ctcvcp', #1
 'Oct19v1_MC_2016_THW_ctcvcp', #2
 'Oct19v1_MC_2016_ttH_ctcvcp', #3

 # TTZ
 'Oct19v1_MC_2016_TTZToLLNuNu_M-10_ext2', #4
 'Oct19v1_MC_2016_TTZToLLNuNu_M-10_ext3', #5
 'Oct19v1_MC_2016_TTZToLL_M-1to10', #6

 # TTW
 'Oct19v1_MC_2016_TTWJetsToLNu', #7

 # TTWW
 'Oct19v1_MC_2016_TTWW', #8

 # TT
 'Oct19v1_MC_2016_ST_s-channel', #9
 'Oct19v1_MC_2016_ST_s-channel_PS', #10
 'Oct19v1_MC_2016_ST_t-channel_top', #11
 'Oct19v1_MC_2016_ST_t-channel_antitop', #12
 'Oct19v1_MC_2016_ST_t-channel_antitop_PS', #13
 'Oct19v1_MC_2016_ST_tW_top', #14
 'Oct19v1_MC_2016_ST_tW_antitop', #15
 'Oct19v1_MC_2016_ST_tWll', #16
 'Oct19v1_MC_2016_TTJets_DiLept_ext1', #17
 'Oct19v1_MC_2016_TTJets_DiLept', #18
 'Oct19v1_MC_2016_TTJets_SingleLeptFromT_ext1', #19
 'Oct19v1_MC_2016_TTJets_SingleLeptFromT', #20
 'Oct19v1_MC_2016_TTJets_SingleLeptFromTbar', #21
 'Oct19v1_MC_2016_TTJets_SingleLeptFromTbar_ext1', #22

 # ggH
 'Oct19v1_MC_2016_GluGluHToTauTau', #23
 'Oct19v1_MC_2016_GluGluHToZZTo4L', #24
 'Oct19v1_MC_2016_GluGluHToWWToLNuQQ', #25
 'Oct19v1_MC_2016_GluGluHToWWTo2L2Nu', #26
 'Oct19v1_MC_2016_GluGluHToMuMu', #27
 'Oct19v1_MC_2016_GluGluHToBB', #28
 'Oct19v1_MC_2016_GluGluHToBB_ext1', #29
 'Oct19v1_MC_2016_GluGluHToGG', #30

 # qqH
 'Oct19v1_MC_2016_VBFHToTauTau', #31
 'Oct19v1_MC_2016_VBF_HToZZTo4L', #32
 'Oct19v1_MC_2016_VBFHToWWToLNuQQ', #33
 'Oct19v1_MC_2016_VBFHToWWTo2L2Nu', #34
 'Oct19v1_MC_2016_VBFHToMuMu', #35
 'Oct19v1_MC_2016_VBFHToBB', #36
 'Oct19v1_MC_2016_VBFHToBB_ext1', #37
 'Oct19v1_MC_2016_VBFHToGG_ext1', #38
 'Oct19v1_MC_2016_VBFHToGG_ext2', #39

 # Rares
 'Oct19v1_MC_2016_WWW', #40
 'Oct19v1_MC_2016_WWZ', #41
 'Oct19v1_MC_2016_WZZ', #42
 'Oct19v1_MC_2016_ZZZ', #43
 'Oct19v1_MC_2016_WZG', #44
 'Oct19v1_MC_2016_WGToLNuG_ext1', #45
 'Oct19v1_MC_2016_WGToLNuG_ext2', #46
 'Oct19v1_MC_2016_WGToLNuG_ext3', #47
 'Oct19v1_MC_2016_ZGTo2LG', #48
 'Oct19v1_MC_2016_TGJets_leptonDecays', #49
 'Oct19v1_MC_2016_TTGJets', #50
 'Oct19v1_MC_2016_TTGJets_ext1', #51
 'Oct19v1_MC_2016_tZq_ll', #52
 'Oct19v1_MC_2016_tZq_ll_PS', #53
 'Oct19v1_MC_2016_WpWpJJ', #54
 'Oct19v1_MC_2016_WWTo2L2Nu_DoubleScattering', #55
 'Oct19v1_MC_2016_TTTT', #56

 # VH
 'Oct19v1_MC_2016_VHToNonbb', #57
 'Oct19v1_MC_2016_ZH_HToBB_ZToLL', #58
 'Oct19v1_MC_2016_ZHToTauTau', #59

 # EWK
 'Oct19v1_MC_2016_DYJetsToLL_M-10to50', #60
 'Oct19v1_MC_2016_DYJetsToLL_M-50', #61
 'Oct19v1_MC_2016_WJetsToLNu_ext2', #62
 'Oct19v1_MC_2016_WJetsToLNu', #63
 'Oct19v1_MC_2016_WWTo2L2Nu', #64
 'Oct19v1_MC_2016_WZTo3LNu', #65
 'Oct19v1_MC_2016_ZZTo4L', #66

 # TTVH
 'Oct19v1_MC_2016_TTWH', #67
 'Oct19v1_MC_2016_TTZH', #68

 # HH
 'Oct19v1_MC_2016_GluGluToHHTo2B2VTo2L2Nu_node_SM', #69
 'Oct19v1_MC_2016_GluGluToHHTo2B2VTo2L2Nu_node_box', #70
 'Oct19v1_MC_2016_GluGluToHHTo2B2VTo2L2Nu_node_1', #71
 'Oct19v1_MC_2016_GluGluToHHTo2B2VTo2L2Nu_node_2', #72
 'Oct19v1_MC_2016_GluGluToHHTo2B2VTo2L2Nu_node_3', #73
 'Oct19v1_MC_2016_GluGluToHHTo2B2VTo2L2Nu_node_4', #74
 'Oct19v1_MC_2016_GluGluToHHTo2B2VTo2L2Nu_node_5', #75
 'Oct19v1_MC_2016_GluGluToHHTo2B2VTo2L2Nu_node_6', #76
 'Oct19v1_MC_2016_GluGluToHHTo2B2VTo2L2Nu_node_7', #77
 'Oct19v1_MC_2016_GluGluToHHTo2B2VTo2L2Nu_node_8', #78
 'Oct19v1_MC_2016_GluGluToHHTo2B2VTo2L2Nu_node_9', #79
 'Oct19v1_MC_2016_GluGluToHHTo2B2VTo2L2Nu_node_10', #80
 'Oct19v1_MC_2016_GluGluToHHTo2B2VTo2L2Nu_node_11', #81
 'Oct19v1_MC_2016_GluGluToHHTo2B2VTo2L2Nu_node_12', #82
 'Oct19v1_MC_2016_GluGluToHHTo2B2Tau_node_SM', #83
 'Oct19v1_MC_2016_GluGluToHHTo2B2Tau_node_box', #84
 'Oct19v1_MC_2016_GluGluToHHTo2B2Tau_node_2', #85
 'Oct19v1_MC_2016_GluGluToHHTo2B2Tau_node_9', #86
 'Oct19v1_MC_2016_GluGluToHHTo2B2Tau_node_10', #87
 'Oct19v1_MC_2016_GluGluToHHTo2B2Tau_node_11', #88
 'Oct19v1_MC_2016_GluGluToHHTo2B2Tau_node_12', #89
 'Oct19v1_MC_2016_GluGluToHHTo2B2Tau_node_13', #90
 'Oct19v1_MC_2016_GluGluToHHTo4Tau_node_SM', #91
 'Oct19v1_MC_2016_GluGluToHHTo4Tau_node_box', #92
 'Oct19v1_MC_2016_GluGluToHHTo4Tau_node_2', #93
 'Oct19v1_MC_2016_GluGluToHHTo4Tau_node_3', #94
 'Oct19v1_MC_2016_GluGluToHHTo4Tau_node_4', #95
 'Oct19v1_MC_2016_GluGluToHHTo4Tau_node_5', #96
 'Oct19v1_MC_2016_GluGluToHHTo4Tau_node_6', #97
 'Oct19v1_MC_2016_GluGluToHHTo4Tau_node_7', #98
 'Oct19v1_MC_2016_GluGluToHHTo4Tau_node_8', #99
 'Oct19v1_MC_2016_GluGluToHHTo4Tau_node_9', #100
 'Oct19v1_MC_2016_GluGluToHHTo4Tau_node_10', #101
 'Oct19v1_MC_2016_GluGluToHHTo4Tau_node_11', #102
 'Oct19v1_MC_2016_GluGluToHHTo4Tau_node_12' #103

    ]

 datasetinputs = [

 # Signal -> datacards
 '/ttHJetToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8_mWCutfix/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM', #0
 '/THQ_ctcvcp_HIncl_M125_TuneCP5_13TeV-madgraph-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', #1
 '/THW_ctcvcp_HIncl_M125_TuneCP5_13TeV-madgraph-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', #2
 '/ttH_4f_ctcvcp_TuneCP5_13TeV_madgraph_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', #3
 
 # TTZ
 '/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/MINIAODSIM', #4
 '/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext3-v1/MINIAODSIM', #5
 '/TTZToLL_M-1to10_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', #6

 # TTW
 '/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/MINIAODSIM', #7

 # TTWW
 '/TTWW_TuneCUETP8M2T4_13TeV-madgraph-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/MINIAODSIM', #8

 # TT
 '/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', #9
 '/ST_s-channel_4f_leptonDecays_13TeV_PSweights-amcatnlo-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', #10
 '/ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', #11
 '/ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', #12
 '/ST_t-channel_antitop_4f_inclusiveDecays_13TeV_PSweights-powhegV2-madspin/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', #13
 '/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/MINIAODSIM', #14
 '/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/MINIAODSIM', #15
 '/ST_tWll_5f_LO_13TeV-MadGraph-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', #16
 '/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM', #17
 '/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #18
 '/TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM', #19
 '/TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #20
 '/TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #21
 '/TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM', #22

  # ggH
 '/GluGluHToTauTau_M125_13TeV_powheg_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v3/MINIAODSIM', #23
 '/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', #24
 '/GluGluHToWWToLNuQQ_M125_13TeV_powheg_JHUGenV628_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/MINIAODSIM', #25
 '/GluGluHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #26
 '/GluGluHToMuMu_M-125_TuneCP5_PSweights_13TeV_powheg_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', #27
 '/GluGluHToBB_M125_13TeV_amcatnloFXFX_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #28
 '/GluGluHToBB_M125_13TeV_amcatnloFXFX_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/MINIAODSIM', #29
 '/GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2/MINIAODSIM', #30

 # qqH
 '/VBFHToTauTau_M125_13TeV_powheg_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #31
 '/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #32
 '/VBFHToWWToLNuQQ_M125_13TeV_powheg_JHUGenV628_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/MINIAODSIM', #33
 '/VBFHToWWTo2L2Nu_M125_13TeV_powheg_JHUgenv628_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #34
 '/VBFHToMuMu_M-125_TuneCP5_PSweights_13TeV_powheg_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', #35
 '/VBFHToBB_M-125_13TeV_powheg_pythia8_weightfix/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', #36
 '/VBFHToBB_M-125_13TeV_powheg_pythia8_weightfix/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/MINIAODSIM', #37
 '/VBFHToGG_M125_13TeV_amcatnlo_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM', #38
 '/VBFHToGG_M125_13TeV_amcatnlo_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2/MINIAODSIM',  #39

 # Rares
 '/WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #40
 '/WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #41
 '/WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #42
 '/ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #43
 '/WZG_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', #44
 '/WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/MINIAODSIM', #45
 '/WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/MINIAODSIM', #46
 '/WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext3-v1/MINIAODSIM', #47
 '/ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/MINIAODSIM', #48
 '/TGJets_leptonDecays_13TeV_amcatnlo_madspin_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #49
 '/TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #50
 '/TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM', #51
 '/tZq_ll_4f_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/MINIAODSIM', #52
 '/tZq_ll_4f_PSweights_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', #53
 '/WpWpJJ_EWK-QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', #54
 '/WWTo2L2Nu_DoubleScattering_13TeV-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #55
 '/TTTT_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #56

 # VH
 '/VHToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #57
 '/ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/MINIAODSIM', #58
 '/ZHToTauTau_M125_13TeV_powheg_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #59

 # EWK
 '/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #60
 '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/MINIAODSIM', #61
 '/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2/MINIAODSIM', #62
 '/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #63
 '/WWTo2L2Nu_13TeV-powheg/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #64
 '/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/MINIAODSIM', #65
 '/ZZTo4L_13TeV_powheg_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', #66

 # TTVH
 '/TTWH_TuneCUETP8M2T4_13TeV-madgraph-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/MINIAODSIM', #67
 '/TTZH_TuneCUETP8M2T4_13TeV-madgraph-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/MINIAODSIM', #68

 # HH
 '/GluGluToHHTo2B2VTo2L2Nu_node_SM_13TeV-madgraph-v2/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #69
 '/GluGluToHHTo2B2VTo2L2Nu_node_box_13TeV-madgraph-v2/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #70
 '/GluGluToHHTo2B2VTo2L2Nu_node_1_13TeV-madgraph-v2/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #71
 '/GluGluToHHTo2B2VTo2L2Nu_node_2_13TeV-madgraph-v2/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v3/MINIAODSIM', #72
 '/GluGluToHHTo2B2VTo2L2Nu_node_3_13TeV-madgraph-v2/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #73
 '/GluGluToHHTo2B2VTo2L2Nu_node_4_13TeV-madgraph-v2/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #74
 '/GluGluToHHTo2B2VTo2L2Nu_node_5_13TeV-madgraph-v2/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #75
 '/GluGluToHHTo2B2VTo2L2Nu_node_6_13TeV-madgraph-v2/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #76
 '/GluGluToHHTo2B2VTo2L2Nu_node_7_13TeV-madgraph-v2/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #77
 '/GluGluToHHTo2B2VTo2L2Nu_node_8_13TeV-madgraph-v2/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #78
 '/GluGluToHHTo2B2VTo2L2Nu_node_9_13TeV-madgraph-v2/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #79
 '/GluGluToHHTo2B2VTo2L2Nu_node_10_13TeV-madgraph-v2/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #80
 '/GluGluToHHTo2B2VTo2L2Nu_node_11_13TeV-madgraph-v2/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #81
 '/GluGluToHHTo2B2VTo2L2Nu_node_12_13TeV-madgraph-v2/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #82
 '/GluGluToHHTo2B2Tau_node_SM_13TeV-madgraph/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #83
 '/GluGluToHHTo2B2Tau_node_box_13TeV-madgraph/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #84
 '/GluGluToHHTo2B2Tau_node_2_13TeV-madgraph/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #85
 '/GluGluToHHTo2B2Tau_node_9_13TeV-madgraph/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #86
 '/GluGluToHHTo2B2Tau_node_10_13TeV-madgraph/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #87
 '/GluGluToHHTo2B2Tau_node_11_13TeV-madgraph/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #88
 '/GluGluToHHTo2B2Tau_node_12_13TeV-madgraph/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #89
 '/GluGluToHHTo2B2Tau_node_13_13TeV-madgraph/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', #90
 '/GluGluToHHTo4Tau_node_SM_13TeV-madgraph/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #91
 '/GluGluToHHTo4Tau_node_box_13TeV-madgraph/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #92
 '/GluGluToHHTo4Tau_node_2_13TeV-madgraph/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #93
 '/GluGluToHHTo4Tau_node_3_13TeV-madgraph/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #94
 '/GluGluToHHTo4Tau_node_4_13TeV-madgraph/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #95
 '/GluGluToHHTo4Tau_node_5_13TeV-madgraph/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #96
 '/GluGluToHHTo4Tau_node_6_13TeV-madgraph/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #97
 '/GluGluToHHTo4Tau_node_7_13TeV-madgraph/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #98
 '/GluGluToHHTo4Tau_node_8_13TeV-madgraph/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #99
 '/GluGluToHHTo4Tau_node_9_13TeV-madgraph/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #100
 '/GluGluToHHTo4Tau_node_10_13TeV-madgraph/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #101
 '/GluGluToHHTo4Tau_node_11_13TeV-madgraph/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', #102
 '/GluGluToHHTo4Tau_node_12_13TeV-madgraph/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM' #103

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
    config.JobType.maxMemoryMB = 2000
    config.JobType.inputFiles=['JECUncertaintySources']

    config.section_("Data")
    config.Data.inputDataset   = datasetinputs[d]
    config.Data.inputDBS = 'global'
    config.JobType.maxMemoryMB = 2000 # Default == 2Gb : maximum guaranteed to run on all sites
    config.Data.inputDBS       = 'global'
    #config.Data.splitting      = 'FileBased'
    config.Data.splitting      = 'Automatic'
    #config.Data.totalUnits     = 40000 #With 'FileBased' splitting tells how many files to analyse
    config.Data.unitsPerJob    = 180
    config.Data.unitsPerJob    = 1
    config.Data.outLFNDirBase  = '/store/user/cmartinp/ttH_Legacy/MC_2016_Oct19v1/'
    config.Data.publication = True
    config.Data.outputDatasetTag = datasetnames[d]    

    print 'multicrab.py: outLFNDirBase = /store/user/cmartinp/ttH_Legacy/MC_2016_Oct19v1/'
    #config.Data.publication = True

    config.section_('Site')
    config.Site.storageSite    = 'T2_FR_GRIF_LLR' #T2_FR_GRIF_IRFU
    print 'multicrab.py: Submitting Jobs'
    submit(config)
