import os
import sys
import re

# NB: you need to source crab environment before lauching this script:
# source /cvmfs/cms.cern.ch/crab3/crab.sh

###################################################################
#### Parameters to be changed for each production


#PROCESS = ["HHBACKGROUNDS"]
#tag = "llrNt_NoSVFit_bkg_27Ago2015"
#datasetsFile = "datasets.txt"

#PROCESS = ["2015RUNBDATA", "2015RUNCDATA"] # 50 ns
#tag = "Data50ns_SVFit_6Ott2015"
#datasetsFile = "datasets.txt"

#PROCESS = ["2015RUNCDATA", "2015RUNDDATA"] # 25 ns
#tag = "Data25ns_SVFit_6Ott2015"
#datasetsFile = "datasets.txt"

#PROCESS = ["2015DATA25NSRESUBMISSION"]
#tag = "Data25ns_SVFit_6Ott2015_resub"
#datasetsFile = "datasets.txt"

#PROCESS = ["2015DATARUND27OTT"]
#tag = "Data25ns_noSVFit_npvFix_30Ott2015_newJson_trigFix2"
#datasetsFile = "datasets.txt"

#PROCESS = ["2015DATAPROMPTRECOONLY"]
#tag = "Data25ns_noSVFit_lumiDiff13Nov2015"
#datasetsFile = "datasets.txt"

#PROCESS = ["MINIAODV2"]
#tag = "MC_NoSVFit_MiniAODV2_13Nov2015"
#datasetsFile = "datasets.txt"

# PROCESS = ["RESUBMINIAODV2"]
# tag = "MC_NoSVFit_MiniAODV2_13Nov2015_DYresub"
# datasetsFile = "datasets.txt"

# PROCESS = ["MINIV2SVFIT"]
# tag = "MC_SVFit_MiniAODV2_22Nov2015"
# datasetsFile = "datasets.txt"

# PROCESS = ["2015DATARUND27OTT"]
# tag = "Data_SVFit_MiniAODV2_22Nov2015_SVFix"
# datasetsFile = "datasets.txt"

# PROCESS = ["MINIV2SVFITDY"]
# tag = "MC_SVFit_MiniAODV2_22Nov2015_DYResub"
# datasetsFile = "datasets.txt"

# PROCESS = ["MINIV2SVFITPLUS"]
# tag = "MC_SVFit_MiniAODV2_22Nov2015_MoreSamples_SVFix"
# datasetsFile = "datasets.txt"

#PROCESS = ["MINIV2SVFITESSENTIAL"]
#tag = "MC_SVFit_MiniAODV2_22Nov2015_EssentialSamples_SVFix"
#datasetsFile = "datasets.txt"

# PROCESS = ["SILVERJSONMC"]
# tag = "MC_SilverJson_SVfit"
# datasetsFile = "datasets.txt"

# PROCESS = ["SILVERJSONDATA"]
# tag = "Data_SilverJson_SVfit"
# datasetsFile = "datasets.txt"

# PROCESS = ["DYNLO"]
# tag = "DY_NLO_noSVfit"
# datasetsFile = "datasets.txt"

# PROCESS = ["MC76X"]
# tag = "MC_76X_15Feb2016"
# datasetsFile = "datasets.txt"

# PROCESS = ["DATA76XZZ"]
# tag = "Data_76X_ZZTolljj_29Mar2016"
# datasetsFile = "datasets.txt"

# PROCESS = ["MC76XZZLUCA"]
# tag = "MC_76X_ZZTolljj_29Mar2016"
# datasetsFile = "datasets.txt"

# PROCESS = ["DATA2016"]
# tag = "data_2016_21Giu"
# datasetsFile = "datasets.txt"

# PROCESS = ["DATA2016RESUB"]
# tag = "data_2016_21Giu_resubSingleEle"
# datasetsFile = "datasets.txt"

# PROCESS = ["MCDY80XRESUB"]
# tag = "MC_2016_24Giu_resubDY_VeryLong"
# datasetsFile = "datasets.txt"

# PROCESS = ["MC80XMOREV2"]
# tag = "MC_2016_27Giu_WAndRadionSub2"
# datasetsFile = "datasets.txt"

# PROCESS = ["DATA2016"]
# tag = "data_2016_22giuJSON_diff_16giuJSON"
# datasetsFile = "datasets.txt"

# PROCESS = ["DATA2016RESUB"]
# tag = "data_2016_22giuJSON_diff_16giuJSON_Tauresub"
# datasetsFile = "datasets.txt"

# PROCESS = ["DATA2016"]
# tag = "data_2016_15lug_NoL1TJSON_diff_8lugJSON_runBrunCrunD"
# datasetsFile = "datasets.txt"

# PROCESS = ["DATA2016"]
# tag = "data_2016_20lug_NoL1TJSON_diff_15lug_NoL1TJSON_runBrunCrunD"
# datasetsFile = "datasets.txt"

# PROCESS = ["DATA2016RESUB"]
# tag = "data_2016_20lug_NoL1TJSON_diff_15lug_NoL1TJSON_runBrunCrunD_resSEle"
# datasetsFile = "datasets.txt"

# PROCESS = ["MC80XSUMMER16LUCA"]
# tag = "MC_Summer16"
# datasetsFile = "datasets.txt"

# PROCESS = ["DATA80XRERECOSETLUCARESUB"]
# tag = "Data_23SepReReco_8Feb2017_res2"
# datasetsFile = "datasets.txt"

# PROCESS = ["DATA80XRERECOSETLUCARESUBMORE"]
# tag = "Data_23SepReReco_8Feb2017_res4"
# datasetsFile = "datasets.txt"

# PROCESS = ["DATA7FEBH"]
# tag = "Data_03FEB2017ReReco_22Feb2017_runH"
# datasetsFile = "datasets.txt"

# PROCESS = ["DATA7FEBBG"]
# tag = "Data_03FEB2017ReReco_22Feb2017_runBG"
# datasetsFile = "datasets.txt"

# PROCESS = ["DATA7FEBBGRESUB"]
# tag = "Data_03FEB2017ReReco_22Feb2017_runBG_res"
# datasetsFile = "datasets.txt"

# PROCESS = ["DATA7FEBHRESUB"]
# tag = "Data_03FEB2017ReReco_22Feb2017_runH_res"
# datasetsFile = "datasets.txt"

# PROCESS = ["SUSYSAMPLE"]
# tag = "MCSUSY_2Apr2017"
# datasetsFile = "datasets.txt"

#PROCESS = ["MC80XGRAVITON"]
#tag = "MC_gravitons_24Apr2017"
#datasetsFile = "datasets.txt"

#PROCESS = ["MC80XRSGRAVITON"]
#tag = "MC_gravitonsRS_29Apr2017"
#datasetsFile = "datasets.txt"

#PROCESS = ["MC_PU12Apr"]
#tag = "MC_PU12Apr"
#datasetsFile = "datasets_Fall17_15May2018.txt"

#PROCESS = ["DATA2017"]
#tag = "SingleMuon2017E_26Jun2018"
#datasetsFile = "datasetsFall17.txt"

#PROCESS = ["SIG_Fall17"]
#tag = "VBFRadion400_29Jun2018"
#datasetsFile = "datasetsFall17.txt"

PROCESS = ["DY_Fall2017"]
tag = "MC_16Oct2018"
datasetsFile = "datasetsFall17.txt"

isMC = True
#twiki page with JSON files info https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmV2015Analysis
#50ns JSON file to be used on 2015B and 2015C PDs - integrated luminosity: 71.52/pb - 18/09/2015
#lumiMaskFileName = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-255031_13TeV_PromptReco_Collisions15_50ns_JSON_v2.txt"

#25ns JSON file to be used on 2015C and 2015D PDs - integrated luminosity: 2.11/fb - 13/11/2015
#lumiMaskFileName = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt"
#lumiMaskFileName = "/home/llr/cms/cadamuro/HiggsTauTauFramework/CMSSW_7_4_7/src/LLRHiggsTauTau/NtupleProducer/test/diffLumiMasks/LumiMask_Diff_2p11fb_minus_1p56_13Nov2015.txt"
#lumiMaskFileName  = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt"

#25ns SILVER JSON : 2.46/fb
# lumiMaskFileName  = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_Silver.txt"

#25ns SILVER JSON : 2.63/fb - dec2016 re-reco
# lumiMaskFileName = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Reprocessing/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_Silver.txt"

# 16 giu Golden JSON 2016
# lumiMaskFileName = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-274443_13TeV_PromptReco_Collisions16_JSON.txt"

# 22 giu Golden JSON 2016 -- 3.99/fb
#lumiMaskFileName = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-275125_13TeV_PromptReco_Collisions16_JSON.txt"

# # 22 giu Golden JSON 2016 MINUS 16 giu Golden JSON 2016
# lumiMaskFileName = "22giuJSON_diff_16giuJSON.txt"

## 8 lug JSON MINUS 22 giu JSON
# compareJSON.py --sub /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-275783_13TeV_PromptReco_Collisions16_JSON.txt /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-275125_13TeV_PromptReco_Collisions16_JSON.txt 8lugJSON_diff_22giuJSON.txt
# lumiMaskFileName = "15lug_NoL1TJSON_diff_8lugJSON.txt"
# lumiMaskFileName = "20lug_NoL1TJSON_diff_15lug_NoL1TJSON.txt"
#lumiMaskFileName = 'Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
#lumiMaskFileName = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-302663_13TeV_PromptReco_Collisions17_JSON.txt'
## 15 May 2018 Golden JSON 2017
# https://twiki.cern.ch/twiki/bin/view/CMS/PdmV2017Analysis
#lumiMaskFileName = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt'
lumiMaskFileName = '/home/llr/cms/amendola/HH2017/CMSSW_9_4_6_patch1/src/LLRHiggsTauTau/NtupleProducer/test/crab3_SingleMuon2017_26Jun2018/crab_SingleMuon_4/results/notFinishedLumis.json'

FastJobs = False # controls number of jobs - true if skipping SVfit, false if computing it (jobs will be smaller)
VeryLong = False # controls time for each job - set to true if jobs contain many real lepton pairs --> request for more grid time
EnrichedToNtuples = False # use only False! Do not create ntuples on CRAB because it is very slow, use tier3
PublishDataset = False # publish dataset; set to false if producing ntuples


###################################################################
#### Automated script starting

# dataset block definition
comment = "#"
sectionBeginEnd = "==="

if EnrichedToNtuples: PublishDataset = False


# check if file with dataset exist
if not os.path.isfile(datasetsFile):
    print "File %s not found!!!" % datasetsFile
    sys.exit()

#check if directory exists
crabJobsFolder = "crab3_" + tag
if os.path.isdir(crabJobsFolder):
    print "Folder %s already exists, please change tag name or delete it" % crabJobsFolder
    sys.exit()

# grep all datasets names, skip lines with # as a comment
# block between === * === are "sections" to be processed

currSection = ""
dtsetToLaunch = []

print " =========  Starting submission on CRAB ========"
print " Parameters: "
print " PROCESS: "
for pr in PROCESS: print "   * " , pr
print " tag: " , tag
print " Fast jobs?: " , FastJobs
print " Publish?: "   , PublishDataset

# READ INPUT FILE
with open(datasetsFile) as fIn:
    for line in fIn:
        line = line.strip() # remove newline at the end and leding/trailing whitespaces

        if not line: #skip empty lines
            continue

        if comment in line:
            continue
        
        #print line        
        words = line.split()
        if len(words) >= 3:
            if words[0] == sectionBeginEnd and words[2] == sectionBeginEnd: 
                currSection = words[1]
        else:
            if currSection in PROCESS:                
                dtsetToLaunch.append(line)

# CREATE CRAB JOBS
os.system ("voms-proxy-init -voms cms")

for name in PROCESS: crabJobsFolder + "_" + name
print crabJobsFolder
os.system ("mkdir %s" % crabJobsFolder)

counter = 1 # appended to the request name to avoid overlaps between datasets with same name e.g. /DoubleEG/Run2015B-17Jul2015-v1/MINIAOD vs /DoubleEG/Run2015B-PromptReco-v1/MINIAOD
outlog = open ((crabJobsFolder + "/submissionLog.txt"), "w")
outlog.write (" =========  Starting submission on CRAB ========\n")
outlog.write (" Parameters: \n")
outlog.write (" PROCESS: \n")
for pr in PROCESS: outlog.write ("   * %s\n" % pr)
outlog.write (" tag: %s\n" % tag)
outlog.write (" Fast jobs?: %s\n" % str(FastJobs))
outlog.write (" Publish?: %s\n"   % str(PublishDataset))
outlog.write (" ===============================================\n\n\n")

print dtsetToLaunch
for dtset in dtsetToLaunch:
    dtsetNames = dtset
    if '/MINIAODSIM' in dtset:
        dtsetNames = dtset.replace('/MINIAODSIM', "")
    elif '/MINIAOD' in dtset:
        dtsetNames = dtset.replace('/MINIAOD', "")
    dtsetNames = dtsetNames.replace('/', "__")
    dtsetNames = dtsetNames.strip("__") # remove leading and trailing double __ 
    shortName = dtset.split('/')[1]

    if (len(shortName) > 95): # requestName not exceed 100 Characters!
        toRemove = len (shortName) - 95
        shortName = shortName[toRemove:]

    #dtSetName = dtsetNames[1]
    command = "crab submit -c crab3_template.py"
    #command += " General.requestName=%s" % (dtsetNames + "_" + tag + "_" + str(counter))
    command += " General.requestName=%s" % (shortName + "_" + str(counter))
    command += " General.workArea=%s" % crabJobsFolder
    command += " Data.inputDataset=%s" % dtset
    #command += " Data.outLFNDirBase=/store/user/lcadamur/HHNtuples/%s/%s" % (tag , str(counter)+"_"+dtsetNames)
    command += " Data.outLFNDirBase=/store/user/camendol/HHNtuples2017/%s/%s" % (tag , str(counter)+"_"+dtsetNames)
    #command += " Data.outLFNDirBase=/store/user/fbrivio/Hhh_1718/%s/%s" % (tag , str(counter)+"_"+dtsetNames) # change to where you want to stage you ntuples
    command += " Data.outputDatasetTag=%s" % (shortName + "_" + tag + "_" + str(counter))
    #command += " Data.splitting='Automatic'"
    if (EnrichedToNtuples): command += " Data.inputDBS=phys03" # if I published the dataset need to switch from global (default)
    if (EnrichedToNtuples): command += " JobType.psetName=ntuplizer.py" # run a different python config for enriched
    if not PublishDataset : command += " Data.publication=False" # cannot publish flat root ntuples
    if (FastJobs)         : command += " Data.unitsPerJob=100000" # circa 50 ev / secondo --> circa 1/2 h ; else leave default of 4000 jobs
    if VeryLong           : command += " JobType.maxJobRuntimeMin=2500" # 32 hours, default is 22 hours -- can do up to 2800 hrs
    if not isMC           : command += " Data.lumiMask=%s" % lumiMaskFileName
    print command ,  "\n"
    os.system (command)
    outlog.write(command + "\n\n")
    counter = counter + 1
