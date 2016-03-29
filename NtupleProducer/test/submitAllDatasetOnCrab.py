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

PROCESS = ["DATA76X"]
tag = "Data_76X_15Feb2016"
datasetsFile = "datasets.txt"

isMC = False
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
lumiMaskFileName = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Reprocessing/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_Silver.txt"


FastJobs = True # controls number of jobs - true if skipping SVfit, false if computing it (jobs will be smaller)
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
    command += " Data.outLFNDirBase=/store/user/lcadamur/HHNtuples/%s/%s" % (tag , str(counter)+"_"+dtsetNames)
    command += " Data.outputDatasetTag=%s" % (dtsetNames + "_" + tag + "_" + str(counter))
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
