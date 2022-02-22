# coding: utf-8

import os
import sys
import re


###################################################################
#### Customization

PROCESS = ["BACKGROUNDS_DY_2017"]
tag = "test_dy"
datasetsFile = "datasets_UL17.txt"
isMC = False

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
    shortName = dtset.split('/')[1][:95]

    #dtSetName = dtsetNames[1]
    command = "crab submit -c crab3_template_uhh.py"
    #command += " General.requestName=%s" % (dtsetNames + "_" + tag + "_" + str(counter))
    command += " General.requestName=%s" % (shortName + "_" + str(counter))
    command += " General.workArea=%s" % crabJobsFolder
    command += " Data.inputDataset=%s" % dtset
    #command += " Data.outLFNDirBase=/store/user/lcadamur/HHNtuples/%s/%s" % (tag , str(counter)+"_"+dtsetNames)
    #command += " Data.outLFNDirBase=/store/user/camendol/HHNtuples2017/%s/%s" % (tag , str(counter)+"_"+dtsetNames)
    #command += " Data.outLFNDirBase=/store/user/fbrivio/Hhh_1718/%s/%s" % (tag , str(counter)+"_"+dtsetNames) # change to where you want to stage you ntuples
    command += " Data.outLFNDirBase=/store/user/mrieger/hbt_resonant_run2/HHNtuples/%s/%s" % (tag , str(counter)+"_"+dtsetNames) # change to where you want to stage you ntuples
    command += " Data.outputDatasetTag=%s" % (shortName + "_" + tag + "_" + str(counter))
    #command += " Data.splitting='Automatic'"
    if (EnrichedToNtuples): command += " Data.inputDBS=phys03" # if I published the dataset need to switch from global (default)
    if (EnrichedToNtuples): command += " JobType.psetName=ntuplizer.py" # run a different python config for enriched
    if not PublishDataset : command += " Data.publication=False" # cannot publish flat root ntuples
    # if (FastJobs)         : command += " Data.unitsPerJob=100000" # circa 50 ev / secondo --> circa 1/2 h ; else leave default of 4000 jobs
    # if VeryLong           : command += " JobType.maxJobRuntimeMin=2500" # 32 hours, default is 22 hours -- can do up to 2800 hrs
    print command ,  "\n"
    os.system (command)
    outlog.write(command + "\n\n")
    counter = counter + 1
