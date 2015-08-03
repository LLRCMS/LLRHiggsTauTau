import os
import sys
import re

# NB: you need to source crab environment before lauching this script:
# source /cvmfs/cms.cern.ch/crab3/crab.sh

###################################################################
#### Parameters to be changed for 

#PROCESS = ["BKGPU40", "OPT"]
PROCESS = ["TEST"]
datasetsFile = "datasets.txt"
tag = "prod_di_test" # a folder with cra3_[tag] is created, ad appended in many places


###################################################################
#### Automated script starting

# dataset block definition
comment = "#"
sectionBeginEnd = "==="

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

# READ INPUT FILE
with open(datasetsFile) as fIn:
    for line in fIn:
        line = line [:-1] # remove newline at the end
        
        if not line: #skip empty lines
            continue

        if comment in line:
            continue
        
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

for dtset in dtsetToLaunch:
    dtsetNames = dtset.replace('/MINIAODSIM', "")
    dtsetNames = dtset.replace('/', "__")
    #dtSetName = dtsetNames[1]
    command = "crab submit -c crab3_template.py"
    command += " General.requestName=%s" % (dtsetNames + "_" + tag)
    command += " General.workArea=%s" % crabJobsFolder
    command += " Data.inputDataset=%s" % dtset
    command += " Data.outLFNDirBase=/store/user/lcadamur/HHNtuples/%s" % tag
    command += " Data.publishDataName=%s" % tag
    #print command
    os.system (command)