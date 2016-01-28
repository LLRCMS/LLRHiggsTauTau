# makes all the file lists fot eh datasets (only files that are finished, could be partial)

import os
import sys
import re
import subprocess

PROCESS = ["PROD_PARTIAL"]
tag = "llrNtuples_partial_5Ago_resub"
datasetsFile = "../datasets_Enriched.txt"


###################################################################
#### Automated script starting

# dataset block definition
comment = "#"
sectionBeginEnd = "==="

# check if file with dataset exist
if not os.path.isfile(datasetsFile):
    print "File %s not found!!!" % datasetsFile
    sys.exit()

# grep all datasets names, skip lines with # as a comment
# block between === * === are "sections" to be processed

currSection = ""
dtsetToLaunch = []

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

# for each dataset call the python script
#usage: source MakeFileListDAS.sh -t "Data_MuMu" -o Data_MuMu.py -p /DoubleMuon/Run2015B-PromptReco-v1/MINIAOD
outDir = "../inputFilesEMiniAOD/" + tag
if os.path.isdir (outDir):
    print "Folder %s already esists, please change tag or delete it" % outDir
    sys.exit()
os.system('mkdir -p %s' % outDir)
makeExe = "chmod u+x MakeFileListDAS.sh"
os.system(makeExe)
for dtset in dtsetToLaunch:
    dtsetModif =  (dtset.split('/')[1]).strip() # might have conflicts if processing two datasets with the same name
    thistag = tag + " ON dtset: " + dtsetModif 
    thisoutput = outDir + "/" + dtsetModif + ".py"
    command = './MakeFileListDAS.sh -t \"%s\" -o %s -p %s -d' % (thistag, thisoutput, dtset)
    print command
    os.system(command)
