# use this script to test the status of the submitted tasks on CRAB

# NB: you need to source crab environment before lauching this script:
# source /cvmfs/cms.cern.ch/crab3/crab.sh

import os, glob
from subprocess import Popen, PIPE

class ProcDataset:
    
    def __init__ (self):
        self.FolderName = ""
        self.DASurl = ""
        self.DatasetName = ""
        self.Nfailed = ""
        self.Nrunning = ""
        self.Nfinished = ""

    def addFolderName (self, name):
        self.FolderName = name

    def addDASurl (self, url):
        self.DASurl = url
    
    def addDatasetName (self, dname):
        self.DatasetName = dname
    
    def setFailed (self, num):
        self.Nfailed = num

    def setRunning (self, num):
        self.Nrunning = num

    def setFinished (self, num):
        self.Nfinished = num

path = "crab3_produzione_MC_3Ago2015"
#path = "crab3_llrNtuples_partial_5Ago"

if not os.path.isdir(path):
    print "Folder %s does not exist" % path
    sys.exit()

filesLevel1 = glob.glob('%s/*' % path)
dirsLevel1 = filter(lambda f: os.path.isdir(f), filesLevel1) # this is the list of folders under "path"

#os.system ("voms-proxy-init -voms cms")

datasets = [] #list of all processed datasets, filled from the list of folders

for dirr in dirsLevel1:
    print dirr
    procdataset = ProcDataset()
    procdataset.addFolderName (dirr)
    command = "crab status -d %s" % dirr
    pipe = Popen(command, shell=True, stdout=PIPE)
    for line in pipe.stdout:
        line = line.strip()
        blocks = line.split(':')
        if blocks [0] == "Output dataset":
            dtsetName = blocks[1].strip()
            procdataset.addDatasetName (dtsetName)
    datasets.append(procdataset)

print " ================================== "
noDataset = []
for dt in datasets:
    name = dt.DatasetName
    if name:
        print name
    else:
        noDataset.append(dt)

print "\nNo Dataset name found for:"
for dt in noDataset:
    print dt.FolderName
