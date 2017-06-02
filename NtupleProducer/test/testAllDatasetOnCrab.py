# use this script to test the status of the submitted tasks on CRAB

# NB: you need to source crab environment before lauching this script:
# source /cvmfs/cms.cern.ch/crab3/crab.sh

import os, glob
import time
import datetime as dt
from subprocess import Popen, PIPE

class ProcDataset:
    
    def __init__ (self):
        self.FolderName = ""
        self.DASurl = ""
        self.DatasetName = ""
        self.Nfailed = ""
        self.Nrunning = ""
        self.Nfinished = ""
        self.TaskStatus = ""

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
    

#path = "crab3_produzione_MC_3Ago2015"
#path = "crab3_produzione_DATA_3Ago2015"
#path = "crab3_llrNtuples_partial_5Ago"

#path = "crab3_llrNt_NoSVFit_bkg_27Ago2015"
#path = "crab3_llrNt_NoSVFit_data_27Ago2015"
#path = "crab3_llrNt_NoSVFit_data_30Ago2015_lumiMaskFix"
#path = "crab3_Data50ns_SVFit_6Ott2015"
#path = "crab3_Data25ns_SVFit_6Ott2015"
#path = "crab3_Data25ns_noSVFit_27Ott2015"
#path = "crab3_Data25ns_noSVFit_npvFix_30Ott2015"
#path = "crab3_Data25ns_noSVFit_npvFix_30Ott2015_newJson"
#path = "crab3_Data25ns_noSVFit_npvFix_30Ott2015_newJson_trigFix"
#path = "crab3_Data25ns_noSVFit_npvFix_30Ott2015_newJson_trigFix2"
#path = "crab3_MC_NoSVFit_MiniAODV2_13Nov2015"

#path = "crab3_MC_SVFit_MiniAODV2_22Nov2015"
#path = "crab3_Data_SVFit_MiniAODV2_22Nov2015"
#path = "crab3_MC_2016_24Giu_resubDY"
#path = "crab3_MC_2016_24Giu_resubDY_VeryLong"
#path = "crab3_MC_2016_27Giu_WAndRadionSub2"
#path = "crab3_data_2016_08lugJSON_diff_22giuJSON_runBrunC"
#path = "crab3_data_2016_15lug_NoL1TJSON_diff_8lugJSON_runBrunCrunD"

#path = "crab3_MC_Summer16"
#path = "crab3_Data_23SepReReco_8Feb2017_res"
path = "crab3_MC_gravitons_24Apr2017"

if not os.path.isdir(path):
    print "Folder %s does not exist" % path
    sys.exit()

filesLevel1 = glob.glob('%s/*' % path)
dirsLevel1 = filter(lambda f: os.path.isdir(f), filesLevel1) # this is the list of folders under "path"

#os.system ("voms-proxy-init -voms cms")

datasets = [] #list of all processed datasets, filled from the list of folders

localtime = time.asctime( time.localtime(time.time()) )
fulllogname = (path + "/fullStatusLog_{}.txt").format( dt.datetime.now().strftime('%Y.%m.%d_%H.%M.%S') )
fulllog = open (fulllogname, "w")
fulllog.write ("Local current time : %s\n" % localtime)

for dirr in dirsLevel1:
    print dirr
    fulllog.write ("=== DIRECTORY TASK IS: %s\n" % dirr)
    procdataset = ProcDataset()
    procdataset.addFolderName (dirr)
    command = "crab status -d %s" % dirr
    pipe = Popen(command, shell=True, stdout=PIPE)
    for line in pipe.stdout:
        fulllog.write(line)
        line = line.strip()
        print line
        blocks = line.split(':', 1) # split at 1st occurrence
        if blocks [0] == "Output dataset":
            dtsetName = blocks[1].strip()
            procdataset.addDatasetName (dtsetName)
        elif blocks [0] == "Task status":
            status = blocks[1].strip()
            procdataset.TaskStatus = status
        elif blocks [0] == "Output dataset DAS URL":
            procdataset.DASurl = blocks[1].strip()
    datasets.append(procdataset)
    print "\n\n"
    fulllog.write("\n\n")

print "\n =================  TASK STATUS  ==================== \n\n"
print "failed: "
for dt in datasets:
    status = dt.TaskStatus
    if status == "FAILED":
        print "crab resubmit -d %s" % dt.FolderName

print "\n =================  PUBLICATION  ==================== \n\n"

noDataset = []
for dt in datasets:
    name = dt.DatasetName
    if name:
        print name
    else:
        noDataset.append(dt)
        print "\n"

print "\nNo Dataset name found for:"
for dt in noDataset:
    print dt.FolderName


print "\n =================  DAS URLS  ==================== \n\n"

#noUrl = []
for dt in datasets:
    url = dt.DASurl
    if url:
        print url
    else:
        #noUrl.append(dt)
        print "\n"

#print "\nNo Dataset name found for:"
#for dt in noDataset:
#    print dt.FolderName

