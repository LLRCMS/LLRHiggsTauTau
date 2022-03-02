# coding: utf-8

import os
import sys
import re
import getpass

from util import TeeStream, get_voms_proxy_user


###################################################################
#### Customization

datasets_file = "datasets_UL17.txt"
processes = ["BACKGROUNDS_DY_2017"]
tag = "test_dy"
lfn_base = "/store/user/{}/hbt_resonant_run2/HHNtuples".format(get_voms_proxy_user())

# dry run, runs all preparations but skips the submission
dry_run = False
# controls number of jobs - true if skipping SVfit, false if computing it (jobs will be smaller)
fast_jobs = False
# controls time for each job - set to true if jobs contain many real lepton pairs --> request for more grid time
very_long = False
# use only False! Do not create ntuples on CRAB because it is very slow, use tier3
enriched_to_ntuples = False
# set to false if producing ntuples
publish_dataset = False


###################################################################
#### Automated script starting

# dataset block definition
comment = "#"
section_begin_end = "==="

if enriched_to_ntuples:
    publish_dataset = False

# check if file with dataset exist
if not os.path.isfile(datasets_file):
    print("file %s not found" % datasets_file)
    sys.exit(1)

#check if directory exists
crab_jobs_folder = "crab3_" + tag
if os.path.isdir(crab_jobs_folder):
    print("folder %s already exists, please change tag name or delete it" % crab_jobs_folder)
    sys.exit(2)
os.makedirs(crab_jobs_folder)

# read datasets
curr_section = ""
datasets = []
with open(datasets_file) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if comment in line:
            continue

        words = line.split()
        if len(words) >= 3:
            if words[0] == section_begin_end and words[2] == section_begin_end:
                curr_section = words[1]
        else:
            if curr_section in processes:
                datasets.append(line)

# start submitting them and log simultaneously
with TeeStream(os.path.join(crab_jobs_folder, "submissionLog.txt")) as log:
    log.write(" =========  Starting submission on CRAB ========")
    log.write(" tag      : %s" % tag)
    log.write(" processes: ")
    for pr in processes:
        log.write("   * %s" % pr)
    log.write(" crab dir : %s" % crab_jobs_folder)
    log.write(" fast jobs: %s" % fast_jobs)
    log.write(" publish  : %s" % publish_dataset)
    log.write("")

    for counter, dataset in enumerate(datasets[:1]):
        sample_name, campaign_name = dataset.split("/")[1:3]
        request_name = "{}_{}".format(sample_name[:95], counter)
        dataset_tag = "{}_{}".format(tag, request_name) if publish_dataset else campaign_name

        # start building the command
        command = "crab submit -c crab3_template_uhh.py"

        command += " General.requestName=%s" % request_name
        command += " General.workArea=%s" % crab_jobs_folder

        command += " Data.inputDataset=%s" % dataset
        command += " Data.outLFNDirBase=%s/%s" % (lfn_base, tag)
        command += " Data.outputDatasetTag=%s" % dataset_tag
        command += " Data.publication=%s" % publish_dataset

        # if enriched_to_ntuples:
        #     command += " Data.inputDBS=phys03"
        #     command += " JobType.psetName=ntuplizer.py"
        # if fast_jobs:
        #     # circa 50 ev / secondo --> circa 1/2 h ; else leave default of 4000 jobs
        #     command += " Data.unitsPerJob=100000"
        # if very_long:
        #     # 32 hours, default is 22 hours -- can do up to 2800 hrs
        #     command += " JobType.maxJobRuntimeMin=2500"

        log.write("submission command:\n%s" % command)
        if dry_run:
            log.write("dry run, skip command execution")
        else:
            os.system(command)
        log.write("\n")
