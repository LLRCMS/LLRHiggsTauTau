# coding: utf-8

import os
import sys
import re
import subprocess
from collections import OrderedDict

from util import TeeStream, hash_truncate, get_wlcg_user


#
# dataset selection
#

datasets_file = "datasets_UL17.txt"
year = 2017
tag = "uhh_2017_v2"
processes = [
    # "DATA_TAU_2017",
    # "DATA_ELE_2017",
    # "DATA_MU_2017",
    # "DATA_MET_2017",
    # "BACKGROUNDS_TT_2017",
    # "BACKGROUNDS_WJETS_2017",
    # "BACKGROUNDS_DY_2017",
    # "BACKGROUNDS_VV_2017",
    # "BACKGROUNDS_VVV_2017",
    # "BACKGROUNDS_ST_2017",
    # "BACKGROUNDS_EWK_2017",
    # "BACKGROUNDS_H_2017",
    # "BACKGROUNDS_TTX_2017",
    # "SIGNALS_GF_NONRES_2017",
    # "SIGNALS_VBF_NONRES_2017",
    # "SIGNALS_GF_NLO_NONRES_2017",
    # "SIGNALS_GF_SPIN0_2017",
    # "SIGNALS_GF_SPIN2_2017",
    # "SIGNALS_VBF_SPIN0_2017",
    # "SIGNALS_VBF_SPIN2_2017",
]


#
# configuration
#

# LFN store path
lfn_base = "/store/user/{}/hbt_resonant_run2/HHNtuples".format(get_wlcg_user())
# dry run, runs all preparations but skips the submission
dry_run = False
# ignore dataset locality rules between CEs and SEs, and define the required whitelist
ignore_locality = False
site_white_list = [
    "T1_DE_KIT", "T2_DE_DESY", "T2_DE_RWTH",
    "T2_CH_CERN",
    "T1_US_FNAL", "T2_US_Caltech", "T2_US_Vanderbilt", "T2_US_Wisconsin",
]
# controls number of jobs - true if skipping SVfit, false if computing it (jobs will be smaller)
fast_jobs = False
# controls time for each job - set to true if jobs contain many real lepton pairs --> request for more grid time
very_long = False
# use only False! Do not create ntuples on CRAB because it is very slow, use tier3
enriched_to_ntuples = False
# set to false if producing ntuples
publish_dataset = False


#
# submission logic
#

# check if file with dataset exist
if not os.path.isfile(datasets_file):
    print("file %s not found" % datasets_file)
    sys.exit(2)

# check if directory exists
crab_jobs_folder = "jobs_crab3_" + tag
if not os.path.isdir(crab_jobs_folder):
    os.makedirs(crab_jobs_folder)

# read datasets for selected processes
curr_section = ""
datasets = []
with open(datasets_file) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        # section or dataset line
        m = re.match(r"^===\s+(.+)\s+===.*$", line)
        if m:
            curr_section = m.group(1)
        elif curr_section and curr_section in processes:
            mc = not curr_section.startswith("DATA_")
            datasets.append((line, mc))

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

    for dataset, mc in datasets:
        # build request names and dataset tags to construct the output directory according to
        # <lfn-prefix>/<primary-dataset>/<publication-name>/<time-stamp>/<counter>[/log]/<file-name>
        sample_name, campaign_name = dataset.split("/")[1:3]
        # build the request name from tag, sample and campaign name
        request_name = hash_truncate("{}__{}__{}".format(tag, sample_name, campaign_name), 100)
        # build the dataset tag that, when defined, takes the role of the publication name
        dataset_tag = campaign_name
        if publish_dataset:
            dataset_tag = "{}__{}__{}".format(tag, sample_name, campaign_name)

        # complain when the job dir is existing
        job_dir = os.path.join(crab_jobs_folder, "crab_" + request_name)
        if os.path.exists(job_dir):
            print("job directory for tag '{}', dataset '{}' already existing at {}, skip\n".format(
                tag, dataset, job_dir))
            continue

        # create crab arguments
        crab_args = OrderedDict()

        crab_args["General.requestName"] = request_name
        crab_args["General.workArea"] = crab_jobs_folder
        crab_args["Data.inputDataset"] = dataset
        crab_args["Data.outLFNDirBase"] = "{}/{}".format(lfn_base, tag)
        crab_args["Data.outputDatasetTag"] = dataset_tag
        crab_args["Data.publication"] = publish_dataset
        crab_args["Data.ignoreLocality"] = ignore_locality
        crab_args["JobType.pyCfgParams"] = ["year={}".format(year), "mc={}".format(mc)]

        if ignore_locality:
            crab_args["Site.whitelist"] = site_white_list

        # for testing
        # crab_args["Data.splitting"] = "FileBased"
        # crab_args["Data.unitsPerJob"] = 1
        # crab_args["Data.totalUnits"] = 1

        # if enriched_to_ntuples:
        #     crab_args["Data.inputDBS"] = "phys03"
        #     crab_args["JobType.psetName"] = "ntuplizer.py"
        # if fast_jobs:
        #     # circa 50 ev / secondo --> circa 1/2 h ; else leave default of 4000 jobs
        #     crab_args["Data.unitsPerJob"] = 100000
        # if very_long:
        #     # 32 hours, default is 22 hours -- can do up to 2800 hrs
        #     crab_args["JobType.maxJobRuntimeMin"] = 2500

        # build the full command
        fmt_list = lambda a: "[{}]".format(",".join(map("'{}'".format, a)))
        cms_run_arg = lambda k, a: "{}=\"{}\"".format(k, fmt_list(a) if isinstance(a, list) else a)
        command = "crab submit -c crab3_template_uhh.py"
        command += " " + " ".join(cms_run_arg(*tpl) for tpl in crab_args.items())

        log.write("submission command:\n%s" % command)
        if dry_run:
            log.write("dry run, skip command execution")
        else:
            subprocess.call(command, shell=True, stderr=subprocess.STDOUT, stdout=log)
        log.write("\n")
