#!/usr/bin/env python
# coding: utf-8

"""
Checks if all files of one or multiple datasets are available on non-blacklisted sites.
"""

import re
import json
import subprocess

import tqdm


# the global black list, please update if necessary
global_black_list = {
    "T3_US_UCR", "T3_FR_IPNL", "T3_CH_PSI", "T3_IN_PUHEP", "T3_RU_FIAN", "T3_CH_CERN_HelixNebula",
    "T3_CN_PKU", "T2_GR_Ioannina", "T3_US_UMiss", "T3_US_ANL", "T3_US_OSG", "T2_US_MIT",
    "T3_US_SDSC", "T2_AT_Vienna", "T2_CH_CERN_HLT", "T3_IN_TIFRCloud", "T3_US_Brown", "T3_US_UCD",
    "T2_RU_ITEP", "T3_US_NERSC", "T2_IT_Pisa", "T2_RU_SINP", "T3_US_Rutgers", "T0_CH_CSCS_HPC",
    "T3_CH_CERNBOX", "T3_US_FNALLPC", "T3_CH_CERN_OpenData", "T0_CH_CERN", "T3_KR_KISTI",
    "T3_US_Anvil", "T3_US_PSC", "T2_BE_IIHE", "T3_IR_IPM", "T3_US_Kansas", "T3_BY_NCPHEP",
    "T3_BG_UNI_SOFIA", "T3_US_CMU", "T2_US_Nebraska", "T3_TW_NTU_HEP", "T3_IT_MIB", "T3_US_TACC",
    "T3_UK_SGrid_Oxford", "T3_US_JHU", "T3_US_MIT", "T3_US_NotreDame", "T3_US_Princeton_ICSE",
    "T3_HR_IRB", "T3_KR_UOS", "T2_UA_KIPT", "T3_GR_IASA",
}


def check_dataset(dataset, white_list=None, via_json=False):
    print("checking dataset {}".format(dataset))

    # get dataset lfns
    print("retrieving LFNs ...")
    lfns = get_dataset_lfns(dataset, via_json=via_json)
    print("done, found {} LFNs".format(len(lfns)))

    # get storage sites
    print("checking sites per LFN ...")
    lfn_sites = {lfn: get_lfn_sites(lfn, via_json=via_json) for lfn in tqdm.tqdm(lfns)}
    print("done")

    # per lfn, check if there is at least one allowed site
    unreachable_lfns = [
        lfn for lfn, sites in lfn_sites.items()
        if not any(site_is_allowed(site, white_list=white_list) for site in sites)
    ]

    # log results
    if unreachable_lfns:
        print("found {} LFNs that are not reacheable through allowed sites:".format(
            len(unreachable_lfns)))
        for i, lfn in enumerate(unreachable_lfns):
            print("{}.".format(i + 1))
            print("    LFN  : {}".format(lfn))
            print("    Sites: {}".format(", ".join(lfn_sites[lfn])))
    else:
        print("all LFNs are reachable throug allowed sites")


def site_is_allowed(site, white_list=None, allow_tape=False):
    # check white list
    if white_list and site in white_list:
        return True

    # check tape
    if not allow_tape and re.match(r"^t\d_[^_]+_[^_]+_.*tape.*$", site.lower()):
        return False

    # last check through black list
    return site not in global_black_list


def get_dataset_lfns(dataset, via_json=False):
    data = run_das_query("file dataset={}".format(dataset), via_json=via_json)

    if via_json:
        lfns = [str(entry["file"][0]["name"]) for entry in data]
    else:
        lfns = data

    return lfns


def get_lfn_sites(lfn, via_json=False):
    data = run_das_query("site file={}".format(lfn), via_json=via_json)

    if via_json:
        sites = [str(entry["site"][0]["name"]) for entry in data]
    else:
        sites = data

    return sites


def run_das_query(query, via_json=False):
    # build and run the command
    cmd = "dasgoclient -query=\"{}\"".format(query)
    if via_json:
        cmd += " -json"

    p = subprocess.Popen(cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE)

    try:
        out, _ = p.communicate()
    except KeyboardInterrupt:
        try:
            p.terminate()
        except:
            pass
        raise

    if via_json:
        data = json.loads(out)
    else:
        data = [line.strip() for line in out.strip().split() if line.strip()]

    return data


def main():
    import argparse

    csv = lambda s: [_s.strip() for _s in s.strip().split(",")]

    # setup arguments
    parser = argparse.ArgumentParser(description="dataset availability check")
    parser.add_argument("dataset", nargs="+", metavar="DATASET", help="the dataset(s) to check")
    parser.add_argument("--white-list", "-w", type=csv, metavar="WHITE_LIST", help="an optional "
        "white list that is evaluated on top of the global black list")
    parser.add_argument("--allow-tape", "-t", action="store_true", help="allow tape storage")
    parser.add_argument("--via-json", "-j", action="store_true", help="run queries via json format")
    args = parser.parse_args()

    # some insightful prints
    print("start dataset checking")
    if args.white_list:
        print("custom white list: {}".format(", ".join(args.white_list)))
    print("tape storage is {}allowed".format("" if args.allow_tape else "not "))
    print("using {} queries".format("json" if args.via_json else "plain"))

    # loop through datasets and peform the check
    for dataset in args.dataset:
        print("")
        check_dataset(dataset, white_list=args.white_list, via_json=args.via_json)
        print(100 * "-")


if __name__ == "__main__":
    main()
