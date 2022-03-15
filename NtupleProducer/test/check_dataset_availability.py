#!/usr/bin/env python3
# coding: utf-8

"""
Checks if all files of one or multiple datasets are available on non-blacklisted sites.
"""

import os
import sys
import re
import json
import subprocess
import urllib
import urllib.request
import json

import tqdm


# dynamic file listing all usable and blacklisted sites
usable_sites_url = "https://cmssst.web.cern.ch/cmssst/analysis/usableSites.json"

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

# cache for the live black list
_live_black_list = None


def check_dataset(dataset, white_list=None, live_black_list=False, via_json=False, das_host=None):
    print(f"checking dataset {dataset}")

    # get dataset lfns
    print("retrieving LFNs ...")
    lfns = get_dataset_lfns(dataset, via_json=via_json, das_host=das_host)
    print(f"done, found {len(lfns)} LFNs")

    # get storage sites
    print("checking sites per LFN ...")
    lfn_sites = {
        lfn: get_lfn_sites(lfn, via_json=via_json, das_host=das_host)
        for lfn in tqdm.tqdm(lfns)
    }
    print("done")

    # per lfn, check if there is at least one allowed site
    unreachable_lfns = [
        lfn for lfn, sites in lfn_sites.items()
        if not any(
            site_is_allowed(site, white_list=white_list, live_black_list=live_black_list)
            for site in sites
        )
    ]

    # log results
    if unreachable_lfns:
        print(f"found {len(unreachable_lfns)} LFNs that are not reacheable through allowed sites:")
        for i, lfn in enumerate(unreachable_lfns):
            print(f"{i + 1}.")
            print(f"    LFN  : {lfn}")
            print(f"    Sites: {str(lfn_sites[lfn])[1:-1]}")
    else:
        print("all LFNs are reachable through allowed sites")


def site_is_allowed(site, white_list=None, live_black_list=False, allow_tape=False):
    # check white list
    if white_list and site in white_list:
        return True

    # check tape
    if not allow_tape and re.match(r"^t\d_[^_]+_[^_]+_.*tape.*$", site.lower()):
        return False

    # last check through black list
    black_list = load_site_black_list() if live_black_list else global_black_list
    return site not in black_list


def get_dataset_lfns(dataset, via_json=False, das_host=None):
    data = run_das_query(f"file dataset={dataset}", via_json=via_json, das_host=das_host)

    if via_json:
        lfns = [str(entry["file"][0]["name"]) for entry in data]
    else:
        lfns = data

    return lfns


def get_lfn_sites(lfn, via_json=False, das_host=None):
    data = run_das_query(f"site file={lfn}", via_json=via_json, das_host=das_host)

    if via_json:
        sites = [str(entry["site"][0]["name"]) for entry in data]
    else:
        sites = data

    return sites


def run_das_query(query, via_json=False, das_host=None):
    # build and run the command
    cmd = f"dasgoclient -query=\"{query}\""
    if das_host:
        cmd += f" -host=\"{das_host}\""
    if via_json:
        cmd += " -json"

    p, out = run_command(cmd, attempts=2, timeout=30)

    # expect the return code to be 0
    if p.returncode != 0:
        raise Exception(f"command failed: {cmd}")

    # catch certain DAS based errors
    if out.startswith("jsonparser failure"):
        raise Exception(f"DAS responded with error: {out}")

    if via_json:
        data = json.loads(out)
    else:
        data = [line.strip() for line in out.strip().split() if line.strip()]

    return data


def run_command(cmd, timeout=None, attempts=1, _attempt=1, **kwargs):
    kwargs["preexec_fn"] = os.setsid
    p = subprocess.Popen(cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, **kwargs)

    try:
        out, _ = p.communicate(timeout=timeout)
        return p, out.decode("utf-8")

    except KeyboardInterrupt:
        p.terminate()
        try:
            pgid = os.getpgid(p.pid)
            if p.poll() is None:
                os.killpg(pgid, signal.SIGTERM)
        except:
            pass
        raise

    except Exception:
        if attempts > _attempt:
            return run_command(cmd, timeout=timeout, attempts=attempts, _attempt=_attempt + 1)

        print(f"command failed after {attempts} attempt(s):\n{cmd}")
        raise


def load_sites():
    content = urllib.request.urlopen(usable_sites_url).read()
    sites = json.loads(content)

    # return a list of dicts {"name": str, "url": str, "black_listed": bool}
    # define sites as black-listed when their "value" is not "usable"
    return [
        {
            "name": str(site["name"]),
            "url": str(site["url"]),
            "black_listed": site.get("value", "not_usable") != "usable",
        }
        for site in sites
    ]


def load_site_black_list():
    global _live_black_list

    if _live_black_list is None:
        _live_black_list = [site["name"] for site in load_sites() if site["black_listed"]]

    return _live_black_list


def main():
    import argparse

    csv = lambda s: [_s.strip() for _s in s.strip().split(",")]

    # setup arguments
    parser = argparse.ArgumentParser(description="dataset availability check")
    parser.add_argument("dataset", nargs="*", metavar="DATASET", help="the dataset(s) to check")
    parser.add_argument("--white-list", "-w", type=csv, metavar="WHITE_LIST", help="an optional "
        "white list that is evaluated on top of the global black list")
    parser.add_argument("--live-black-list", "-l", action="store_true", help="use the live site "
        "black list from the site status board")
    parser.add_argument("--allow-tape", "-t", action="store_true", help="allow tape storage")
    parser.add_argument("--via-json", "-j", action="store_true", help="run queries via json format")
    parser.add_argument("--das-host", help="a custom DAS hostname to use, see 'dasgoclient --help'")
    parser.add_argument("--print-live-black-list", action="store_true", help="print the live site "
        "black list and exit")
    args = parser.parse_args()

    # maybe just print the live black list
    if args.print_live_black_list:
        print("\n".join(load_site_black_list()))
        sys.exit(0)

    # at least one dataset is required
    if not args.dataset:
        print("at least one DATASET is required")
        sys.exit(1)

    # some insightful prints
    print(f"number of datasets   : {len(args.dataset)}")
    if args.white_list:
        print(f"custom white list    : {str(args.white_list)[1:-1]}")
    print(f"use live black list  : {args.live_black_list}")
    print(f"tape storage allowed : {args.allow_tape}")
    print(f"use json queries     : {args.via_json}")
    if args.das_host:
        print(f"custom DAS host      : {args.das_host}")

    # loop through datasets and peform the check
    for i, dataset in enumerate(args.dataset):
        print("")
        check_dataset(dataset, white_list=args.white_list, live_black_list=args.live_black_list,
            via_json=args.via_json, das_host=args.das_host or None)
        if i != len(args.dataset) - 1:
            print(100 * "-")


if __name__ == "__main__":
    main()
