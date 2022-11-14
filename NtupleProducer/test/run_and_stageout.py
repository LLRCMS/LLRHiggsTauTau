#!/usr/bin/env python3
# coding: utf-8

"""
Script that runs an analyzer config with a certain input file and stages out the output file
to a remote storage element. Example:

> ./run_and_stageout.py -c -f \
      /store/data/Run2017B/MET/MINIAOD/UL2017_MiniAODv2-v1/260000/DC7B6B63-FFBC-FA4C-A8C0-28C1DC70463D.root \
      gsiftp://dcache-door-cms04.desy.de//pnfs/desy.de/cms/tier2/store/user/mrieger/hbt_resonant_run2/HHNtuples/uhh_2017_v2/MET/Run2017B-UL2017_MiniAODv2-v1/220311_220230/0000/HTauTauAnalysis_131.root \
      year=2017 mc=False
"""

import os
import re
import tempfile
import shutil

from check_dataset_availability import run_command


thisdir = os.path.dirname(os.path.abspath(__file__))


def run_and_stageout(input_file, output_file, analyzer_path, analyzer_options=None,
        change_dir=False, fetch_input=False, dry_run=False):
    # input checks
    input_scheme = get_scheme(input_file)
    input_is_remote = input_file.startswith("/store/") or (input_scheme and input_scheme != "file")
    if not input_scheme and not input_is_remote:
        input_file = "file://" + os.path.abspath(input_file)
    if not get_scheme(output_file):
        output_file = "file://" + os.path.abspath(output_file)
    analyzer_path = os.path.abspath(analyzer_path)
    if not os.path.isfile(analyzer_path):
        raise IOError(f"no such file: {analyzer_path}")
    analyzer_dir = os.path.dirname(analyzer_path)
    tmp_dir = tempfile.mkdtemp()

    # some logs
    print(f"input     : {input_file}")
    print(f"output    : {output_file}")
    print(f"analyzer  : {analyzer_path}")
    print(f"options   : {analyzer_options}")
    print(f"change dir: {change_dir}")
    print(f"tmp dir   : {tmp_dir}")

    # define a temporary output file
    _output_file = os.path.join(tmp_dir, "output.root")

    # optinally fetch the input file first
    _input_file = input_file
    if fetch_input and input_is_remote:
        print_section("fetch input")
        # assume the INFN redirector
        remote_file = input_file if input_scheme else f"root://xrootd-cms.infn.it/{input_file}"
        _input_file = "file://" + os.path.join(tmp_dir, "input.root")
        cmd = f"xrdcp {remote_file} {_input_file}"
        print(f"cmd: {cmd}")
        if dry_run:
            print("dry run, skip")
        else:
            print("")
            p, _ = run_command(cmd)
            if p.returncode:
                raise Exception("xrdcp failed")

    # build the cmsRun command
    cmd = "cmsRun"
    cwd = os.getcwd()
    if change_dir:
        cwd = analyzer_dir
        cmd += f" {os.path.basename(analyzer_path)}"
    else:
        cmd += f" {analyzer_path}"
    cmd += f" inputFiles={_input_file}"
    cmd += f" outputFile=file://{_output_file}"
    if analyzer_options:
        cmd += " " + " ".join(analyzer_options)

    # run it
    print_section("run analyzer")
    print(f"cwd: {cwd}")
    print(f"cmd: {cmd}")
    if dry_run:
        print("dry run, skip")
    else:
        print("")
        p, _ = run_command(cmd, cwd=cwd)
        if p.returncode:
            raise Exception("cmdRun failed")

    # stageout
    print_section("stage out")
    cmd = f"gfal-copy file://{_output_file} {output_file}"
    print(f"cmd: {cmd}")
    if dry_run:
        print("dry run, skip")
    else:
        print("")
        p, _ = run_command(cmd)
        if p.returncode:
            raise Exception("stageout failed")

    # cleanup the tmp dir
    if os.path.exists(tmp_dir):
        try:
            shutil.rmtree(tmp_dir)
        except BaseException:
            pass


def get_scheme(path):
    m = re.match(r"^(\w+)\:\/\/.*$", path)
    return m.group(1) if m else None


def print_section(name):
    print("")
    print(f" {name} ".center(100, "-"))


def main():
    import argparse

    # setup arguments
    parser = argparse.ArgumentParser(description="run a local job and perform output stage out")
    parser.add_argument("input_file", help="the input file to process")
    parser.add_argument("output_file", help="the full location of the output file, e.g. "
        "'root://...', 'gsiftp://', or '/path/to/local/file'")
    parser.add_argument("options", nargs="*", help="additional analyzer options")
    parser.add_argument("--analyzer", "-a", default=os.path.join(thisdir, "analyzer.py"),
        help="the analyzer config")
    parser.add_argument("--change-dir", "-c", action="store_true", help="change into the directory "
        "of the analyzer")
    parser.add_argument("--fetch-input", "-f", action="store_true", help="fetch the input file "
        "first before running")
    parser.add_argument("--dry-run", "-d", action="store_true", help="perform a dry run")
    args = parser.parse_args()

    # start
    run_and_stageout(args.input_file, args.output_file, analyzer_path=args.analyzer,
        analyzer_options=args.options, change_dir=args.change_dir, fetch_input=args.fetch_input,
        dry_run=args.dry_run)


if __name__ == "__main__":
    main()