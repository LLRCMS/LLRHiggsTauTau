# coding: utf-8

import os
import sys
import re
import subprocess


def get_voms_proxy_file():
    """
    Returns the path to the voms proxy file.
    """
    default_file = "/tmp/x509up_u{}".format(os.getuid())
    return os.getenv("X509_USER_PROXY", default_file)


def _voms_proxy_info(args=None, proxy_file=None, silent=False):
    cmd = ["voms-proxy-info"] + (args or [])

    # when proxy_file is None, get the default, then add it
    if proxy_file is None:
        proxy_file = get_voms_proxy_file()
    if proxy_file:
        proxy_file = os.path.expandvars(os.path.expanduser(proxy_file))
        cmd.extend(["--file", proxy_file])

    try:
        out = subprocess.check_output(" ".join(cmd), shell=True, executable="/bin/bash",
            stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as err:
        raise Exception("voms-proxy-info failed: {}".format(err))

    return out


def get_voms_proxy_user(proxy_file=None):
    """
    Returns the owner of the voms proxy. When *proxy_file* is *None*, it defaults to the result of
    :py:func:`get_voms_proxy_file`. Otherwise, when it evaluates to *False*, ``voms-proxy-info`` is
    queried without a custom proxy file.
    """
    out = _voms_proxy_info(args=["--identity"], proxy_file=proxy_file).strip()
    cns = re.findall(r"/CN=[a-z]+", out.replace(" ", ""))
    if not cns:
        raise Exception("no valid identity found in voms proxy: {}".format(out))

    # extract actual names
    names = [cn[4:] for cn in cns]

    # return the shortest name
    return sorted(names, key=len)[0]


class TeeStream(object):
    """
    Stream class that writes simultaneously into a file and stdout.
    """

    def __init__(self, log_file, mode="w"):
        super(TeeStream, self).__init__()

        self.log_file = log_file
        self.stdout = sys.stdout
        self.mode = mode
        self.f = None
        self.flush_on_write = True

    @property
    def closed(self):
        return self.f is None

    def __del__(self):
        self.close()

    def __enter__(self):
        self.open()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def close(self):
        if self.closed:
            return
        self.flush()
        self.f.close()
        self.f = None

    def open(self):
        if not self.closed:
            return
        self.f = open(self.log_file, self.mode)

    def flush(self):
        if self.closed:
            return
        self.f.flush()
        self.stdout.flush()

    def write(self, s):
        if self.closed:
            return
        line = "{}\n".format(s)
        self.f.write(line)
        self.stdout.write(line)
        if self.flush_on_write:
            self.flush()
