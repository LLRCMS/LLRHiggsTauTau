# coding: utf-8

import os
import sys
import hashlib

try:
    from CRABClient.UserUtilities import getUsername as _crab_get_username
except ImportError:
    _crab_get_username = None


def create_hash(inp, l=10, algo="sha256", to_int=False):
    h = getattr(hashlib, algo)(str(inp)).hexdigest()[:l]
    return int(h, 16) if to_int else h


def hash_truncate(s, max_length, hash_length=8, **kwargs):
    if hash_length > max_length:
        raise ValueError("the hash length must not be smaller than the maximum string length")
    if len(s) > max_length:
        s = s[:max_length - hash_length] + create_hash(s, l=hash_length, **kwargs)
    return s


def get_wlcg_user(*args, **kwargs):
    if callable(_crab_get_username):
        return _crab_get_username(*args, **kwargs)

    if "WLCG_USER" in os.environ:
        return os.environ["WLCG_USER"]

    raise Exception("cannot determine wlcg user name; either source crab or set WLCG_USER")


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

    def fileno(self):
        return self.stdout.fileno()

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
