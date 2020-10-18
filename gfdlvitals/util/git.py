""" Tools to work with Git Repository """

import os
import subprocess
import sys

__ALL__ = ["retrieve_commit", "status"]


def is_clean(path):
    """ Determines if git repo is clean """
    return_code = status(path)
    return bool(return_code == 0)


def retrieve_commit(path):
    """ Obtain current commit """
    cwd = os.getcwd()
    os.chdir(path)
    cmd = "git rev-parse HEAD"
    commit = subprocess.check_output(cmd.split(" "))
    commit = str(commit.decode("utf-8")).replace("\n", "")
    os.chdir(cwd)
    return commit


def status(path):
    """ Gets the current status of the local repo """
    cwd = os.getcwd()
    os.chdir(path)
    cmd = "git diff-index --quiet HEAD"
    result = subprocess.run(cmd.split(" "), stdout=subprocess.DEVNULL, check=False)
    os.chdir(cwd)
    return result.returncode


if __name__ == "__main__":
    if len(sys.argv) > 1:
        _path = os.path.abspath(sys.argv[1])
    else:
        _path = os.getcwd()
    COMMIT = retrieve_commit(_path)
    print(COMMIT)
    _return_code = status(_path)
    print(_return_code)
    _clean = is_clean(_path)
    print(_clean)
