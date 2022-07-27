from io import StringIO
import os
from difflib import context_diff
import sys
from subprocess import Popen, PIPE
from clang_format import _get_executable as get_executable

clang_formatter = get_executable("clang-format")

verbose = False
nochange = False

# ANSI colors
red = "\033[31;1m"
green = "\033[32;1m"
reset = "\033[0m"

class SafeWrite:
    def __init__(self, fname, do_format=False):
        self.fname = fname
        self.fd = None
        self.do_format = do_format
    def __enter__(self):
        self.fd = StringIO()
        return self.fd
    def __exit__(self, ty, val, tb):
        print("Checking",self.fname,end="...")
        newcontent = self.fd.getvalue()
        if self.do_format:
            pipe = Popen([clang_formatter],stdout=PIPE,stdin=PIPE,universal_newlines=True)
            out, err = pipe.communicate(newcontent)
            newcontent = out
        if os.path.exists(self.fname):
            with open(self.fname) as fd:
                oldcontent = fd.read()
            do_write = newcontent.strip() != oldcontent.strip()
            if do_write and verbose:
                print("Diff for:",self.fname)
                oldlines=[line+"\n" for line in oldcontent.strip().split("\n")]
                newlines=[line+"\n" for line in newcontent.strip().split("\n")]
                sys.stdout.writelines(context_diff(oldlines,newlines,fromfile='before',tofile='after'))
        else:
            do_write = True
        if do_write:
            assert nochange == False
            print("Write:",self.fname)
            with open(self.fname, "w") as fd:
                fd.write(newcontent)
            print(" "+red+"[written]"+reset)
        else:
            print(" "+green+"[no changes]"+reset)
