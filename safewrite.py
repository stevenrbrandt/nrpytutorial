from io import StringIO
import os
from colored import colored

class SafeWrite:
    def __init__(self, fname):
        self.fname = fname
        self.fd = None
    def __enter__(self):
        self.fd = StringIO()
        return self.fd
    def __exit__(self, ty, val, tb):
        newcontent = self.fd.getvalue()
        if os.path.exists(self.fname):
            with open(self.fname) as fd:
                oldcontent = fd.read()
            do_write = newcontent.strip() != oldcontent.strip()
            if do_write:
                print("<<",oldcontent.strip(),">>",sep='')
                print("<<",newcontent.strip(),">>",sep='')
        else:
            do_write = True
        if do_write:
            print("Write:",colored(self.fname,"cyan"))
            with open(self.fname, "w") as fd:
                fd.write(newcontent)
                raise Exception()
