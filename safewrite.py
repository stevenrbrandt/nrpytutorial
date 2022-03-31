from io import StringIO
import os

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
        else:
            do_write = True
        if do_write:
            print("Write:",self.fname)
            with open(self.fname, "w") as fd:
                fd.write(newcontent)
