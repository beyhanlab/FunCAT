#!/usr/bin/env python3

from os import chdir
import os.path
from tempfile import mkdtemp
from shutil import rmtree

class TempDir:
    """Convenience class for creating and implicitly cleaning up a temp dir

    Usage:
    
    with TempDir() as tmp_dir:
       # do stuff relative to tmp_dir
    """
    def __init__(self, move_to_tmpdir = False):
        self.move_to_tmpdir = move_to_tmpdir
    def __enter__(self):
        self.tmp_dir = mkdtemp()
        if(self.move_to_tmpdir):
            self.olddir = os.path.abspath(os.curdir)
            chdir(self.tmp_dir)
        return self.tmp_dir
    def __exit__(self, type, value, traceback):
        if(self.move_to_tmpdir):
            chdir(self.olddir)
        rmtree(self.tmp_dir)

if(__name__ == "__main__"):
    import os.path
    with TempDir() as temp_dir:
        fname = os.path.join(temp_dir,"test")
        s = "Hello, world\n"
        open(fname,"w").write(s)
        assert(open(fname).read() == s)
        print("Passed")

