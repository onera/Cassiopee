#!/usr/bin/env python
# Usage: wc . 
# Compte toutes les lignes dans les C/py/.. files from .
# Usage: wc.py .
import os.path as OP
import re

regBuild = re.compile('build')

def subfunction(args, dir, file):
  nlines = args[0]
  for f in file:
    (root,ext) = OP.splitext(f)
    tot = '%s/%s'%(dir,f)
    t = OP.islink(tot)
    m = True
    if regBuild.search(dir) is not None: m = False
    if (ext in [".C",".h",".hpp",".for",".f90",".html",".htm",".c",".cpp",".cxx",".py",".scons",".pyx",".rst"] and not t and m):
      f1 = open("%s/%s"%(dir,f),'r')
      l1 = f1.readlines()
      f1.close()
      nlines[0] += len(l1)

def walkOn(rootdir, nlines):
  OP.walk(rootdir, subfunction, [nlines])

if __name__ == "__main__":
  import sys
  nlines = [0]
  walkOn(sys.argv[1], nlines)
  print(nlines[0])
