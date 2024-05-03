#!/usr/bin/env python
# remplace recursivement expression par replacement
# tree_replace . expression replacement

import os
import os.path as OP
import re

def subfunction(args, dir, file):
  arg0 = args[0]; arg1 = args[1]; arg2 = args[2]
  for f in file:
    (root,ext) = OP.splitext(f)
    tot = '%s/%s'%(dir,f)
    t = OP.islink(tot)
    if (ext in [".C",".h",".hxx",".hpp",".for",".f90",".html",".htm",".c",".cpp",".cxx",".py",".scons",".pyx",".pxi",".pxd",".rst"] and not t):

      f1 = open("%s/%s"%(dir,f), 'rb')
      l1 = f1.readlines()
      f1.close()
      
      i = 0
      modified = False
      for l in l1:
        if arg0.search(l) is not None:
          l1[i] = l.replace(arg1, arg2)
          modified = True
        i += 1
      if modified:
        f1 = open("%s/%s"%(dir,f),'wb')
        f1.writelines(l1)
        f1.close()

def walkOn(rootdir, myre, string, replaceString):
  try: # python 3
    for root, dirs, files in os.walk(rootdir):
      subfunction([myre, string, replaceString], root, files)
    
  except: # python 2
    OP.walk(rootdir, subfunction, [myre, string, replaceString])

if __name__ == "__main__":
  import sys
  exp = bytes(sys.argv[2], encoding='utf-8')
  myre = re.compile(exp)
  walkOn(sys.argv[1], myre, exp, bytes(sys.argv[3], encoding='utf-8'))
