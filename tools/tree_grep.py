# Usage: tree_grep . toto
# Find "toto" string in all C/py/.. files from .
# Usage: tree_grep . " ; "
# Find " ; " string in all files from .
import os.path as OP
import os
import re

regBuild = re.compile('build')

def subfunction(args, dir, file):
  regexp = args[0]
  exts = args[1]
  if exts is None:
    exts = [".C",".h",".hpp",".hxx",".F",".for",".f90",".html",".htm",".c",".cpp",".cxx",".py",".scons",".pyx",".pxi",".pxd",".rst"]
  for f in file:
    (root,ext) = OP.splitext(f)
    tot = '%s/%s'%(dir,f)
    t = OP.islink(tot)
    m = True
    if regBuild.search(dir) is not None: m = False
    if ext in exts and not t and m:
      
      f1 = open("%s/%s"%(dir,f), 'rb')
      l1 = f1.readlines()
      f1.close()
      
      a = 1
      for l in l1:
        if re.search(regexp, l) is not None:
          if isinstance(l, bytes):
            try: l = l.decode()
            except: pass
          print(tot, "line", a, l)
        a += 1
        
def walkOn(rootdir, myre, exts):
  try: # python 3
    for root, dirs, files in os.walk(rootdir):
        subfunction([myre, exts], root, files)
    
  except: # python 2
    OP.walk(rootdir, subfunction, [myre, exts])

if __name__ == "__main__":
  import sys
  exts = None
  if len(sys.argv) == 4: 
    if sys.argv[3] != '':
      exts = sys.argv[3]
      if exts[0] == '*': exts = exts[1:]
      exts = [exts]
  exp = bytes(sys.argv[2], encoding='utf-8')
  #exp = bytes(sys.argv[2])
  myre = re.compile(exp)
  walkOn(sys.argv[1], myre, exts)
