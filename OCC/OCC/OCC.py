"""OpenCascade definition module.
"""
__version__ = '3.1'
__author__ = "Sam Landier"

from . import occ

def convertIGES2Arrays(fileName, h=0., chordal_err=0., growth_ratio=0.):
    """Read file and create arrays containing file data.
    Usage: a = convertIGES2Arrays(fileName, options)"""
    try: file = open(fileName, 'r')
    except: raise IOError("convertIGES2Arrays: file not found.")
    file.close()
    import os.path
    ext = os.path.splitext(fileName)
    extension = ext[len(ext)-1]
    if extension == '.iges' or extension == '.igs' or extension == '.IGES' or extension == '.IGS':
      fmt = "iges"
    elif extension == '.step' or extension == '.stp' or extension == '.STEP' or extension == '.STP':
      fmt = "step"
    else:
      raise IOError("convertIGES2Arrays: file with extension %s is not supported."%extension)
      return
        
    try:
      return  occ.convertIGES2Arrays(fileName, fmt, h, chordal_err, growth_ratio)
    except: pass
