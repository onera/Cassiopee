"""OpenCascade definition module.
"""
__version__ = '3.1'
__author__ = "Sam Landier"

from . import occ

def convertCAD2Arrays(fileName, format='fmt_iges', 
                      h=0., chordal_err=0., growth_ratio=0., 
                      deflection=1., algo=0):
    """Read file and create arrays containing file data.
    Usage: a = convertCAD2Arrays(fileName, options)"""
    if algo == 0: return  occ.convertCAD2Arrays1(fileName, format, h, chordal_err, growth_ratio)
    else: return occ.convertCAD2Arrays2(fileName, format, "None", "None", deflection)
