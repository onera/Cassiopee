"""OpenCascade definition module.
"""
__version__ = '3.1'
__author__ = "Sam Landier"

from . import occ

def convertCAD2Arrays(fileName, format='fmt_iges', 
                      h=0., chordal_err=0., growth_ratio=0., algo=0):
    """Convert a CAD (IGES or STEP) file to arrays.
    Usage: a = convertCAD2Arrays(fileName, options)"""
    if algo == 0: 
    	return  occ.convertCAD2Arrays1(fileName, format, h, chordal_err, growth_ratio)
    else:
        if chordal_err == 0.: chordal_err = 1.
        return occ.convertCAD2Arrays2(fileName, format, "None", "None", chordal_err)
