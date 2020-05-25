"""OpenCascade definition module.
"""
__version__ = '3.1'
__author__ = "Sam Landier"

from . import occ

def convertCAD2Arrays(fileName, format='fmt_iges', 
                      h=0., chordal_err=0., growth_ratio=0., algo=1):
    """Convert a CAD (IGES or STEP) file to arrays.
    Usage: a = convertCAD2Arrays(fileName, options)"""
    if algo == 0: # pure OCC
        if chordal_err == 0.: chordal_err = 1.
        return occ.convertCAD2Arrays0(fileName, format, "None", "None", chordal_err)
    elif algo == 1: # OCC+T3Mesher
    	return  occ.convertCAD2Arrays1(fileName, format, h, chordal_err, growth_ratio)
    else: # OCC+T3Mesher v2
    	return  occ.convertCAD2Arrays2(fileName, format, h, chordal_err, growth_ratio)