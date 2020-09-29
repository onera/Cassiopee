"""OpenCascade definition module.
"""
__version__ = '3.2'
__author__ = "Sam Landier"

from . import occ

# algo=0: mailleur open cascade (chordal_error)
# algo=1: algorithme T3mesher (h, chordal_error, growth_ratio)
# algo=2: algorithme T3mesher (h, chordal_error, growth_ratio, merge_tol)
def convertCAD2Arrays(fileName, format='fmt_iges', 
                      h=0., chordal_err=0., growth_ratio=0., 
                      merge_tol=-1, algo=1):
    """Convert a CAD (IGES or STEP) file to arrays.
    Usage: a = convertCAD2Arrays(fileName, options)"""
    if algo == 0: # pure OCC
        if chordal_err == 0.: chordal_err = 1.
        a = occ.convertCAD2Arrays0(fileName, format, "None", "None", chordal_err)
        try: import Generator; a = Generator.close(a)
        except: pass
        return a
    elif algo == 1: # OCC+T3Mesher
    	return  occ.convertCAD2Arrays1(fileName, format, h, chordal_err, growth_ratio)
    else: # OCC+T3Mesher v2
    	return  occ.convertCAD2Arrays2(fileName, format, h, chordal_err, growth_ratio, merge_tol)