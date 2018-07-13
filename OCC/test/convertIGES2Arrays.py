# - convertIGES2Arrays (IGES file format) -
import Converter as C
import OCC

A = OCC.convertIGES2Arrays('hammer.iges', h=0., chordal_err=0.)
C.convertArrays2File(A, 'hammer.plt')
