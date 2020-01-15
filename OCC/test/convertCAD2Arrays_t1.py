# - convertCAD2Arrrays (arrays) -
import Converter as C
import OCC
import KCore.test as T

#~ coarse =  O.convertCAD2Arrays("hammer.iges", h=500., chordal_err=8.)
#~ C.convertArrays2File(coarse, 'hammer_coarsec.plt')

A = OCC.convertCAD2Arrays("hammer.iges", h=0., chordal_err=0.,algo=0)
T.testA(A,1)

#~ fine = OCC.convertCAD2Arrays("hammer.iges", h=20., chordal_err=0.4)
#~ C.convertArrays2File(fine, 'hammer_fine.plt')

A = OCC.convertCAD2Arrays("hammer.iges", deflection=1. ,algo=1)
T.testA(A,2)
