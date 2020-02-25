# - convertCAD2Arrrays (arrays) -
import Converter as C
import OCC
import KCore.test as test

#~ coarse =  O.convertCAD2Arrays("hammer.iges", h=500., chordal_err=8.)
#~ C.convertArrays2File(coarse, 'hammer_coarsec.plt')

#~ fine = OCC.convertCAD2Arrays("hammer.iges", h=20., chordal_err=0.4)
#~ C.convertArrays2File(fine, 'hammer_fine.plt')

A = OCC.convertCAD2Arrays("hammer.iges", h=0., chordal_err=0., algo=0)
test.testA(A,1)

A = OCC.convertCAD2Arrays("hammer.iges", deflection=1., algo=1)
test.testA(A,2)
