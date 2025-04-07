# - convertCAD2Arrays (arrays) -
import OCC
import KCore.test as test

A = OCC.convertCAD2Arrays("hammer.iges", chordal_err=1., algo=0)
test.testA(A,1)

A = OCC.convertCAD2Arrays("hammer.iges", h=0., chordal_err=0., algo=1)
test.testA(A,2)

A = OCC.convertCAD2Arrays("hammer.iges", h=0., chordal_err=0., algo=2)
test.testA(A,3)
