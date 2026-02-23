# - convertCAD2PyTree (PyTree) -
import OCC.PyTree as OCC
import KCore.test as test

default = OCC.convertCAD2PyTree("hammer.iges", chordal_err=1., algo=0)
test.testT(default,1)
