# - convertCAD2PyTree (PyTree) -
import Converter.PyTree as C
import OCC.PyTree as OCC
import KCore.test as test

# fmt_iges
#~ coarse =  OCC.convertCAD2PyTree("hammer.iges", h=500., chordal_err=8.)
#~ C.convertPyTree2File(coarse, 'hammer_coarsec.plt')

default = OCC.convertCAD2PyTree("hammer.iges", h=0., chordal_err=0., algo=0)
test.testT(default,1)

#~ fine = OCC.convertCAD2PyTree("hammer.iges", h=20., chordal_err=0.4)
#~ C.convertPyTree2File(fine, 'hammer_fine.plt')
