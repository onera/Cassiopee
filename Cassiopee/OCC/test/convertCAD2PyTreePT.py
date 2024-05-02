# - convertIGES2PyTree (PyTree) -
import OCC.PyTree as OCC
import Converter.PyTree as C

t = OCC.convertCAD2PyTree('hammer.iges', h=0., chordal_err=0., algo=1)
C.convertPyTree2File(t, 'out.cgns')
