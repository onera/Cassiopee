# - convertIGES2PyTree (IGES file format) -
import OCC.PyTree as OCC
import Converter.PyTree as C

t = OCC.convertIGES2PyTree('hammer.iges', h=0., chordal_err=0.)
C.convertPyTree2File(t, 'out.cgns')
