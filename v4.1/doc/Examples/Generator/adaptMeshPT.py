# - adaptMesh (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C

# HEXA but no adaptation in the 3rd direction (2D adaptation)
a = G.cartHexa((0,0,0),(0.1,0.1,0.1),(11,11,2))
C._fillEmptyBCWith(a, 'nref', 'BCFarfield', dim=2)
C._initVars(a,'{centers:indicator}=({centers:CoordinateX})>0.5')
a = G.adaptMesh(a, indicator="indicator", hook=None, dim=2, conformize=False)
C.convertPyTree2File(a, "out.cgns")
