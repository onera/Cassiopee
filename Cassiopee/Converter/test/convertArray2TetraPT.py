# - convertArray2Tetra (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

# 2D : triangles
a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
b = C.convertArray2Tetra(a)

# 3D : tetrahedras
a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
b = C.convertArray2Tetra(a)
C.convertPyTree2File([a,b], 'out.cgns')
