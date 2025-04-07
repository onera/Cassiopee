# - mergeConnectivity (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cartHexa((0,0,0), (1,1,1), (10,10,10))
b = G.cartHexa((0,0,0), (1,1,1), (10,10,1))
c = C.mergeConnectivity(a, b, boundary=1)
# c contains now a volume HEXA connectivity and a QUAD boundary connectivity.
C.convertPyTree2File(c, 'out0.cgns')

a = G.cartHexa((0,0,0), (1,1,1), (10,10,10))
b = G.cartTetra((0,0,0), (1,1,1), (10,10,10))
c = C.mergeConnectivity(a, b, boundary=0)
# c is now a multiple-element zone containing a volume HEXA connectivity and
# a volume TETRA connectivity.

C.convertPyTree2File(c, 'out.cgns')
