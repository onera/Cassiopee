# - breakConnectivity (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cartHexa((0,0,0), (1,1,1), (10,10,10))
b = G.cartTetra((9,0,0), (1,1,1), (10,10,5))
c = C.mergeConnectivity(a, b, boundary=0)
# c is a zone with two connectivities (one HEXA and one TETRA)
t = C.newPyTree(['Base',c])
t = C.breakConnectivity(t)
# t contains now two zones (one pure HEXA, one pure TETRA)
C.convertPyTree2File(t, 'out.cgns')

# You can directly break a zone
A = C.breakConnectivity(c)
# A contains 2 zones
C.convertPyTree2File(A, 'out2.cgns')
