# - bboxOfCells (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as test

# STRUCT
z = G.cart((0.,0.,0.),(0.1,0.1,0.1),(10,10,10))
z1 = G.bboxOfCells(z)
test.testT(z1, 11)
G._bboxOfCells(z)
test.testT(z, 1)

# NGONv3
z = G.cartNGon((0.,0.,0.),(0.1,0.1,0.1),(10,10,10), api=1)
z1 = G.bboxOfCells(z)
test.testT(z1, 21)
G._bboxOfCells(z)
test.testT(z, 2)

# NGONv4
z = G.cartNGon((0.,0.,0.),(0.1,0.1,0.1),(10,10,10), api=3)
z1 = G.bboxOfCells(z)
test.testT(z1, 31)
G._bboxOfCells(z)
test.testT(z, 3)

# BE (HEXA)
z = G.cartHexa((0.,0.,0.),(0.1,0.1,0.1),(10,10,10))
z1 = G.bboxOfCells(z)
test.testT(z1, 41)
G._bboxOfCells(z)
test.testT(z, 4)

# ME (HEXA+TETRA)
a = G.cartHexa((0.,0.,0.),(0.1,0.1,0.1),(10,10,10))
b = G.cartTetra((0.9,0.,0.),(0.1,0.1,0.1),(10,10,10))
z = C.mergeConnectivity(a, b, boundary=0)
z1 = G.bboxOfCells(z)
test.testT(z1, 51)
G._bboxOfCells(z)
test.testT(z, 5)
