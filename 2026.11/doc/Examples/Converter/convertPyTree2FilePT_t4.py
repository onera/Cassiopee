# - convertPyTree2File (pyTree) -
# - binary tecplot -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

LOCAL = test.getLocal()

# Struct
z = G.cart((0,0,0), (1,1,1), (10,10,10))
C.convertPyTree2File(z, LOCAL+'/out1.plt')
test.testF(LOCAL+'/out1.plt', 1)

# BE
z = G.cartHexa((0,0,0), (1,1,1), (10,10,10))
C.convertPyTree2File(z, LOCAL+'/out2.plt')
test.testF(LOCAL+'/out2.plt', 2)

# ME - output HEXA
a = G.cartHexa((0,0,0), (1,1,1), (10,10,10))
b = G.cartTetra((9,0,0), (1,1,1), (10,10,10))
z = C.mergeConnectivity(a, b, boundary=0)
C.convertPyTree2File(z, LOCAL+'/out3.plt')
test.testF(LOCAL+'/out3.plt', 3)

# NGON3
z = G.cartNGon((0,0,0), (1,1,1), (10,10,10), api=1)
C.convertPyTree2File(z, LOCAL+'/out4.plt')
test.testF(LOCAL+'/out4.plt', 4)

# NGON4
z = G.cartNGon((0,0,0), (1,1,1), (10,10,10), api=3)
C.convertPyTree2File(z, LOCAL+'/out5.plt')
test.testF(LOCAL+'/out5.plt', 5)
