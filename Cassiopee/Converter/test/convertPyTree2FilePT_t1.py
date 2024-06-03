# - convertPyTree2File (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

LOCAL = test.getLocal()

# Struct - OK
z = G.cart((0,0,0), (1,1,1), (10,10,10))
C.convertPyTree2File(z, LOCAL+'/out1.tp')
test.testF(LOCAL+'/out1.tp', 1)

# BE - OK
z = G.cartHexa((0,0,0), (1,1,1), (10,10,10))
C.convertPyTree2File(z, LOCAL+'/out2.tp')
test.testF(LOCAL+'/out2.tp', 2)

# ME - output HEXA
a = G.cartHexa((0,0,0), (1,1,1), (10,10,10)) 
b = G.cartTetra((9,0,0), (1,1,1), (10,10,10))
z = C.mergeConnectivity(a, b, boundary=0)
C.convertPyTree2File(z, 'out3.tp')
test.testF(LOCAL+'/out3.tp', 3)

# NGON3
z = G.cartNGon((0,0,0), (1,1,1), (10,10,10), api=1)
C.convertPyTree2File(z, 'out4.tp')
test.testF(LOCAL+'/out4.tp', 4)

# NGON4
z = G.cartNGon((0,0,0), (1,1,1), (10,10,10), api=3)
C.convertPyTree2File(z, 'out5.tp')
test.testF(LOCAL+'/out5.tp', 5)
