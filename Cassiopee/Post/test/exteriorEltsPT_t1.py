# - exteriorElts (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G
import KCore.test as test

a = G.cartTetra((0,0,0), (1,1,1), (10,10,10))
C._addVars(a,'Density'); C._addVars(a,'centers:cellN')
b = P.exteriorElts(a)
C.convertPyTree2File(b, "gg.cgns")
test.testT(b, 1)

"""
a = G.cartNGon((0,0,0), (1,1,1), (10,10,10), api=1)
C._addVars(a,'Density'); C._addVars(a,'centers:cellN')
b = P.exteriorElts(a)
test.testT(b, 10)
"""
