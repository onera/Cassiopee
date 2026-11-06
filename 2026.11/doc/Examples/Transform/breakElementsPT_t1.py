# - breakElements (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import KCore.test as test

a1 = G.cartNGon((-1.,0.,0.),(0.1,0.1,0.1),(11,3,3)) # 3D
a3 = G.cartNGon((-1.,0.,0.2),(0.1,0.1,0.1),(1,1,3)) # 1D
a4 = G.cartTetra((-0.4,0.2,0),(0.1,0.1,1),(5,5,1)) # 2D
a5 = G.cartHexa((-0.4,0.6,0),(0.1,0.1,1),(5,5,1)) # 2D
a6 = G.cartTetra((-0.4,1.,0),(0.1,0.1,0.1),(5,5,5)) # 3D
a7 = G.cartPenta((-0.4,1.4,0),(0.1,0.1,0.1),(5,5,5)) # 3D
A = [a1,a6,a7]
A = C.convertArray2NGon(A); a = T.join(A)
C._initVars(a,'{F}={CoordinateX}+{CoordinateY}**2')
# sur une zone
res = T.breakElements(a)
t = C.newPyTree(['Base', res]); t[2][1][2]+=res
test.testT(t,1)
# sur une liste de zones
res = T.breakElements([a])
t = C.newPyTree(['Base']); t[2][1][2]+=res
test.testT(t,2)
# sur une base
t = C.newPyTree(['Base']); t[2][1][2]+=[a]
res = T.breakElements(t[2][1])
t = C.newPyTree(['Base']); t[2][1][2]+=res
test.testT(t,3)
# sur un arbre
t = C.newPyTree(['Base']); t[2][1][2]+=[a]
res = T.breakElements(t)
t = C.newPyTree(['Base']); t[2][1][2]+=res
test.testT(t,4)
