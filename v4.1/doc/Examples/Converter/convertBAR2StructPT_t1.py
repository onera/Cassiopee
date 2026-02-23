# - convertBAR2Struct (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Geom.PyTree as D
import Transform.PyTree as T
import Post.PyTree as P
import KCore.test as test

# cas simple
a = D.circle((0.,0.,0.),1.)
a = C.convertArray2Hexa(a); a = G.close(a)
a = C.initVars(a, '{F}=3*{CoordinateX}+4*{CoordinateY}')
a = C.initVars(a, '{centers:D}={centers:CoordinateX}')
b = C.convertBAR2Struct(a)
t = C.newPyTree(['Base',1,b])
test.testT(t, 1)
#
# cas ou la bar n'est pas ordonnee - pas de boucle
#
a1 = G.cart((0.,0.,0.),(0.1,0.1,0.1),(11,11,1))
a2 = G.cart((0.,1.,0.),(0.1,0.1,0.1),(11,1,11))
a2 = T.rotate(a2,(0.,1.,0.),(1,0,0),-55)
A = [a1,a2]; A = C.convertArray2Hexa(A); a = T.join(A)
a = T.reorder(a,(1,))
a = P.sharpEdges(a,30.)[0];
a = G.close(a)
a = C.initVars(a,'{F}=3*{CoordinateX}+4*{CoordinateY}')
a = C.initVars(a,'centers:D',1.)
b = C.convertBAR2Struct(a)
t = C.newPyTree(['Base',1]);t[2][1][2] += [b]
test.testT(t,2)
