# - convertBAR2Struct (array) -
import Converter as C
import Generator as G
import Geom as D
import Transform as T
import Post as P
import KCore.test as test
def F(x,y): return 3*x + 4*y
# cas simple
a = D.circle((0.,0.,0.),1.)
a = C.convertArray2Hexa(a); a = G.close(a)
a = C.initVars(a,'F', F, ['x','y'])
b = C.convertBAR2Struct(a)
test.testA([b], 1)
#
# cas ou la bar n'est pas ordonnee - pas de boucle
#
a1 = G.cart((0.,0.,0.),(0.1,0.1,0.1),(11,11,1))
a2 = G.cart((0.,1.,0.),(0.1,0.1,0.1),(11,1,11))
a2 = T.rotate(a2,(0.,1.,0.),(1,0,0),-55)
A = [a1,a2]; A = C.convertArray2Hexa(A); a = T.join(A)
a = T.reorder(a,(1,))
c = P.sharpEdges(a,30.)[0]; c = G.close(c)
c = C.convertBAR2Struct(c)
test.testA([c],2)
