# - bbox (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as test

# bbox d une zone
a = G.cart((0.,0.,0.),(0.1,0.1,1.),(20,20,20))
a = C.addBC2Zone(a, 'wall1','BCWall','jmin')
a = C.addBC2Zone(a, 'match1','BCMatch','imin',a,'imax',[1,2,3])
a = C.addBC2Zone(a, 'match2','BCMatch','imax',a,'imin',[1,2,3])
sol = G.bbox(a); test.testO(sol,1)

# bbox d un arbre
t = C.newPyTree(['Base'])
b = G.cart((0.5,0.,0.),(0.1,0.1,1.),(20,20,20)); b[0] = 'cart2'
t[2][1][2] += [a,b]
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
sol = G.bbox(t); test.testO(sol,2)
