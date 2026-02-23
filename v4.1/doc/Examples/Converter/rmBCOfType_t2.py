# - rmBCOfType (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

# Connectivite sur une zone
a1 = G.cart((0.,0.,0.),(0.1,0.1,0.1),(11,11,11))
a1 = C.initVars(a1,'F',1.); a1 = C.initVars(a1,'centers:G',2.)
a1 = C.addBC2Zone(a1,'nearmatch1','BCNearMatch','imax',a1,'imin',[1,2,3])
a1 = C.rmBCOfType(a1,'BCNearMatch')
test.testT(a1,1)
a1 = C.addBC2Zone(a1,'match1','BCMatch','imax',a1,'imin',[1,2,3])
a1 = C.rmBCOfType(a1,'BCMatch')
test.testT(a1,2)
a1 = C.addBC2Zone(a1,'overlap1','BCOverlap','imax')
a1 = C.rmBCOfType(a1,'BCOverlap')
test.testT(a1,3)

# sur un arbre
a1 = G.cart((0.,0.,0.),(0.1,0.1,0.1),(11,11,11)); a1[0] = 'cart1'
a2 = G.cart((1.,0.,0.),(0.2,0.2,0.2),(6,6,6)); a2[0] = 'cart2'
a1 = C.addBC2Zone(a1,'nearmatch1','BCNearMatch','imax',a2,'imin',[1,2,3])
a2 = C.addBC2Zone(a2,'nearmatch2','BCNearMatch','imin',a1,'imax',[1,2,3])
a3 = G.cart((-1.,0.,0.),(0.1,0.1,0.1),(11,11,11)); a3[0] = 'cart3'
a1 = C.addBC2Zone(a1,'match1','BCMatch','imin',a3,'imax',[1,2,3])
a3 = C.addBC2Zone(a3,'match1','BCMatch','imax',a1,'imin',[1,2,3])
t = C.newPyTree(['Base']); t[2][1][2] += [a1,a2]
t = C.initVars(t, 'F', 1.); t = C.initVars(t, 'centers:G', 2.)
t = C.fillEmptyBCWith(t,'ov','BCOverlap')
t = C.rmBCOfType(t,'BCNearMatch')
test.testT(t,4)
t = C.rmBCOfType(t,'BCMatch')
test.testT(t,5)
t = C.rmBCOfType(t,'BCOverlap')
test.testT(t,6)
