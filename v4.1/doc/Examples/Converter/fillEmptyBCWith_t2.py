# - fillEmptyBC (pyTree) -
# avec des raccords partiels
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

a = G.cart((0.,0.,0.),(0.1,0.1,0.1),(10,10,10))
a = C.addBC2Zone(a,'wall1','BCWall',[1,1,5,7,1,10])
a = C.addBC2Zone(a,'wall2','BCWall',[10,10,5,7,1,10])
a = C.addBC2Zone(a,'wall3','BCWall',[1,4,1,1,1,10])
a = C.addBC2Zone(a,'wall4','BCWall',[1,4,10,10,1,10])
a = C.fillEmptyBCWith(a,'wallv','BCWallViscous')
t = C.newPyTree(['Base',a])
test.testT(t)
