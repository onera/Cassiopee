# - checkMesh (pyTree) -
import Generator.PyTree as G
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (10,10,10))
infos = G.checkMesh(a)
test.testO(infos, 1)
