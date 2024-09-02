# - checkMesh (pyTree) -
import Generator.PyTree as G
import Connector.PyTree as X
import Geom.PyTree as D
import KCore.test as test

t = D.sphere6((0,0,0),2,20)
t = X.connectMatch(t, tol=1.e-6, dim=2)
infos = G.checkMesh(t, addGC=True)
test.testO(infos, 1)
