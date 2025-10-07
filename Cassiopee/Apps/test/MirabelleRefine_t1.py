# - Mesh.Mirabelle3 -
# propagate a refinement
import Converter.PyTree as C
import Transform.PyTree as T
import Generator.PyTree as G
import Connector.PyTree as X
import Geom.PyTree as D
import Converter.Internal as Internal
import Apps.Mesh.Mirabelle3 as Mirabelle
import KCore.test as test

#a = D.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (50,50,30))
a = D.cylinder((0,0,0), 1., 10.,N=20)
distrib = G.cart((0.,0.,0.),(0.05,1,1),(11,1,1))
a = G.addNormalLayers(a, distrib,niter=50)
for z in Internal.getZones(a):
    C._addBC2Zone(z, "wall", "BCWall", 'kmin')
    C._addBC2Zone(z, "ovst", "BCOverlap", 'kmax')

t = C.newPyTree(['Base',Internal.getZones(a)])
t = X.connectMatch(t)

zonesToRefine = ["cyl-part11"]
dirs = [1]

Mirabelle._refine(t, zonesToRefine, dirs, refined={}, factor=2)
test.testT(t, 1)