# - eikonal (array) -
import Dist2Walls
import Generator as G
import Converter as C
import Geom as D
import Connector as X
import KCore.test as test

# Bloc cartesien
N = 512
N = 128
h = 2./N
a = G.cart((-1.,-1.,-1.),(h,h,h),(N,N,N))
cellN = C.initVars(a, 'cellN', 1)
cellN = C.extractVars(cellN, ['cellN'])

# Init wall
sphere = D.sphere((0.,0.,0.), 0.1, 30)
sphere = C.convertArray2Tetra(sphere)
sphere = G.close(sphere)
cellN = X.blankCellsTri([a], [cellN], [sphere], blankingType=0)

a = C.addVars([a, cellN[0]])
a = C.initVars(a,'speed', 1./h)
a = C.initVars(a, '{Phi}=1.e12*({cellN}>0.)')

# Eikonal
b = Dist2Walls.eikonal(a)
test.testA([b],1)
