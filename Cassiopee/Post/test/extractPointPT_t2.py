# - extractPoint (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P
import KCore.test as test

# Maillage en noeuds
ni = 11; nj = 11; nk = 1
a = G.cart((0,0,0), (1./(ni-1),1./(nj-1),1.), (ni,nj,nk))
a = C.addBC2Zone(a,'overlap','BCOverlap','imin')
a = C.fillEmptyBCWith(a, 'nref','BCFarfield')

# Create a function
def F(x,y,z): return 2*x*x*x*x*x + 2.*y*y*z + z*z
val0 = F(0.55,0.38,0.) # reference
# init by function
a = C.initVars(a, 'F', F, ['CoordinateX','CoordinateY','CoordinateZ'])
val = P.extractPoint(a, (0.55, 0.38, 0.))
test.testO(val)

# sur un arbre
t = C.newPyTree(['Base',2]); t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t[2][1][2].append(a)
t = C.initVars(t, 'F', F, ['CoordinateX','CoordinateY','CoordinateZ'])
t = C.initVars(t, 'centers:G', 3.)
val = P.extractPoint(t, (0.55, 0.38, 0.))
test.testO(val,2)

# non structure
a = C.convertArray2Tetra(a)
t = C.newPyTree(['Base',2]); t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t[2][1][2].append(a)
t = C.initVars(t, 'F', F, ['CoordinateX','CoordinateY','CoordinateZ'])
val = P.extractPoint(t, (0.55, 0.38, 0.))
test.testO(val,3)
