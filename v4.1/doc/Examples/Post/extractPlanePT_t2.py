# - extractPlane (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P
import KCore.test as test

def f(u,v,w): return 3*u*v + 2*w

# test structure
m = G.cylinder((0.,0.,0.), 0., 1., 0, 360, 1., (50,20,2))
m = C.addBC2Zone(m,'overlap','BCOverlap','imin')
m = C.fillEmptyBCWith(m, 'nref','BCFarfield')
m = C.initVars(m, 'F', f, ['CoordinateX','CoordinateY','CoordinateZ'])
m = C.initVars(m, 'centers:G', 3.)
#m = T.subzone(m,(1,20,1),(50,20,2)); m = G.close(m) : pb : voir extractPlanePT_t1.py
t = C.newPyTree(['Base']); t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t[2][1][2].append(m)
t2 = P.extractPlane(t,(0.,0.,1.,-0.5))
test.testT(t2,1)

# test non structure tetra
m2 = C.convertArray2Tetra(m)
m2 = C.initVars(m2, 'F', f, ['CoordinateX','CoordinateY','CoordinateZ'])
m2 = C.initVars(m2, 'centers:G', 3.)
t = C.newPyTree(['Base']); t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t[2][1][2].append(m2)
t2 = P.extractPlane(t,(0.,0.,1.,-0.5))
test.testT(t2,2)

# test non structure hexa
m2 = C.convertArray2Hexa(m)
m2 = C.initVars(m2, 'F', f, ['CoordinateX','CoordinateY','CoordinateZ'])
m2 = C.initVars(m2, 'centers:G', 3.)
t = C.newPyTree(['Base']); t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t[2][1][2].append(m2)
t2 = P.extractPlane(t,(0.,0.,1.,-0.5))
test.testT(t2,3)
