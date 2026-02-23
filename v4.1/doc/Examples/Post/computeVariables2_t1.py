# - computeVariables2 (array) -
import Converter as C
import Post as P
import Generator as G
import KCore.test as test

# Variables a calculer
# --------------------
vars = ['Pressure','VelocityX','VelocityZ','VelocityMagnitude',
        'Temperature','Entropy','Enthalpy','Mach','ViscosityMolecular',
        'PressureStagnation','TemperatureStagnation']

# test sur un array avec variables conservatives
# ----------------------------------------------
ni = 30; nj = 40
m  = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,2))

c1 = C.array('ro,rou,rov,row,roE', ni, nj, 2)
c1 = C.initVars(c1, 'ro', 1.5)
c1 = C.initVars(c1, 'rou', 2.)
c1 = C.initVars(c1, 'rov', 3.)
c1 = C.initVars(c1, 'row', 6.)
c1 = C.initVars(c1, 'roE', 1.e5)
a1 = C.addVars([m,c1])

P._computeVariables2(a1,vars,rgp=12.5)
test.testA([a1], 1)

# test sur un array avec variables primitives
# -------------------------------------------
c2 = C.array('ro,u,v,w,T', ni, nj, 2)
c2 = C.initVars(c2, 'ro', 1.)
c2 = C.initVars(c2, 'u', 1.)
c2 = C.initVars(c2, 'v', 0.)
c2 = C.initVars(c2, 'w', 0.)
c2 = C.initVars(c2, 'Temp', 1.)
a2 = C.addVars([m,c2])

b2 = P.computeVariables2(a2, vars)
test.testA([b2], 2)
