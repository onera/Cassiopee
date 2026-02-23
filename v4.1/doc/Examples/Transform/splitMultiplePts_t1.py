# - splitMultiplePts (array)
import Generator as G
import Transform as T
import Converter as C
import KCore.test as test

# Structure 2D
z0 = G.cart((0.,0.,0.),(0.1,0.1,1.),(10,10,1))
z0 = C.initVars(z0, 'F',1)
z1 = T.subzone(z0,(1,1,1),(5,10,1))
z2 = T.subzone(z0,(5,1,1),(10,5,1))
z3 = T.subzone(z0,(5,5,1),(10,10,1))
zones = [z1,z2,z3]
zones2 = T.splitMultiplePts(zones,dim=2)
test.testA(zones2, 1)

# Structure 3D
nk = 2
z0 = G.cart((0.,0.,0.),(0.1,0.1,1.),(10,10,nk))
z0 = C.initVars(z0, 'F',1)
z1 = T.subzone(z0,(1,1,1),(5,10,nk))
z2 = T.subzone(z0,(5,1,1),(10,5,nk))
z3 = T.subzone(z0,(5,5,1),(10,10,nk))
z4 = G.cart((-0.9,0.9,0.),(0.1,0.1,1.),(19,5,nk))
C._addVars(z4, 'F')
z0 = T.translate(z0,(-0.9,0.,0.))
zones = [z0,z1,z2,z3,z4]
zones2 = T.splitMultiplePts(zones,dim=3)
test.testA(zones2,2)
