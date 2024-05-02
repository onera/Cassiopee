# - splitMultiplePts (pyTree) -
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.PyTree as C
import KCore.test as test

# Structure 2D
z0 = G.cart((0.,0.,0.),(0.1,0.1,1.),(10,10,1))
z0 = C.initVars(z0, 'F', 1); z0 = C.initVars(z0, 'centers:G', 1)
z1 = T.subzone(z0,(1,1,1),(5,10,1)); z1[0] = 'cart1'
z2 = T.subzone(z0,(5,1,1),(10,5,1)); z2[0] = 'cart2'
z3 = T.subzone(z0,(5,5,1),(10,10,1)); z3[0] = 'cart3'
t = C.newPyTree(['Base',2])
t[2][1][2] += [z1,z2,z3]
t = C.initVars(t, 'F', 1.); t = C.initVars(t, 'centers:G', 2.)
t = T.splitMultiplePts(t,dim=2)
test.testT(t,1)

# Structure 3D
nk = 10
z0 = G.cart((0.,0.,0.),(0.1,0.1,1.),(10,10,nk))
z0 = C.initVars(z0, 'F', 1); z0 = C.initVars(z0, 'centers:G', 1)
z1 = T.subzone(z0,(1,1,1),(5,10,nk)); z1[0] = 'cart1'
z2 = T.subzone(z0,(5,1,1),(10,5,nk)); z2[0] = 'cart2'
z3 = T.subzone(z0,(5,5,1),(10,10,nk)); z3[0] = 'cart3'
z0 = T.translate(z0,(-0.9,0.,0.)); z0[0] = 'cart0'
z4 = G.cart((-0.9,0.9,0.),(0.1,0.1,1.),(19,5,nk)); z4[0] = 'cart4'
t = C.newPyTree(['Base']); t[2][1][2] += [z0,z1,z2,z3,z4]
t = C.initVars(t,'F',1.); t = C.initVars(t,'centers:G',2.)
t2 = T.splitMultiplePts(t,dim=3)
test.testT(t2,2)

# Sur une liste de bases
t2[2][1:] = T.splitMultiplePts(t[2][1:])
test.testT(t2,3)
# Sur une base
t2[2][1] = T.splitMultiplePts(t[2][1])
test.testT(t2,4)
# Sur une liste de zones
t2[2][1][2] = T.splitMultiplePts(t[2][1][2])
test.testT(t2,5)
# Sur une zone
zones = T.splitMultiplePts(t[2][1][2][0])
t2 = C.newPyTree(['SplitZones']); t2[2][1][2] = zones
test.testT(t2,6)
