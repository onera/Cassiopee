# test : importVariables (pyTree)
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P
import KCore.test as test

# z2 sans solutions
t1 = C.newPyTree(['Base']); t2 = C.newPyTree(['Base'])
z1 = G.cart((0.,0.,0.),(0.1,0.1,0.1),(3,3,3))
z1 = C.addBC2Zone(z1,'overlap','BCOverlap','imin')
z1 = C.fillEmptyBCWith(z1, 'nref','BCFarfield')
t1[2][1][2].append(z1); t2[2][1][2].append(z1)
t1 = C.initVars(t1,'centers:cellN',1.)
t2 = C.initVars(t1,'centers:Density',1.)
C._initVars(t1,'Pressure',10.)
C._initVars(t1,'Density',1.)
t1[2][1] = C.addState(t1[2][1], 'Mach', 0.6)
t2[2][1] = C.addState(t2[2][1], 'Mach', 0.6)
t2 = P.importVariables(t1,t2)
test.testT(t2,1)

# z2 avec cellN = 0
C._initVars(t2,'centers:cellN',0.)
t2 = P.importVariables(t1,t2)
test.testT(t2,2)

# z1 en centres avec z2 sans solution
t2 = C.newPyTree(['Base',z1])
t1 = C.node2Center(t1)
t2 = P.importVariables(t1,t2)
test.testT(t2,3)

# test t1 en centres/ t2 en noeuds
t2 = C.newPyTree(['Base']); t2[2][1][2].append(z1)
C._initVars(t2,'centers:cellN',0.)
t2 = P.importVariables(t1,t2)
test.testT(t2,4)

