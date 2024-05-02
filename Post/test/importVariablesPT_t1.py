# - importVariables (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P
import Converter.Internal as Internal
import KCore.test as test

# z2 sans solutions
z1 = G.cart((0.,0.,0.),(0.1,0.1,0.1),(3,3,3))
C._addBC2Zone(z1,'overlap','BCOverlap','imin')
C._fillEmptyBCWith(z1, 'nref','BCFarfield')
t1 = C.newPyTree(['Base',z1]); t2 = Internal.copyRef(t1)
C._initVars(t1,'centers:cellN',1.)
C._initVars(t1,'Pressure',10.)
C._initVars(t1,'Density',1.)
C._addState(t1[2][1], 'Mach', 0.6)
C._addState(t2[2][1], 'Mach', 0.6)
t2 = P.importVariables(t1,t2)
test.testT(t2,1)

# z2 avec cellN = 0
C._initVars(t2,'centers:cellN',0.)
t2 = P.importVariables(t1,t2)
test.testT(t2,2)

# z1 en centres avec z2 sans solution
Internal._rmNodesByType(t2,'FlowSolution_t')
t1 = C.node2Center(t1)
t2 = P.importVariables(t1,t2)
test.testT(t2,3)

# test t1 en centres/ t2 en noeuds
Internal._rmNodesByType(t2,'FlowSolution_t')
C._initVars(t2,'centers:cellN',0.)
t2 = P.importVariables(t1,t2)
test.testT(t2,4)

# EXTRA GRID : pas de meme nom-memes dimensions
z11 = G.cart((0.,0.,0.),(0.1,0.1,0.1),(3,3,3))
z12 = G.cart((10.,0.,0.),(0.1,0.1,0.1),(3,3,3))
z21 = G.cart((20.,0.,0.),(0.1,0.1,0.1),(3,3,3))
z22 = G.cart((30.,0.,0.),(0.1,0.1,0.1),(3,3,3))
z11[0]=z21[0]; z12[0]='toto'
t1 = C.newPyTree(['Base', [z11,z12]]); t2 = C.newPyTree(['Base', [z21,z22]])
C._initVars(t1,'Pressure',10.)
C._initVars(t1,'centers:cellN',1.)
t2 = P.importVariables(t1, t2)
test.testT(t2,5)

# EXTRA GRID : pas de meme nom-dims en centres de t1 correspondent aux noeuds de t2
z11 = G.cart((0.,0.,0.),(0.1,0.1,0.1),(2,2,2))
z12 = G.cart((10.,0.,0.),(0.1,0.1,0.1),(2,2,2))
z21 = G.cart((20.,0.,0.),(0.1,0.1,0.1),(3,3,3))
z22 = G.cart((30.,0.,0.),(0.1,0.1,0.1),(3,3,3))
z11[0]=z21[0]; z12[0]='toto'
t1 = C.newPyTree(['Base', [z11,z12]]); t2 = C.newPyTree(['Base', [z21,z22]])
C._initVars(t1,'Pressure',10.)
C._initVars(t1,'centers:cellN',1.)
C._initVars(t2,'cellN',0.)
t2 = P.importVariables(t1, t2)
test.testT(t2,6)
