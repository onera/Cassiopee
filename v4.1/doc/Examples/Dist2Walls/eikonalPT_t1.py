# - eikonal (pyTree) -
import Dist2Walls.PyTree as Dist2Walls
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D
import Connector.PyTree as X
import numpy
import KCore.test as test
# Bloc cartesien
a = G.cart((0.,0.,0.),(0.1,0.1,0.1),(128,128,1))
C._initVars(a,'cellN', 1)

# Init wall
sphere = D.sphere((6.4,6.4,0), 1., 30)
sphere = C.convertArray2Tetra(sphere)
sphere = G.close(sphere)
t = C.newPyTree(['Base']); t[2][1][2] = [a]
t = X.blankCellsTri(t, [[sphere]], numpy.array([[1]]), blankingType='node_in')
C._initVars(t,'{Phi}=1.e12*({cellN}>0.)')
C._initVars(t,'{speed}=%f'%(1./0.1))

# Eikonal
t = Dist2Walls.eikonal(t)
C._rmVars(t,['speed'])
test.testT(t,1)

# loc='centers'
a = G.cart((0.,0.,0.),(0.1,0.1,0.1),(128,128,1))
C._initVars(a,'centers:cellN', 1)
t = C.newPyTree(['Base']); t[2][1][2] = [a]
t = X.blankCellsTri(t, [[sphere]], numpy.array([[1]]), blankingType='center_in')
C._initVars(t,'{centers:speed}=%f'%(1./0.1))
# Eikonal
C._initVars(t,'{centers:flag}=({centers:cellN})<1.') # pts sources
C._initVars(t,'{centers:Phi}=1.e12*({centers:cellN}>0.)')
t = Dist2Walls.eikonal(t,loc='centers')
C._rmVars(t,['centers:speed'])
test.testT(t,2)
