# - dist2WallsEikonal (pyTree) -
import Converter.PyTree as C
import Connector.PyTree as X
import Converter.Internal as Internal
import Dist2Walls.PyTree as DTW
import Geom.PyTree as D
import Generator.PyTree as G
import numpy
import KCore.test as test
DEPTH = 2
snear = 0.4; vmin = 21

# Init wall
body = D.circle((0,0,0),1.,N=60)
res = G.octree([body],[snear], dfar=5., balancing=1)
res = G.octree2Struct(res, vmin=vmin, ext=DEPTH+1,merged=1)
t = C.newPyTree(['Base', res])

# Mark solid and fluid points
X._applyBCOverlaps(t,depth=DEPTH,loc='nodes')
tc = Internal.copyRef(t)
tc = X.setInterpData(t,tc,loc='nodes',storage="inverse")
C._initVars(t,"cellN",1.)
t = X.blankCells(t, [[body]], numpy.array([[1]]), blankingType='node_in')
X._setHoleInterpolatedPoints(t,depth=1,loc='nodes')
C._initVars(t,'{flag}=({cellN}>1.)')
t = DTW.distance2WallsEikonal(t,body,tc=tc,DEPTH=DEPTH,nitmax=10)
test.testT(t,1)
# aux centres
t = C.newPyTree(['Base']); t[2][1][2] = res
Internal._rmNodesFromName(t,"FlowSolution")
# Mark solid and fluid points
X._applyBCOverlaps(t,depth=DEPTH,loc='centers')
tc = C.node2Center(t)
tc = X.setInterpData(t,tc,loc='centers',storage="inverse")
C._initVars(t,"centers:cellN",1.)
t = X.blankCells(t, [[body]], numpy.array([[1]]), blankingType='center_in')
X._setHoleInterpolatedPoints(t,depth=1)
C._initVars(t,'{centers:flag}=({centers:cellN}>1.)')
t = DTW.distance2WallsEikonal(t,body,tc=tc,DEPTH=DEPTH,loc='centers',nitmax=10)
test.testT(t,2)
