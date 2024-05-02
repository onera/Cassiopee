# - addGC + BCDataSet (pyTree) -
import Converter.Internal as Internal
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test
N = 10; dim = 2
a = G.cart((0,0,0), (1,1,1), (N,N,1))
a = C.addBC2Zone(a, 'wall', 'BCWall', 'imin')
a = C.addBC2Zone(a, 'wall', 'BCWall', 'jmax')
nobc = 1
for b in Internal.getNodesFromName(a, 'wall*'):
    d = Internal.newBCDataSet(name='BCDataSet', value='UserDefined',
                              gridLocation='FaceCenter', parent=b)
    d = Internal.newBCData('BCNeumann', parent=d)
    if dim == 2: nfaces = (N-1)
    else: nfaces = (N-1)**2
    d1 = Internal.newDataArray('Density', value=nfaces*[1.*nobc], parent=d)
    d2 = Internal.newDataArray('G', value=nfaces*[-1.*nobc], parent=d)
    nobc += 1
C._initVars(a,'centers:Density',2.)
C._initVars(a,'centers:G', 2.)
t = C.newPyTree(['Base',a,dim])
t = Internal.addGhostCells(t,t,2,adaptBCs=0)
test.testT(t,1)
#
N = 10; dim = 3
a = G.cart((0,0,0), (1,1,1), (N,N,N))
C._addBC2Zone(a, 'wall', 'BCWall', 'imin')
C._addBC2Zone(a, 'wall', 'BCWall', 'jmax')
nobc = 1
for b in Internal.getNodesFromName(a, 'wall*'):
    d = Internal.newBCDataSet(name='BCDataSet', value='UserDefined',
                              gridLocation='FaceCenter', parent=b)
    d = Internal.newBCData('BCNeumann', parent=d)
    if dim == 2: nfaces = (N-1)
    else: nfaces = (N-1)**2
    d1 = Internal.newDataArray('Density', value=nfaces*[1.*nobc], parent=d)
    d2 = Internal.newDataArray('G', value=nfaces*[-1.*nobc], parent=d)
    nobc+=1
C._initVars(a,'centers:Density',2.)
C._initVars(a,'centers:G', 2.)
t = C.newPyTree(['Base',a,dim])
t = Internal.addGhostCells(t,t,2,adaptBCs=0)
test.testT(t,2)
#
# located at nodes: nothing is done yet
#
N = 10; dim = 3
a = G.cart((0,0,0), (1,1,1), (N,N,N))
C._addBC2Zone(a, 'wall', 'BCWall', 'imin')
b = Internal.getNodeFromName(a, 'wall')
d = Internal.newBCDataSet(name='BCDataSet', value='UserDefined',
                          gridLocation='Vertex', parent=b)
d = Internal.newBCData('BCNeumann', parent=d)
d1 = Internal.newDataArray('Density', value=nfaces*[1.*nobc], parent=d)
nobc+=1
C._initVars(a,'{centers:Density}={centers:CoordinateX}')
C._initVars(a,'centers:G', 2.)
t = C.newPyTree(['Base',a,dim])
t = Internal.addGhostCells(t,t,2,adaptBCs=0)
test.testT(t,3)
