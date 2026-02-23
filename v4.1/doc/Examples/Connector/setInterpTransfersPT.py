# - setInterpTransfers (pyTree) -
import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G

a = G.cylinder( (0,0,0), 1, 2, 0, 360, 1, (60, 20, 3) )
b = G.cylinder( (0,0,0), 1, 2, 3, 160, 1, (30, 10, 3) )
a = C.addBC2Zone(a, 'wall', 'BCWall', 'jmin')
a = C.addBC2Zone(a, 'match', 'BCMatch', 'imin', a, 'imax', trirac=[1,2,3])
a = C.addBC2Zone(a, 'match', 'BCMatch', 'imax', a, 'imin', trirac=[1,2,3])
b = C.addBC2Zone(b, 'wall', 'BCWall', 'jmin')
b = C.addBC2Zone(b, 'overlap', 'BCOverlap', 'imin')
b = C.addBC2Zone(b, 'overlap', 'BCOverlap', 'imax')
t1 = C.newPyTree(['Base']); t1[2][1][2] = [a];
t2 = C.newPyTree(['Base']); t2[2][1][2] = [b]
t1 = C.fillEmptyBCWith(t1, 'nref', 'BCFarfield')
t2 = C.fillEmptyBCWith(t2, 'nref', 'BCFarfield')

t1 = C.initVars(t1, '{Density}= 1.')
t2 = C.initVars(t2, '{Density}=-1.')
t2 = C.node2Center(t2, ['Density'])

t2 = X.applyBCOverlaps(t2, depth=1)
t1 = X.setInterpData(t2, t1, double_wall=1, loc='centers',
                     storage='inverse', order=3)
t2 = X.setInterpTransfers(t2, t1, variables=['Density'])
C.convertPyTree2File(t2, 'out.cgns')
