# - boolean diffSurf (PyTree) -
import Generator.PyTree as G
import Geom.PyTree as D
import Converter.PyTree as C
import Intersector.PyTree as XOR

# octree
s = D.sphere((0,0,0), 1., 100); snear = 0.1
t = G.octree([s],[snear], dfar=5., balancing=1,ratio=2)

# ngon converion
t = C.convertArray2NGon(t)
# ngon conformization
t = C.conformizeNGon(t); t = G.close(t)
# ngon close cells
t = XOR.closeCells(t)
#t = XOR.reorientExternalFaces(t)

# ngon converion of the sphere
s = C.convertArray2NGon(s)
# ngon converion to the nuga format
s = XOR.convertNGON2DToNGON3D(s)
#s = XOR.reorientExternalFaces(s)

# Boolean operation
x = XOR.diffSurf(t, s, tol=0., preserve_right=1, agg_mode=2) # agg_mode=2 : full mode aggregation

x = XOR.agglomerateSmallCells(x, vmin=0., vratio=0.1)

t = C.newPyTree(['Base',2]); t[2][1][2].append(x)
C.convertPyTree2File(t, 'diffs.cgns')
