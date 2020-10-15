# - writePyTreeFromPaths (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Filter as Filter
import KCore.test as test

LOCAL = test.getLocal()

t = C.newPyTree(['Base'])
C.convertPyTree2File(t, LOCAL+'/out.hdf')

a = G.cart((0,0,0), (1,1,1), (10,10,10))
t[2][1][2] += [a]
Filter.writePyTreeFromPaths(t, LOCAL+'/out.hdf', 'CGNSTree/Base/cart')

a = C.convertFile2PyTree(LOCAL+'/out.hdf')
test.testT(a, 1)
