# - convertFile2SkeletonTree (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Distributed as Distributed
import KCore.test as test

LOCAL = test.getLocal()

a = G.cart((0.,0.,0.),(0.1,0.1,0.1),(11,11,11))
t = C.newPyTree(['Base', a])
C.convertPyTree2File(t, LOCAL+'/in.cgns')
t1 = Distributed.convertFile2SkeletonTree(LOCAL+'/in.cgns')
test.testT(t1, 1)
