# - readZones (pyTree) -
import Converter.PyTree as C
import Converter.Distributed as Distributed
import Distributor2.PyTree as Distributor2
import Generator.PyTree as G
import KCore.test as test

LOCAL = test.getLocal()

# Cree le fichier test HDF
a = G.cart((0,0,0), (1,1,1), (10,10,10))
b = G.cart((12,0,0), (1,1,1), (10,10,10))
t = C.newPyTree(['Base',a,b])
C.convertPyTree2File(t, LOCAL+'/test.hdf')

# Relit des zones par procs
t = Distributed.convertFile2SkeletonTree(LOCAL+'/test.hdf')
(t, dic) = Distributor2.distribute(t, NProc=2, algorithm='fast')
t = Distributed.readZones(t, LOCAL+'/test.hdf', rank=0)
t = Distributed.readZones(t, LOCAL+'/test.hdf', rank=1)
test.testT(t, 1)

# Cree le fichier test ADF
a = G.cart((0,0,0), (1,1,1), (10,10,10))
b = G.cart((12,0,0), (1,1,1), (10,10,10))
t = C.newPyTree(['Base',a,b])
C.convertPyTree2File(t, LOCAL+'/test.adf')

# Relit des zones par procs
t = Distributed.convertFile2SkeletonTree(LOCAL+'/test.adf')
(t, dic) = Distributor2.distribute(t, NProc=2, algorithm='fast')
t = Distributed.readZones(t, LOCAL+'/test.adf', rank=0)
t = Distributed.readZones(t, LOCAL+'/test.adf', rank=1)
test.testT(t, 2)
