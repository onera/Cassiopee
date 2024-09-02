# - optimizeOverlap (pyTree) -
import Converter.PyTree as C
import Converter.Mpi as Cmpi
import Generator.PyTree as G
import Transform.PyTree as T
import Connector.PyTree as X
import Connector.Mpi as Xmpi
import Distributor2.PyTree as Distributor2
import KCore.test as test

LOCAL = test.getLocal()

# Cree le fichier test
if Cmpi.rank == 0:
    Ni = 50; Nj = 50; Nk = 2
    a = G.cart((0,0,0),(1./(Ni-1), 1./(Nj-1),1), (Ni,Nj,Nk))
    b = G.cart((0,0,0),(2./(Ni-1), 2./(Nj-1),1), (Ni,Nj,Nk)); b[0] = 'cart2'
    a = T.rotate(a, (0,0,0), (0,0,1), 10.); a = T.translate(a, (0.5,0.5,0))
    t = C.newPyTree(['Base1','Base2'])
    t[2][1][2].append(a); t[2][2][2].append(b)
    t = C.fillEmptyBCWith(t, 'overlap', 'BCOverlap', dim=2)
    t = C.addVars(t,'Density'); t = C.addVars(t,'centers:Pressure')
    C.convertPyTree2File(t, LOCAL+'/in.cgns')
Cmpi.barrier()

t = Cmpi.convertFile2SkeletonTree(LOCAL+'/in.cgns')
(t, dic) = Distributor2.distribute(t, NProc=Cmpi.size, algorithm='fast')
t = Cmpi.readZones(t, LOCAL+'/in.cgns', rank=Cmpi.rank)

t2 = Xmpi.optimizeOverlap(t)
if Cmpi.rank == 0: test.testT(t2, 1)
# priorite sur Base2
#t2 = X.optimizeOverlap(t, priorities=['Base1',1,'Base2',0])
#if Cmpi.rank == 0: test.testT(t2,2)
