# - setInterpData2 (pyTree) -
import Generator.PyTree as G
import Connector.Mpi as Xmpi
import Converter.PyTree as C
import Converter.Filter as Filter
import Converter.Mpi as Cmpi
import KCore.test as test

LOCAL = test.getLocal()

if Cmpi.rank==0:
    aD = G.cart((0,0,0),(1,1,1), (11,11,11))
    Cmpi._setProc(aD,0)
    aR = G.cart((0,0,0),(0.5,0.5,0.5), (21,21,21))
    Cmpi._setProc(aR,1)
    t = C.newPyTree(['Base']); t[2][1][2] = [aR,aD]
    C.convertPyTree2File(aR, LOCAL+'/rcv.cgns')
    C.convertPyTree2File(aD, LOCAL+'/dnr.cgns')

Cmpi.barrier()

hR = Filter.Handle(LOCAL+'/rcv.cgns')
hD = Filter.Handle(LOCAL+'/dnr.cgns')
aR  = hR.loadFromProc()
aD = hD.loadFromProc()

Xmpi._setInterpData2(aR, aD, loc='centers', cartesian=False)
if Cmpi.rank==0: test.testT(aD,1)
