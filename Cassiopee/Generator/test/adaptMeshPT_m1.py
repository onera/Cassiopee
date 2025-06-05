# - adaptMesh (pyTree) -
import Generator.PyTree as G
import Generator.Mpi as Gmpi
import Converter.PyTree as C
import Converter.Filter2 as Filter2
import Converter.Internal as Internal
import Converter.Mpi as Cmpi
import KCore.test as test

LOCAL = test.getLocal()

a = G.cartNGon((0,0,0),(0.1,0.1,0.1),(11,11,2))
Internal._adaptNGon32NGon4(a)
C._fillEmptyBCWith(a, 'nref', 'BCFarfield')
C._initVars(a, 'F', 1.)
C._initVars(a, '{centers:indicator}=({centers:CoordinateX})>0.35')
if Cmpi.rank == 0: C.convertPyTree2File(a, LOCAL+'/tmp.cgns')
Cmpi.barrier()

a, res = Filter2.loadAndSplit(LOCAL+'/tmp.cgns')
splitInfos={}
splitInfos["graph"]=res[1]
splitInfos["cellGlobalIndex"]=res[5]
splitInfos["faceGlobalIndex"]=res[6]

a2 = Gmpi.adaptMesh(a, indicator="indicator", hook=None, dim=2,
                    conformize=False, splitInfos=splitInfos)
test.testT(a2, Cmpi.rank)

a2 = Gmpi.adaptMesh(a, indicator="indicator", hook=None, dim=2,
                    conformize=True, splitInfos=splitInfos)
test.testT(a2,Cmpi.rank+Cmpi.size)

a = G.cartNGon((0,0,0),(0.1,0.1,0.1),(11,11,11))
Internal._adaptNGon32NGon4(a)
C._fillEmptyBCWith(a, 'nref', 'BCFarfield')
C._initVars(a, "F", 1.)
C._initVars(a, '{centers:indicator}=({centers:CoordinateX})>0.35')
if Cmpi.rank == 0: C.convertPyTree2File(a, LOCAL+"/tmp.cgns")
Cmpi.barrier()

a, res = Filter2.loadAndSplit(LOCAL+'/tmp.cgns')
splitInfos = {}
splitInfos["graph"] = res[1]
splitInfos["cellGlobalIndex"] = res[5]
splitInfos["faceGlobalIndex"] = res[6]

a2 = Gmpi.adaptMesh(a, indicator="indicator", hook=None, dim=3,
                    conformize=False, splitInfos=splitInfos)
test.testT(a2, Cmpi.rank+2*Cmpi.size)

a2 = Gmpi.adaptMesh(a, indicator="indicator", hook=None, dim=3,
                    conformize=True, splitInfos=splitInfos)
test.testT(a2, Cmpi.rank+3*Cmpi.size)
#Cmpi.convertPyTree2File(a2, "out.cgns")
#C.convertPyTree2File(a2, "out_%d.cgns"%(Cmpi.rank))
