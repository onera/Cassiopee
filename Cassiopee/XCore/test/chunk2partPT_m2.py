# - chunk2part (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal
import Converter.Mpi as Cmpi
import KCore.test as test
import XCore.PyTree as X

LOCAL = test.getLocal()

rank = Cmpi.rank
fileName = LOCAL+'/case.cgns'

# 1 - Make the case
if rank == 0:
    a = G.cartTetra((0,0,0),(1,1,1),(5,7,11))
    a = C.convertArray2NGon(a)
    a = C.initVars(a, '{centers:Density} = {centers:CoordinateX} + sin({centers:CoordinateY}) + cos({centers:CoordinateZ})')
    a = C.initVars(a, '{centers:Pressure} = {centers:CoordinateX} + cos({centers:CoordinateY}) + sin({centers:CoordinateZ})')
    a = C.initVars(a, '{Density} = {CoordinateX} + sin({CoordinateY}) + cos({CoordinateZ})')
    Internal._adaptNGon32NGon4(a) # NGONv4
    C.convertPyTree2File(a, fileName)
Cmpi.barrier()

# 2 - Load
t, res = X.loadAndSplitNGon(fileName)

#Cmpi.convertPyTree2File(t, 'out.cgns')
if Cmpi.rank == 0: test.testT(t, 1)
