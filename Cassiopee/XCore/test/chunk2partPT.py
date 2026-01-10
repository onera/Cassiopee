# - chunk2part (pyTree) -
# mpirun -np 4 python3 chunk2partPT.py

import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal
import Converter.Mpi as Cmpi
import Converter.Filter2 as Filter2
import XCore.xcore
import sys

rank = Cmpi.rank
fileName = 'case.cgns'

# 1 - Make the case
if rank == 0:
    a = G.cartTetra((0,0,0),(1,1,1),(5,7,11))
    a = C.convertArray2NGon(a)
    a = C.initVars(a, '{centers:Density} = {centers:CoordinateX} + sin({centers:CoordinateY}) + cos({centers:CoordinateZ})')
    a = C.initVars(a, '{centers:Pressure} = {centers:CoordinateX} + cos({centers:CoordinateY}) + sin({centers:CoordinateZ})')
    a = C.initVars(a, '{Density} = {CoordinateX} + sin({CoordinateY}) + cos({CoordinateZ})')
    Internal._adaptNGon32NGon4(a)
    C.convertPyTree2File(a, fileName)
Cmpi.barrier()

# 2 - Load
distTree = Filter2.loadAsChunks(fileName)
Cmpi._setProc(distTree, Cmpi.rank)

# 3 - chunk2part
arrays = []
for z in Internal.getZones(distTree):
    cx = Internal.getNodeFromName2(z, 'CoordinateX')[1]
    cy = Internal.getNodeFromName2(z, 'CoordinateY')[1]
    cz = Internal.getNodeFromName2(z, 'CoordinateZ')[1]

    ngon = Internal.getNGonNode(z)
    ngonc = Internal.getNodeFromName1(ngon, 'ElementConnectivity')[1]
    ngonso = Internal.getNodeFromName1(ngon, 'ElementStartOffset')[1]

    nface = Internal.getNFaceNode(z)
    nfacec = Internal.getNodeFromName1(nface, 'ElementConnectivity')[1]
    nfaceso = Internal.getNodeFromName1(nface, 'ElementStartOffset')[1]

    fsolc = Internal.getNodeFromName2(z, 'FlowSolution#Centers')
    denc = Internal.getNodeFromName1(fsolc, 'Density')[1]
    prec = Internal.getNodeFromName1(fsolc, 'Pressure')[1]

    fsol = Internal.getNodeFromName2(z, 'FlowSolution')
    den = Internal.getNodeFromName1(fsol, 'Density')[1]

    arrays.append([cx,cy,cz,ngonc,ngonso,nfacec,nfaceso,[],[]])

#comm_data = list of [neighbor proc (int), interproc faces (array),
#                     corresponding global neighbor ids (array)]
RES = XCore.xcore.chunk2partNGon(arrays)

cells = RES[0]
comm_data = RES[1]
mesh = RES[2]
solc = RES[3]
sol = RES[4]

if rank == 0:
    print(solc)

print('rank', rank, '-> interproc patches:', len(comm_data))

Cmpi.barrier()

name = 'part' + str(Cmpi.rank) + '.cgns'
z1 = Internal.newZone('Zone1')
t1 = C.newPyTree(['Base', z1])
new_zones = Internal.getZones(t1)
C.setFields([mesh], new_zones[0], 'nodes')

C.convertPyTree2File(t1, name)
