# - writePyTreeFromFilter (pyTree) -
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.PyTree as C
import Converter.Filter as Filter
import Converter.Mpi as Cmpi
import KCore.test as test

LOCAL = test.getLocal()

# Ecrit la zone en entier
t = C.newPyTree(['Base'])
a = G.cart((0,0,0), (1,1,1), (10,10,10))
t[2][1][2] += [a]
if Cmpi.rank == 0:
    C.convertPyTree2File(t, LOCAL+'/out.hdf')
Cmpi.barrier()

# Prend une subzone et la remplace dans le fichier
t = C.newPyTree(['Base'])
if Cmpi.rank == 0:
    b = T.subzone(a, (2,2,2), (5,5,5)); b[0] = 'cart'
    C._initVars(b, 'CoordinateX', 1.)
    t[2][1][2] += [b]
elif Cmpi.rank == 1:
    b = T.subzone(a, (7,7,7), (9,9,9)); b[0] = 'cart'
    C._initVars(b, 'CoordinateX', 1.)
    t[2][1][2] += [b]

if Cmpi.rank == 0:
    DataSpaceMMRY = [[0,0,0], [1,1,1], [4,4,4], [1,1,1]]
    DataSpaceFILE = [[2,2,2], [1,1,1], [4,4,4], [1,1,1]]
    DataSpaceGLOB = [[10,10,10]]
elif Cmpi.rank == 1:
    DataSpaceMMRY = [[0,0,0], [1,1,1], [3,3,3], [1,1,1]]
    DataSpaceFILE = [[7,7,7], [1,1,1], [3,3,3], [1,1,1]]
    DataSpaceGLOB = [[10,10,10]]

f = {}
f['/Base/cart/GridCoordinates/CoordinateX'] = DataSpaceMMRY+DataSpaceFILE+DataSpaceGLOB
# skelData != None car le noeud existe deja
Filter.writePyTreeFromFilter(t, LOCAL+'/out.hdf', f, skelData=[])
Cmpi.barrier()

if Cmpi.rank == 0:
    r = C.convertFile2PyTree(LOCAL+'/out.hdf')
    test.testT(r, 1)
