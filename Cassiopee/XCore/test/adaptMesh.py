# - adaptMesh (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import XCore.PyTree as X
import Converter.Mpi as Cmpi
import Converter.Internal as I
import numpy as np

a1 = 0.37**2
a2 = 0.45**2
xc = 0.5
yc = 0.5

def Func(x, y, z):
    x2 = (x-xc)**2
    y2 = (y-yc)**2
    return x2+y2 > a1 and x2+y2 < a2

if Cmpi.rank == 0:
    a = G.cartHexa((0,0,0),(0.1,0.1,0.1),(11,11,2))
    a = C.convertArray2NGon(a)
    a = C.fillEmptyBCWith(a, 'wall', 'BCWall', dim=3)
    I._adaptNGon32NGon4(a)
    C.convertPyTree2File(a, "case.cgns")

Cmpi.barrier()

t, res = X.loadAndSplitNGon("case.cgns")
gcells = res[5]
gfaces = res[6]
comm = res[1]
normal2D = np.array([0.0,0.0,1.0])
normal2D = None

AM = X.AdaptMesh_Init(t, normal2D, comm, gcells, gfaces)

itermax = 7 # 7

for iter in range(itermax):
    if Cmpi.rank == 0:
        print("\niter:", iter, flush=True)

    C._initVars(t, 'centers:F', Func, ['centers:CoordinateX', 'centers:CoordinateY', 'centers:CoordinateZ'])
    f = I.getNodeFromName(t, 'F')[1]
    REF = f.astype(dtype=I.E_NpyInt)

    X.AdaptMesh_AssignRefData(AM, REF)

    '''
    t = X.AdaptMesh_ExtractMesh(AM, conformize=1)
    I._adaptNGon42NGon3(t)
    Cmpi.convertPyTree2File(t, "unbalanced%d.cgns"%iter)
    '''

    X.AdaptMesh_LoadBalance(AM)

    '''
    t = X.AdaptMesh_ExtractMesh(AM, conformize=1)
    I._adaptNGon42NGon3(t)
    Cmpi.convertPyTree2File(t, "balanced%d.cgns"%iter)
    '''

    X.AdaptMesh_Adapt(AM)

    t = X.AdaptMesh_ExtractMesh(AM, conformize=1)
    I._adaptNGon42NGon3(t)
    Cmpi.convertPyTree2File(t, "refined%d.cgns"%iter)

X.AdaptMesh_Exit(AM)
