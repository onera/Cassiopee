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
    a = G.cartHexa((0,0,0),(0.1,0.1,0.1),(4,3,2))
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

itermax = 7 

for iter in range(itermax):
    if Cmpi.rank == 0:
        print("\niter:", iter, flush=True)
    
    C._initVars(t, 'centers:F', Func, ['centers:CoordinateX', 'centers:CoordinateY', 'centers:CoordinateZ'])
    f = I.getNodeFromName(t, 'F')[1]
    #REF = np.empty(len(f), dtype=I.E_NpyInt)
    #REF[:] = f.astype(int)
    REF = f.astype(dtype=np.int32)
    
    X.AdaptMesh_AssignRefData(AM, REF)

    #X.AdaptMesh_LoadBalance(AM)
    
    X.AdaptMesh_Adapt(AM)

    t = X.AdaptMesh_ExtractMesh(AM, conformize=1)

X.AdaptMesh_Exit(AM)

Cmpi.convertPyTree2File(t, "refined"+str(Cmpi.size)+".cgns")

