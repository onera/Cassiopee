# - connectMatch 3D MPI (pyTree) -
import Generator.PyTree    as G
import Converter.PyTree    as C
import Converter.Mpi       as Cmpi
import Connector.Mpi       as Xmpi
import Converter.Filter    as Filter
import KCore.test          as test

LOCAL = test.getLocal()

# Cree le fichier test
if Cmpi.rank == 0:
    a = G.cart((0.,0.,0.), (0.1, 0.1, 0.1), (11, 21, 2))

    # partiellement coincident
    # --- CL
    a2 = G.cart((1., 0.5, 0.), (0.1, 0.1, 0.1), (11, 21, 2))
    a2 = C.addBC2Zone(a2, 'wall1', 'BCWall', 'jmax')
    a2 = C.addBC2Zone(a2, 'overlap1', 'BCOverlap', 'jmin')
    # --- champ aux centres
    C._initVars(a2, 'centers:Density', 1.)
    # --- champ aux noeuds
    C._initVars(a2, 'cellN', 2.)
    t = C.newPyTree(['Base',a,a2])
    # --- Equation state
    t[2][1] = C.addState(t[2][1], 'EquationDimension', 3)
    C.convertPyTree2File(t, LOCAL+'/in.cgns')
Cmpi.barrier()

h  = Filter.Handle(LOCAL+'/in.cgns')
a  = h.loadAndDistribute()
a  = Xmpi.connectMatch(a)

if Cmpi.rank == 0:
    test.testT(t,1)
