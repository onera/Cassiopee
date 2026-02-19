# - meshDeviation (pyTree) -
import OCC.PyTree as OCC
import Converter.PyTree as C
import Post.PyTree as P
import Converter.Internal as Internal

# meshDeviation on a mesh circle
hook = OCC.createEmptyCAD()
OCC._addCircle(hook, (0,0,0), (0,0,1), 2.)
m = OCC.meshAll(hook, hmin=0.05, hmax=0.05)
OCC._meshDeviation(hook, m, loc="centers")

# meshDeviation on a mesh sphere
hook = OCC.createEmptyCAD()
OCC._addSphere(hook, (0.,0.,0.), 1.)
m = OCC.meshAll(hook, hmin=0.05, hmax=0.05)

# doesnt work on edges? why?
Internal._rmNodesFromName(m, 'EDGES')
OCC._meshDeviation(hook, m, loc="nodes")

C.convertPyTree2File(m, 'out.cgns')
