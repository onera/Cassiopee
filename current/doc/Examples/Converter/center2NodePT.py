# - center2Node (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

# center2Node: create a new zone
a = G.cart((0,0,0), (1,1,1), (30,40,1))
C._initVars(a, 'centers:Density', 1.)
b = C.center2Node(a); b[0] = a[0]+'_nodes'
C.convertPyTree2File(b, 'out0.cgns')

# center2Node: modify a variable
a = G.cart((0,0,0), (1,1,1), (30,40,1))
C._initVars(a, 'centers:Density', 1.)
a = C.center2Node(a, 'centers:Density')

# center2Node: modify a container
a = C.center2Node(a, 'FlowSolution#Centers')

C.convertPyTree2File(a, 'out.cgns')
