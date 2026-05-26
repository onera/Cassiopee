# - rmNodes (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
a = G.cart((0,0,0),(1,1,1),(10,10,10))
a = C.addVars(a, ['Density', 'centers:cellN', 'rou', 'rov', 'Hx', 'Hy'])
b = C.rmNodes(a, 'FlowSolution#Centers')
t = C.newPyTree(['Base']); t[2][1][2].append(b)
C.convertPyTree2File(t, 'out.cgns')
