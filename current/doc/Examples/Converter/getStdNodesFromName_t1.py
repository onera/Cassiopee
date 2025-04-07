# - getStandardNodeFromName (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

# zone
a = G.cart((0.,0.,0.),(1.,1.,1.),(10,10,10))
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imax')
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',2.)
n = C.getStdNodesFromName(a, 'GridCoordinates')
#test.testT(n) # ce n'est pas un arbre!

# arbre
b = G.cart( (-1,0,0), (1,1,1), (2,2,1) ); b[0] = 'cart2'
t = C.newPyTree(['Base',2]); t[2][1][2] += [a,b]
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
nodes = []
for z in t[2][1][2]:
    n = C.getStdNodesFromName(z, 'GridCoordinates')
    nodes.append(n)
#test.testO(nodes,2)
