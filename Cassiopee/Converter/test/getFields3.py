# - getFields (pyTree) -
# api=3
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal

# STRUCT
#z = G.cart((0.,0.,0.),(1.,1.,1.),(10,10,10))
#res = C.getFields('GridCoordinates', z, api=3)
#z = C.initVars(z, 'F', 1)
#z = C.initVars(z, 'centers:G', 2)
#res = C.getFields('GridCoordinates', z, api=3)
#res = C.getFields('FlowSolution', z, api=3)
#res = C.getFields('FlowSolution#Centers', z, api=3)

# BE
#z = G.cartTetra((0.,0.,0.),(1.,1.,1.),(10,10,10))
#z = C.initVars(z, 'F', 1)
#z = C.initVars(z, 'centers:G', 2)
#res = C.getFields('GridCoordinates', z, api=3)
#res = C.getFields('FlowSolution', z, api=3)
#res = C.getFields('FlowSolution#Centers', z, api=3)

# ME
z1 = G.cartTetra((0.,0.,0.),(1.,1.,1.),(10,10,10))
z2 = G.cartHexa((9,0.,0.),(1.,1.,1.),(10,10,10))
z = C.mergeConnectivity(z1, z2)
C._addVars(z, 'F')
C._addVars(z, 'centers:G')
Internal.printTree(z)
C._initVars(z, 'F', 1)
##C._initVars(z, 'centers:G', 1) - A TESTER
#res = C.getFields('GridCoordinates', z, api=3)
#res = C.getFields('FlowSolution', z, api=3)
#res = C.getFields('FlowSolution#Centers', z, api=3)
#print(res)
#C.convertPyTree2File(z, 'out.cgns')
#import sys; sys.exit() # END UP


# NGON
#z = G.cartNGon((0.,0.,0.),(1.,1.,1.),(10,10,10), api=3)
#C.convertPyTree2File(z, 'out.cgns')
#print(Internal.getZoneDim(z))
#res = C.getFields('GridCoordinates', z, api=3)
#C._addVars(z, 'F')
#C._addVars(z, 'centers:G')
##C._initVars(z, 'F', 1)
##z = C.initVars(z, 'centers:G', 2)
#res = C.getFields('FlowSolution', z, api=3)
#res = C.getFields('FlowSolution#Centers', z, api=3)
