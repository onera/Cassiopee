# - rmOrphans (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import KCore.test as test

# -- BE/ME
# Create a 3D ME
a = G.cartPyra((0.4,0.4,0.), (0.1,0.1,0.1), (5,5,5))
b = G.cartPenta((0.,0.,0.), (0.1,0.1,0.1), (5,5,5))
c = G.cartHexa((0.4,0.,0.), (0.1,0.1,0.1), (5,5,5))
a = C.mergeConnectivity([a, b, c], None)
C._initVars(a, '{F}={CoordinateX}')

# Add BCs
sz1 = G.cartHexa((0.4,0.4,0.), (0.1,0.1,0.), (5,5,1))
sz2 = G.cartHexa((0.4,0.8,0.), (0.1,0.,0.1), (5,1,5))
sz3 = G.cartHexa((0.4,0.,0.), (0.1,0.,0.1), (5,1,5))
C._addBC2Zone(a, 'wall1', 'BCWall', subzone=sz1)
C._addBC2Zone(a, 'wall2', 'BCWall', subzone=sz2)
#C._addBC2Zone(a, 'wall3', 'BCWall', subzone=sz3)  # TODO fail

# Delete the second Elements_t node and update the ERs
Internal._rmNodesFromName(a, 'GridElements-2')
Internal._updateElementRange(a)

#import numpy as np
#elts = Internal.getElementBoundaryNodes(a)
#for elt in elts:
#    ec = Internal.getNodeFromName(elt, 'ElementConnectivity')[1]
#    print(elt[0], np.min(ec), np.max(ec))

# There are now orphan nodes, delete them
status = {}
indices = []
#bcInfo = C.getBCs(z)
G._rmOrphans(a, indices=indices, status=status)  # 389 vertices before, 289 after
# recoverBCs # TODO revisit this test once BCs can be added by subzone to a ME

#elts = Internal.getElementBoundaryNodes(a)
#for elt in elts:
#    ec = Internal.getNodeFromName(elt, 'ElementConnectivity')[1]
#    print(elt[0], np.min(ec), np.max(ec))

#C.convertPyTree2File(a, "out.cgns")
test.testT(a, 1)
test.testO(indices, 2)
test.testO(status, 3)
