# - addChild (pyTree) -
import Converter.Internal as Internal
import Generator.PyTree as G
import KCore.test as test

node = Internal.createBaseNode('Base', 3)
child1 = G.cart((0.,0.,0.),(1,1,1),(10,10,10))
child2 = G.cart((0.,0.,0.),(1,1,1),(10,10,10))
child3 = G.cart((0.,0.,0.),(1,1,1),(10,10,10))

# add child nodes
Internal.addChild(node, child1, pos=-1) # at the end
Internal.addChild(node, child2, pos=0) # first
test.testT(node, 1)

# add children nodes
node = Internal.createBaseNode('Base', 3)
Internal.addChild(node, child1) # at the end
Internal.addChild(node, [child2,child3], pos=-1) # at the end
test.testT(node, 2)

node = Internal.createBaseNode('Base', 3)
Internal.addChild(node, child1) # at the end
Internal.addChild(node, [child2,child3], pos=0) # first
test.testT(node, 3)
