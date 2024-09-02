# - addChild (pyTree) -
import Converter.Internal as Internal

node = Internal.createNode('myNode', 'DataArray_t', value=1.)
child1 = Internal.createNode('child1', 'DataArray_t', value=2.)
child2 = Internal.createNode('child2', 'DataArray_t', value=3.)

# add children nodes to node
Internal.addChild(node, child1, pos=-1) # at the end
Internal.addChild(node, child2, pos=0) # first
print(node)
#>> ['myNode', array([ 1.]), [['child2', array([ 3.]), [], 'DataArray_t'], ['child1', array([ 2.]), [], 'DataArray_t']], 'DataArray_t']
