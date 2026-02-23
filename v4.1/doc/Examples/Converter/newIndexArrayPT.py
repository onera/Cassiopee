# - newIndexArray (pyTree) -
import Converter.Internal as Internal

# Create a node
n = Internal.newIndexArray(name='Index', value=[101,51,22,1024,78]); Internal.printTree(n)
#>> ['Index',array(shape=(5,),dtype='int32',order='F'),[0 son],'IndexArray_t']

# Attach it to a parent node
d = Internal.newFlowSolution('FlowSolution', 'Vertex')
Internal.newIndexArray('Index', [101,51,22,1024,78], parent=d)
