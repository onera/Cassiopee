# - newReferenceState (pyTree) -
import Converter.Internal as Internal

# Create a ReferenceState node
n = Internal.newReferenceState(name='ReferenceState'); Internal.printTree(n)
#>> ['ReferenceState',None,[0 son],'ReferenceState_t']

# Create a ReferenceState node and attach it to base
t = Internal.newCGNSTree()
b = Internal.newCGNSBase('Base', 3, 3, parent=t)
n = Internal.newReferenceState(name='ReferenceState', parent=b); Internal.printTree(n)
