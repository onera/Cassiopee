# - newConvergenceHistory (pyTree) -
import Converter.Internal as Internal

# Create a ConvergenceHistory node
n = Internal.newConvergenceHistory(name='ZoneConvergenceHistory', value=100); Internal.printTree(n)
#>> ['ZoneConvergenceHistory',array([100],dtype='int32'),[0 son],'ConvergenceHistory_t']

# Create a ConvergenceHistory node and attach it to base
t = Internal.newCGNSTree()
b = Internal.newCGNSBase('Base', 3, 3, parent=t)
n = Internal.newConvergenceHistory(name='GlobalConvergenceHistory', value=100, parent=b)
