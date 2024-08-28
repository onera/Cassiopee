# - newFamily (pyTree) -
import Converter.Internal as Internal

# Create a Family node
n = Internal.newFamily(name='Wing'); Internal.printTree(n)
#>> ['Wing',None,[0 son],'Family_t']

# Create a Family node and attach it to base
t = Internal.newCGNSTree()
b = Internal.newCGNSBase('Base', 3, 3, parent=t)
n = Internal.newFamily(name='Wing', parent=b)
