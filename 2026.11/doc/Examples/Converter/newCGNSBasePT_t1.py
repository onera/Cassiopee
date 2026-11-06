# - newCGNSBase (pyTree) -
import Converter.Internal as Internal
import KCore.test as test

# Create a base node
b = Internal.newCGNSBase('Base')
test.testT(b, 1)

# Create a base node and attach it to tree
t = Internal.newCGNSTree()
Internal.newCGNSBase('Base', 3, 3, parent=t)
test.testT(t, 2)
