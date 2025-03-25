# - copyNode (pyTree) -
import Converter.Internal as Internal
import KCore.test as test

a = Internal.newDataArray(value=[1,2,3])
# Copy only the numpy of this node
b = Internal.copyNode(a)
# Modify numpy of b
b[1][0]=5
test.testO([a,b],1)
