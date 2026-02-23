# - isStdNode (pyTree) -
import Converter.Internal as Internal
import numpy

# This is a standard node
a = ['toto', numpy.zeros(12), [], 'DataArray_t']
print(Internal.isStdNode(a))
#>> -1

# This is not a standard node
b = ['toto', 'tata']
print(Internal.isStdNode(b))
#>> -2

# This is a list of standard nodes
c = ['titi', numpy.zeros(13), [], 'DataArray_t']
print(Internal.isStdNode([a,c]))
#>> 0
