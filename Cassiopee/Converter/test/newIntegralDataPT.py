# - newIntegralData (pyTree) -
import Converter.Internal as Internal

# Create an integral data node
n = Internal.newIntegralData(name='IntegralData'); Internal.printTree(n)
#>> ['IntegralData',None,[0 son],'IntegralData_t']

# Attach it to base
b = Internal.newCGNSBase('Base', 3, 3)
Internal.newIntegralData('IntegralData', parent=b)
