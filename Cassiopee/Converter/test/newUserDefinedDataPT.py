# - newUserDefinedData (pyTree) -
import Converter.Internal as Internal

# Create a UserDefinedData node
n = Internal.newUserDefinedData(name='UserDefined', value=None); Internal.printTree(n)
#>> ['UserDefined',None,[0 son],'UserDefinedData_t']

# Create a UserDefinedData node and attach it to a Family node (just an example)
f = Internal.newFamily(name='FamInjection')
n = Internal.newFamilyBC(value='BCInflow', parent=f)
n = Internal.newUserDefinedData(name='.Solver#BC', value=None, parent=n)
