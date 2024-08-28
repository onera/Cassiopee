# - newGravity (pyTree) -
import Converter.Internal as Internal

# Create a gravity node
n = Internal.newGravity(value=[0.,0.,9.81]); Internal.printTree(n)
#>> ['Gravity',None,[1 son],'Gravity_t']
#>>    |_['GravityVector',array(shape=(3,),dtype='float64',order='F'),[0 son],'DataArray_t']

# Create a Gravity node and attach it to base
t = Internal.newCGNSTree()
b = Internal.newCGNSBase('Base', 3, 3, parent=t)
n = Internal.newGravity(value=[0.,0.,9.81], parent=b)
