# - newBaseIterativeData (pyTree) -
import Converter.Internal as Internal

# Create a zone node
n = Internal.newBaseIterativeData(name='BaseIterativeData', nsteps=100, itype='IterationValues'); Internal.printTree(n)
#>> ['BaseIterativeData',array([100],dtype='int32'),[1 son],'BaseIterativeData_t']
#>>    |_['IterationValues',array(shape=(100,),dtype='int32',order='F'),[0 son],'DataArray_t']

# Create a node and attach it to parent
t = Internal.newCGNSTree()
b = Internal.newCGNSBase('Base', 3, 3, parent=t)
n = Internal.newBaseIterativeData(name='BaseIterativeData', nsteps=100, itype='IterationValues', parent=b)
