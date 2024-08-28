# - newElements (pyTree) -
import Converter.Internal as Internal

# Create an elements node
b = Internal.newElements(name='Elements', etype='UserDefined', econnectivity=None, eboundary=0); Internal.printTree(b)
#>> ['Elements',array(shape=(2,),dtype='int32',order='F'),[2 sons],'Element_t']
#>>   |_['ElementConnectivity',None,[0 son],'DataArray_t']
#>>   |_['ElementRange',None,[0 son],'IndexRange_t']

# Attach it to zone
z = Internal.newZone('Zone', [[10],[2],[0]], 'Unstructured')
Internal.newElements(name='Elements', etype='UserDefined', econnectivity=None, eboundary=0, parent=z)
