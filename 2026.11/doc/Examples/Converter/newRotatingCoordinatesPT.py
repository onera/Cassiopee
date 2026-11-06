# - newRotatingCoordinates (pyTree) -
import Converter.Internal as Internal

# Create an rotating coordinates node
n = Internal.newRotatingCoordinates(rotationCenter=[0.,0.,0.], rotationRateVector=[0.,0.,0.]); Internal.printTree(n)
#>> ['RotatingCoordinates',None,[2 sons],'RotatingCoordinates_t']
#>>    |_['RotationCenter',array(shape=(3,),dtype='float64',order='F'),[0 son],'DataArray_t']
#>>    |_['RotationRateVector',array(shape=(3,),dtype='float64',order='F'),[0 son],'DataArray_t']

# Attach it to base
b = Internal.newCGNSBase('Base', 3, 3)
Internal.newRotatingCoordinates([0.,0.,0.], [0.,0.,0.], parent=b)
