# - newPeriodic (pyTree) -
import Converter.Internal as Internal

# Create a node
n = Internal.newPeriodic(rotationCenter=[0.,0.,0.], rotationAngle=[0.,0.,0.], translation=[0.,0.,0.]); Internal.printTree(n)
#>> ['Periodic',None,[3 sons],'Periodic_t']
#>>    |_['RotationCenter',array(shape=(3,),dtype='float64',order='F'),[0 son],'DataArray_t']
#>>    |_['RotationAngle',array(shape=(3,),dtype='float64',order='F'),[0 son],'DataArray_t']
#>>    |_['Translation',array(shape=(3,),dtype='float64',order='F'),[0 son],'DataArray_t']

# Attach it to a parent node
d = Internal.newGridConnectivityProperty()
Internal.newPeriodic(rotationCenter=[0.,0.,0.], rotationAngle=[0.,0.,0.], translation=[0.,0.,0.], parent=d)
