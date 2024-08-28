# - newRigidGridMotion (pyTree) -
import Converter.Internal as Internal

# Create a node
n = Internal.newRigidGridMotion(name='Motion', origin=[0.,0.,0.], mtype='ConstantRate'); Internal.printTree(n)
#>> ['Motion',None,[2 sons],'RigidGridMotion_t']
#>>    |_['OriginLocation',array(shape=(3,),dtype='float64',order='F'),[0 son],'DataArray_t']
#>>    |_['RigidGridMotionType',array('ConstantRate',dtype='|S1'),[0 son],'RigidGridMotionType_t']

# Attach it to a parent node
z = Internal.newZone('Zone', zsize=[[10],[2],[0]], ztype='Structured')
Internal.newRigidGridMotion(name='Motion', origin=[0.,0.,0.], mtype='ConstantRate', parent=z)
