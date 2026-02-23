# - newRigidGridMotionType (pyTree) -
import Converter.Internal as Internal

# Create a node
n = Internal.newRigidGridMotionType(value='ConstantRate'); Internal.printTree(n)
#>> ['RigidGridMotionType',array('ConstantRate',dtype='|S1'),[0 son],'RigidGridMotionType_t']

# Attach it to a parent node
z = Internal.newRigidGridMotion(name='Motion', origin=[0.,0.,0.], mtype='Null')
Internal.newRigidGridMotionType(value='ConstantRate', parent=z)
