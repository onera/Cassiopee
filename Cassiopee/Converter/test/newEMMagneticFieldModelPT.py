# - newEMMagneticFieldModel (pyTree) -
import Converter.Internal as Internal

# Create a node
n = Internal.newEMMagneticFieldModel(value='Null'); Internal.printTree(n)
#>> ['EMMagneticFieldModel',array('Null',dtype='|S1'),[0 son],'EMMagneticFieldModel_t']

# Create a node and attach it to parent
t = Internal.newFlowEquationSet()
z = Internal.newEMMagneticFieldModel(value='Interpolated', parent=t)
