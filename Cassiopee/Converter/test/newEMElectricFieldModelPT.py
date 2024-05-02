# - newEMElectricFieldModel (pyTree) -
import Converter.Internal as Internal

# Create a node
n = Internal.newEMElectricFieldModel(value='Null'); Internal.printTree(n)
#>> ['EMElectricFieldModel',array('Null',dtype='|S1'),[0 son],'EMElectricFieldModel_t']

# Create a node and attach it to parent
t = Internal.newFlowEquationSet()
n = Internal.newEMElectricFieldModel(value='Voltage', parent=t)
