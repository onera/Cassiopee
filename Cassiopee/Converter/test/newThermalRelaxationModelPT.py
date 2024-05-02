# - newThermalRelaxationModel (pyTree) -
import Converter.Internal as Internal

# Create a node
n = Internal.newThermalRelaxationModel(value='Null'); Internal.printTree(n)
#>> ['ThermalRelaxationModel',array('Null',dtype='|S1'),[0 son],'ThermalRelaxationModel_t']

# Create a node and attach it to parent
t = Internal.newFlowEquationSet()
n = Internal.newThermalRelaxationModel(value='ThermalNonequilib', parent=t)
