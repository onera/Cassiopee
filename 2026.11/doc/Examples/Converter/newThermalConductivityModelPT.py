# - newThermalConductivityModel (pyTree) -
import Converter.Internal as Internal

# Create a node
n = Internal.newThermalConductivityModel(value='Null'); Internal.printTree(n)
#>> ['ThermalConductivityModel',array('Null',dtype='|S1'),[0 son],'ThermalConductivityModel_t']

# Create a node and attach it to parent
t = Internal.newFlowEquationSet()
n = Internal.newThermalConductivityModel(value='SutherlandLaw', parent=t); Internal.printTree(t)
