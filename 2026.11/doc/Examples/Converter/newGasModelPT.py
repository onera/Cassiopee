# - newGasModel (pyTree) -
import Converter.Internal as Internal

# Create a node
z = Internal.newGasModel(value='Ideal'); Internal.printTree(z)
#>> ['GasModel',array('Ideal',dtype='|S1'),[0 son],'GasModel_t']

# Create a node and attach it to parent
t = Internal.newFlowEquationSet()
n = Internal.newGasModel(value='Ideal', parent=t)
