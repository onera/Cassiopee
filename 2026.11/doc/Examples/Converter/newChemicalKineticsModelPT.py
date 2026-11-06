# - newChemicalKineticsModel (pyTree) -
import Converter.Internal as Internal

# Create a node
n = Internal.newChemicalKineticsModel(value='Null'); Internal.printTree(n)
#>> ['ChemicalKineticsModel',array('Null',dtype='|S1'),[0 son],'ChemicalKineticsModel_t']

# Create a node and attach it to parent
t = Internal.newFlowEquationSet()
n = Internal.newChemicalKineticsModel(value='ChemicalNonequilib', parent=t)
