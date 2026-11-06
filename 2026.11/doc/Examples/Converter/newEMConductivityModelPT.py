# - newEMConductivityModel (pyTree) -
import Converter.Internal as Internal

# Create a zone node
n = Internal.newEMConductivityModel(value='Null'); Internal.printTree(n)
#>> ['EMConductivityModel',array('Null',dtype='|S1'),[0 son],'EMConductivityModel_t']

# Create a node and attach it to parent
t = Internal.newFlowEquationSet()
z = Internal.newEMConductivityModel(value='Chemistry_LinRessler', parent=t)
