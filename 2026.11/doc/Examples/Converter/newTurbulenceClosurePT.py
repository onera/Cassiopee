# - newTurbulenceClosure (pyTree) -
import Converter.Internal as Internal

# Create a node
n = Internal.newTurbulenceClosure(value='Null'); Internal.printTree(n)
#>> ['TurbulenceClosure',array('Null',dtype='|S1'),[0 son],'TurbulenceClosure_t']

# Create a node and attach it to parent
t = Internal.newFlowEquationSet()
n = Internal.newTurbulenceClosure(value='ReynoldsStress', parent=t)
