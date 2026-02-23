# - newGoverningEquation (pyTree) -
import Converter.Internal as Internal

# Create a node
n = Internal.newGoverningEquations(value='Euler'); Internal.printTree(n)
#>> ['GoverningEquations',array('Euler',dtype='|S1'),[0 son],'GoverningEquations_t']

# Create a node and attach it to parent
t = Internal.newFlowEquationSet()
n = Internal.newGoverningEquations(value='NSTurbulent', parent=t)
