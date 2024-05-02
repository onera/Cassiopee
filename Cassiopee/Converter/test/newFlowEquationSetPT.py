# - newFlowEquationSet (pyTree) -
import Converter.Internal as Internal

# Create a node
n = Internal.newFlowEquationSet(); Internal.printTree(n)
#>> ['FlowEquationSet',None,[0 son],'FlowEquationSet_t']

# Create a node and attach it to parent
t = Internal.newCGNSTree()
b = Internal.newCGNSBase('Base', 3, 3, parent=t)
n = Internal.newFlowEquationSet(parent=b)
