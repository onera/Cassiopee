# - newTurbulenceModel (pyTree) -
import Converter.Internal as Internal

# Create a node
n = Internal.newTurbulenceModel(value='TwoEquation_MenterSST'); Internal.printTree(n)
#>> ['TurbulenceModel',array('TwoEquation_MenterSST',dtype='|S1'),[0 son],'TurbulenceModel_t']

# Create a node and attach it to parent
t = Internal.newFlowEquationSet()
n = Internal.newTurbulenceModel(value='OneEquation_SpalartAllmaras', parent=t)
