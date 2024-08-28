# - newViscosityModel (pyTree) -
import Converter.Internal as Internal

# Create a node
n = Internal.newViscosityModel(value='Null'); Internal.printTree(n)
#>> ['ViscosityModel',array('Null',dtype='|S1'),[0 son],'ViscosityModel_t']

# Create a node and attach it to parent
t = Internal.newFlowEquationSet()
n = Internal.newViscosityModel(value='SutherlandLaw', parent=t)
