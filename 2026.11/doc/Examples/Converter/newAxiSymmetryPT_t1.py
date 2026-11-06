# - newAxiSymmetry (pyTree) -
import Converter.Internal as Internal
import KCore.test as test

# Create an axisymmetry node
n = Internal.newAxiSymmetry(referencePoint=[0.,0.,0.], axisVector=[0.,0.,0.])
test.testT(n, 1)

# Attach it to base
b = Internal.newCGNSBase('Base', 3, 3)
Internal.newAxiSymmetry([0.,0.,0.], [0.,0.,0.], parent=b)
test.testT(b, 2)
