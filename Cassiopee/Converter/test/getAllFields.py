# - getFields (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

z = G.cart((0.,0.,0.),(1.,1.,1.),(10,10,10))
print(C.getAllFields(z, 'nodes'))

z = C.addVars(z, 'centers:cellN')
print(C.getAllFields(z, 'centers'))
