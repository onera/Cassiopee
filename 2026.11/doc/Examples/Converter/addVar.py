# - addVars (array) -
import Converter as C
import Generator as G

a = G.cart((0,0,0), (1,1,1), (10,10,11))

# Add a variable defined by a string
a = C.addVars(a, 'ro')
a = C.addVars(a, 'cellN')
C.convertArrays2File(a, 'out1.plt')

# Add variables defined by a list of varNames
a = C.addVars(a, ['rou','rov'])
C.convertArrays2File(a, 'out2.plt')
