# - cartr1 (array) -
import Generator as G
import Converter as C

a = G.cartr1((0,0,0), (1.,1.,1.), (1.1,1.2,1.), (10,10,10))
C.convertArrays2File(a, 'out.plt')
