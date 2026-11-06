# - cartr2 (array) -
import Generator as G
import Converter as C

a = G.cartr2((10,5,1), (1,1,1), (1.5,1.3,1.), (200.,100.,100.))
C.convertArrays2File(a, 'out.plt')
