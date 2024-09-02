# - createHook (array) -
import Converter as C
import Generator as G

a = G.cart((0,0,0), (1,1,1), (10,10,10))
hook = C.createHook(a, function='nodes')
