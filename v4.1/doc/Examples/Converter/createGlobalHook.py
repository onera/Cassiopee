# - createGlobalHook (array) -
import Converter as C
import Generator as G

a = G.cart((0,0,0), (1,1,1), (10,10,10))
b = G.cart((9,0,0), (1,1,1), (10,10,10))
hook = C.createGlobalHook([a,b], function='nodes')
