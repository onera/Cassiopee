# - createGlobalHook (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cart((0,0,0), (1,1,1), (10,10,10))
b = G.cart((9,0,0), (1,1,1), (10,10,10))
hook, indir = C.createGlobalHook([a,b], function='nodes', indir=1)
