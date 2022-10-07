# - trace (pyTree) -
import Generator.PyTree as G
import Converter.Mpi as Cmpi

Cmpi.trace('>>> Start')
a = G.cart((0,0,0), (1,1,1), (100,100,100))
Cmpi.trace('>>> After cart')
