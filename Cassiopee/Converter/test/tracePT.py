# - trace (pyTree) -
import Generator.PyTree as G
import Converter.Mpi as Cmpi

Cmpi.trace('>>> Start', fileName="stdout")
a = G.cart((0,0,0), (1,1,1), (100,100,100))
Cmpi.trace('>>> After cart', fileName="stdout")
#>> 0: >>> Start [0 secs][8.5 MB]
#>> 0: >>> After cart [0.0060176 secs][9.2 MB]