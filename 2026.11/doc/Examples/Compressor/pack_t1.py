# - pack -
import Compressor
import KCore.test as test
import Generator.PyTree as G
a = G.cart((0,0,0), (1,1,1), (10,10,10))
b = Compressor.pack(a)
test.testO(b, 1)
