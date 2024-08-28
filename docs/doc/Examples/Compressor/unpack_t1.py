# - unpack -
import Compressor
import Generator.PyTree as G
import KCore.test as test
a = G.cart((0,0,0), (1,1,1), (30,30,20))
b = Compressor.pack(a)
c = Compressor.unpack(b)
test.testT(a, 1)
