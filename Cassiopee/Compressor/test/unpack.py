# - unpack -
import Compressor
import Generator.PyTree as G
a = G.cart((0,0,0), (1,1,1), (1000,100,100))
b = Compressor.pack(a)
c = Compressor.unpack(b)
