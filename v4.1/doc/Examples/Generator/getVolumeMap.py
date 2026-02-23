# - getVolumeMap (array) -
import Generator as G
import Converter as C

a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,3))
vol = G.getVolumeMap(a)
vol = C.center2Node(vol); vol = C.addVars([a,  vol])
C.convertArrays2File(vol, "out.plt")
