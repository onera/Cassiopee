# - TFI 3D (array) -
import Generator as G
import Converter as C
import Transform as T

xo = 0.; yo = 0.; zo = 0.
ni = 21; nj = 35; nk = 40
hi = 1./(ni-1); hj = 1./(nj-1); hk = 1./(nk-1)

# grilles z = cste
fkmin = G.cart((xo,yo,zo), (hi,hj,1.), (ni,nj,1))
fkmax = T.translate(fkmin, (0.,0.,1.))

# grilles x = cste
fimin = G.cart((xo,yo,zo),(1,hj,hk),(1,nj,nk))
fimin = T.reorder(fimin,(3,1,2))
fimax = T.translate(fimin, (1.,0.,0.))

# grilles y = cste
fjmin = G.cart((xo,yo,zo),(hi,1.,hk),(ni,1,nk))
fjmin = T.reorder(fjmin,(1,3,2))
fjmax = T.translate(fjmin, (0.,1.,0.))
r = [fimin, fimax, fjmin, fjmax, fkmin, fkmax]
m = G.TFI(r)
C.convertArrays2File(r+[m],"out.plt","bin_tp")
