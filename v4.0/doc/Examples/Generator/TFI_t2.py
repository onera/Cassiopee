# - TFI 3D structure -
import Generator as G
import Transform as T
import KCore.test as test

xo = 0.; yo = 0.; zo = 0.
nx = 21; ny = 21; nz = 21
hx = 1./(nx-1); hy = 1./(ny-1); hz = 1./(nz-1)

# grilles z = cste
fzmin = G.cart((xo,yo,zo), (hx,hy,1.), (nx,ny,1))
fzmax = T.translate(fzmin, (0.,0.,1.))

# grilles x = cste
fxmin = G.cart((xo,yo,zo),(1,hy,hz),(1,ny,nz))
fxmin = T.reorder(fxmin,(3,1,2))
fxmax = T.translate(fxmin, (1.,0.,0.))

# grilles y = cste
fymin = G.cart((xo,yo,zo),(hx,1.,hz),(nx,1,nz))
fymin = T.reorder(fymin,(1,3,2))
fymax = T.translate(fymin, (0.,1.,0.))

r = [fxmin,fxmax,fymin,fymax,fzmin,fzmax]
m = G.TFI(r)
test.testA([m],1)
