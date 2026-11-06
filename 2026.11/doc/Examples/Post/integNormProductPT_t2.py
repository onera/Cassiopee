# - integNormProduct (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P
import KCore.test as test

N = 20

#-------------------------------
# 2D
#-------------------------------

# STRUCT
z = G.cart((0,0,0), (1./(N-1),1./(N-1),1), (N,N,1))
z = C.initVars(z, 'VeloxityX', 2.)
z = C.initVars(z, 'VeloxityY', 2.)
z = C.initVars(z, 'VeloxityZ', 2.)
z = C.initVars(z, 'centers:MomentumX', 4.)
z = C.initVars(z, 'centers:MomentumY', 4.)
z = C.initVars(z, 'centers:MomentumZ', 4.)
res = [P.integNormProduct(z, ['VeloxityX','VeloxityY','VeloxityZ']), P.integNormProduct(z, ['centers:MomentumX','centers:MomentumY','centers:MomentumZ'])]
test.testO(res,1)

# BE (TRI)
z = G.cartTetra((0,0,0), (1./(N-1),1./(N-1),1), (N,N,1))
z = C.initVars(z, 'VeloxityX', 2.)
z = C.initVars(z, 'VeloxityY', 2.)
z = C.initVars(z, 'VeloxityZ', 2.)
z = C.initVars(z, 'centers:MomentumX', 4.)
z = C.initVars(z, 'centers:MomentumY', 4.)
z = C.initVars(z, 'centers:MomentumZ', 4.)
res = [P.integNormProduct(z, ['VeloxityX','VeloxityY','VeloxityZ']), P.integNormProduct(z, ['centers:MomentumX','centers:MomentumY','centers:MomentumZ'])]
test.testO(res,2)

# BE (QUAD)
z = G.cartHexa((0,0,0), (1./(N-1),1./(N-1),1), (N,N,1))
z = C.initVars(z, 'VeloxityX', 2.)
z = C.initVars(z, 'VeloxityY', 2.)
z = C.initVars(z, 'VeloxityZ', 2.)
z = C.initVars(z, 'centers:MomentumX', 4.)
z = C.initVars(z, 'centers:MomentumY', 4.)
z = C.initVars(z, 'centers:MomentumZ', 4.)
res = [P.integNormProduct(z, ['VeloxityX','VeloxityY','VeloxityZ']), P.integNormProduct(z, ['centers:MomentumX','centers:MomentumY','centers:MomentumZ'])]
test.testO(res,3)

# ME (QUAD+TRI)
a = G.cartHexa((0.,0.,0.), (1./(N-1),1./(N-1),1), (N,N,1))
b = G.cartTetra((1.,0.,0.), (1./(N-1),1./(N-1),1), (N,N,1))
z = C.mergeConnectivity(a, b, boundary=0)
z = C.initVars(z, 'VeloxityX', 2.)
z = C.initVars(z, 'VeloxityY', 2.)
z = C.initVars(z, 'VeloxityZ', 2.)
z = C.initVars(z, 'centers:MomentumX', 4.)
z = C.initVars(z, 'centers:MomentumY', 4.)
z = C.initVars(z, 'centers:MomentumZ', 4.)
res = [P.integNormProduct(z, ['VeloxityX','VeloxityY','VeloxityZ']), P.integNormProduct(z, ['centers:MomentumX','centers:MomentumY','centers:MomentumZ'])]
test.testO(res,4)