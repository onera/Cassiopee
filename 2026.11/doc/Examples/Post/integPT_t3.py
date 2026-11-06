# - integ (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P
import KCore.test as test

N = 20

#-------------------------------
# 1D
#-------------------------------

# STRUCT
z = G.cart((0,0,0), (1./(N-1),1,1), (N,1,1))
z = C.initVars(z, 'VeloxityX', 2.)
z = C.initVars(z, 'centers:VelocityY', 4.)
res = P.integ2(z,'VeloxityX') + P.integ2(z,'centers:VelocityY')
test.testO(res,1)

# BE (BAR)
z = G.cartTetra((0,0,0), (1./(N-1),1,1), (N,1,1))
z = C.initVars(z, 'VeloxityX', 2.)
z = C.initVars(z, 'centers:VelocityY', 4.)
res = P.integ2(z,'VeloxityX') + P.integ2(z,'centers:VelocityY')
test.testO(res,2)

#-------------------------------
# 2D
#-------------------------------

# STRUCT
z = G.cart((0,0,0), (1./(N-1),1./(N-1),1), (N,N,1))
z = C.initVars(z, 'VeloxityX', 2.)
z = C.initVars(z, 'centers:VelocityY', 4.)
res = P.integ2(z,'VeloxityX') + P.integ2(z,'centers:VelocityY')
test.testO(res,3)

# BE (TRI)
z = G.cartTetra((0,0,0), (1./(N-1),1./(N-1),1), (N,N,1))
z = C.initVars(z, 'VeloxityX', 2.)
z = C.initVars(z, 'centers:VelocityY', 4.)
res = P.integ2(z,'VeloxityX') + P.integ2(z,'centers:VelocityY')
test.testO(res,4)

# BE (QUAD)
z = G.cartHexa((0,0,0), (1./(N-1),1./(N-1),1), (N,N,1))
z = C.initVars(z, 'VeloxityX', 2.)
z = C.initVars(z, 'centers:VelocityY', 4.)
res = P.integ2(z,'VeloxityX') + P.integ2(z,'centers:VelocityY')
test.testO(res,5)

# ME (QUAD+TRI)
a = G.cartHexa((0.,0.,0.), (1./(N-1),1./(N-1),1), (N,N,1))
b = G.cartTetra((1.,0.,0.), (1./(N-1),1./(N-1),1), (N,N,1))
z = C.mergeConnectivity(a, b, boundary=0)
z = C.initVars(z, 'VeloxityX', 2.)
z = C.initVars(z, 'centers:VelocityY', 4.)
res = P.integ2(z,'VeloxityX') + P.integ2(z,'centers:VelocityY')
test.testO(res,6)