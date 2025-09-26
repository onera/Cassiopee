# - integMoment (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P
import KCore.test as test

N = 20

def F1D(x,y,z): return (0*y, 9*x, 0) # int_[0,1] OM x F dS = (0,0,3/3) = (0,0,3)
def F2D(x,y,z): return (2*x, 3*y, 0) # int_[0,1]*[0,1] OM x F dS = (0,0,(3-2)/4) = (0,0,0.25)

#-------------------------------
# 1D
#-------------------------------

# STRUCT
z = G.cart((0,0,0), (1./(N-1),1,1), (N,1,1))
z = C.initVars(z, ['VelocityX', 'VelocityY', 'VelocityZ'], F1D, ['CoordinateX','CoordinateY','CoordinateZ'], isVectorized=True)
z = C.initVars(z, 'centers:MomentumX', 0.)
z = C.initVars(z, 'centers:MomentumY', 4.)
z = C.initVars(z, 'centers:MomentumZ', 0.)
res = [P.integMoment(z, center=(0.,0.,0.), vector=['VelocityX','VelocityY','VelocityZ']),
       P.integMoment(z, center=(0.,0.,0.), vector=['centers:MomentumX','centers:MomentumY','centers:MomentumZ'])]
test.testO(res,1)

# BE (BAR)
z = G.cartTetra((0,0,0), (1./(N-1),1,1), (N,1,1))
z = C.initVars(z, ['VelocityX', 'VelocityY', 'VelocityZ'], F1D, ['CoordinateX','CoordinateY','CoordinateZ'], isVectorized=True)
z = C.initVars(z, 'centers:MomentumX', 0.)
z = C.initVars(z, 'centers:MomentumY', 4.)
z = C.initVars(z, 'centers:MomentumZ', 0.)
res = [P.integMoment(z, center=(0.,0.,0.), vector=['VelocityX','VelocityY','VelocityZ']),
       P.integMoment(z, center=(0.,0.,0.), vector=['centers:MomentumX','centers:MomentumY','centers:MomentumZ'])]
test.testO(res,2)

#-------------------------------
# 2D
#-------------------------------

# STRUCT
z = G.cart((0,0,0), (1./(N-1),1./(N-1),1), (N,N,1))
z = C.initVars(z, ['VelocityX', 'VelocityY', 'VelocityZ'], F2D, ['CoordinateX','CoordinateY','CoordinateZ'], isVectorized=True)
z = C.initVars(z, 'centers:MomentumX', 4.)
z = C.initVars(z, 'centers:MomentumY', 0.)
z = C.initVars(z, 'centers:MomentumZ', 0.)
res = [P.integMoment(z, center=(0.,0.,0.), vector=['VelocityX','VelocityY','VelocityZ']),
       P.integMoment(z, center=(0.,0.,0.), vector=['centers:MomentumX','centers:MomentumY','centers:MomentumZ'])]
test.testO(res,3)

# BE (TRI)
z = G.cartTetra((0,0,0), (1./(N-1),1./(N-1),1), (N,N,1))
z = C.initVars(z, ['VelocityX', 'VelocityY', 'VelocityZ'], F2D, ['CoordinateX','CoordinateY','CoordinateZ'], isVectorized=True)
z = C.initVars(z, 'centers:MomentumX', 4.)
z = C.initVars(z, 'centers:MomentumY', 0.)
z = C.initVars(z, 'centers:MomentumZ', 0.)
res = [P.integMoment(z, center=(0.,0.,0.), vector=['VelocityX','VelocityY','VelocityZ']),
       P.integMoment(z, center=(0.,0.,0.), vector=['centers:MomentumX','centers:MomentumY','centers:MomentumZ'])]
test.testO(res,4)

# BE (QUAD)
z = G.cartHexa((0,0,0), (1./(N-1),1./(N-1),1), (N,N,1))
z = C.initVars(z, ['VelocityX', 'VelocityY', 'VelocityZ'], F2D, ['CoordinateX','CoordinateY','CoordinateZ'], isVectorized=True)
z = C.initVars(z, 'centers:MomentumX', 4.)
z = C.initVars(z, 'centers:MomentumY', 0.)
z = C.initVars(z, 'centers:MomentumZ', 0.)
res = [P.integMoment(z, center=(0.,0.,0.), vector=['VelocityX','VelocityY','VelocityZ']),
       P.integMoment(z, center=(0.,0.,0.), vector=['centers:MomentumX','centers:MomentumY','centers:MomentumZ'])]
test.testO(res,5)

# ME (QUAD+TRI)
a = G.cartHexa((0.,0.,0.), (0.5/(N-1),1./(N-1),1), (N,N,1))
b = G.cartTetra((0.5,0.,0.), (0.5/(N-1),1./(N-1),1), (N,N,1))
z = C.mergeConnectivity(a, b, boundary=0)
z = C.initVars(z, ['VelocityX', 'VelocityY', 'VelocityZ'], F2D, ['CoordinateX','CoordinateY','CoordinateZ'], isVectorized=True)
z = C.initVars(z, 'centers:MomentumX', 4.)
z = C.initVars(z, 'centers:MomentumY', 0.)
z = C.initVars(z, 'centers:MomentumZ', 0.)
res = [P.integMoment(z, center=(0.,0.,0.), vector=['VelocityX','VelocityY','VelocityZ']),
       P.integMoment(z, center=(0.,0.,0.), vector=['centers:MomentumX','centers:MomentumY','centers:MomentumZ'])]
test.testO(res,6)
