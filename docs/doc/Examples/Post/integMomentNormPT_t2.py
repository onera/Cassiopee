# - integMomentNorm (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P
import KCore.test as test

N = 20

def F(x,y,z): return y # int_[0,1]*[0,1] OM x (F.n) dS = (1/3, -0.25, 0)

#-------------------------------
# 2D
#-------------------------------

# STRUCT
z = G.cart((0,0,0), (1./(N-1),1./(N-1),1), (N,N,1))
z = C.initVars(z, 'VelocityX', F, ['CoordinateX', 'CoordinateY', 'CoordinateZ'])
z = C.initVars(z, 'centers:VelocityY', 2.)
res = [P.integMomentNorm(z, center=(0.,0.,0.), var='VelocityX'),
       P.integMomentNorm(z, center=(0.,0.,0.), var='centers:VelocityY')]
test.testO(res,1)

# BE (TRI)
z = G.cartTetra((0,0,0), (1./(N-1),1./(N-1),1), (N,N,1))
z = C.initVars(z, 'VelocityX', F, ['CoordinateX', 'CoordinateY', 'CoordinateZ'])
z = C.initVars(z, 'centers:VelocityY', 2.)
res = [P.integMomentNorm(z, center=(0.,0.,0.), var='VelocityX'),
       P.integMomentNorm(z, center=(0.,0.,0.), var='centers:VelocityY')]
test.testO(res,2)

# BE (QUAD)
z = G.cartHexa((0,0,0), (1./(N-1),1./(N-1),1), (N,N,1))
z = C.initVars(z, 'VelocityX', F, ['CoordinateX', 'CoordinateY', 'CoordinateZ'])
z = C.initVars(z, 'centers:VelocityY', 2.)
res = [P.integMomentNorm(z, center=(0.,0.,0.), var='VelocityX'),
       P.integMomentNorm(z, center=(0.,0.,0.), var='centers:VelocityY')]
test.testO(res,3)

# ME (QUAD+TRI)
a = G.cartHexa((0.,0.,0.), (0.5/(N-1),1./(N-1),1), (N,N,1))
b = G.cartTetra((0.5,0.,0.), (0.5/(N-1),1./(N-1),1), (N,N,1))
z = C.mergeConnectivity(a, b, boundary=0)
z = C.initVars(z, 'VelocityX', F, ['CoordinateX', 'CoordinateY', 'CoordinateZ'])
z = C.initVars(z, 'centers:VelocityY', 2.)
res = [P.integMomentNorm(z, center=(0.,0.,0.), var='VelocityX'),
       P.integMomentNorm(z, center=(0.,0.,0.), var='centers:VelocityY')]
test.testO(res,4)
