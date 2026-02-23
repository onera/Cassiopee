# - getUV (array) -
import Geom as D
import Generator as G
import Converter as C
import Post as P

a = G.cartTetra((0,0,0), (1,1,1), (10,10,10))
a = P.exteriorFaces(a)
a = C.initVars(a, '{VelocityX} = {x}')
a = C.initVars(a, '{VelocityY} = 0.1')
a = C.initVars(a, '{VelocityZ} = 0.')

(a, color, normal) = D.getUV(a, 2., 1920, fields=['VelocityX'])

# model with uv
C.convertArrays2File(a, 'out.plt')
# Texture Image in uv space
C.convertArrays2File(color, 'color.png')
# Bump map in uv space
C.convertArrays2File(normal, 'bump.png')
