import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G

a = G.cartNGon((0.,0.,0.), (1./10.,1./10.,1./10.), (11,11,3))

def F(x, y, z): return 3.*x*x + 2.*y*y + z

a = C.initVars(a, 'centers:F', F, ['centers:CoordinateX', 'centers:CoordinateY', 'centers:CoordinateZ'])

# dim = 2 (pur 2D), 3 (2D extrudé ou 3D)
dim = 3
a = P.computeGradLSQ(a, 'centers:F', dim)

C.convertPyTree2File(a, 'out.cgns')
