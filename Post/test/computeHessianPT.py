import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G

a = G.cartNGon((0.,0.,0.), (1./10.,1./10.,1./10.), (11,11,2))

def F(x, y): return 3.*x*x + 2.*y*y

a = C.initVars(a, 'centers:F', F, ['centers:CoordinateX', 'centers:CoordinateY'])

# dim = 2 (pur 2D), 3 (2D extrudé ou 3D)
dim = 2
a = P.computeHessian(a, 'centers:F', dim)

C.convertPyTree2File(a, 'out.cgns')
