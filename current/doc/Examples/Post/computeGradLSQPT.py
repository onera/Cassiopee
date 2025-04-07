import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G
import Converter.Internal as I
import KCore.test as test

a = G.cartHexa((0.,0.,0.), (1./10.,1./10.,1./10.), (11,11,3))
a = C.convertArray2NGon(a)
I._adaptNGon32NGon4(a)

def f(x, y, z): return 3.*x + 2.*y + z
def g(x, y, z): return 4.*x + 3.*y + 2.*z

a = C.initVars(a, 'centers:f', f, ['centers:CoordinateX', 'centers:CoordinateY', 'centers:CoordinateZ'])
a = C.initVars(a, 'centers:g', g, ['centers:CoordinateX', 'centers:CoordinateY', 'centers:CoordinateZ'])

a = C.makeParentElements(a)
a = P.computeGradLSQ(a, ['g', 'f'])
C.convertPyTree2File(a, 'out.cgns')
