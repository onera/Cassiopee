# - extractPn (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Post.ExtraVariables2 as PE

a = G.cart((0.,0.,0.), (13./100.,13./100.,1.), (100,100,1))
for n in ['Pressure']:
    C._initVars(a, '{centers:%s} = 1.'%n)
PE._extractPn(a)
C.convertPyTree2File(a, 'out.cgns')
