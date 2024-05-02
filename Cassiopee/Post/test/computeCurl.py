# - computeCurl (array) -
import Converter as C
import Post as P
import Generator as G

def F(x,y,z):
    return 12*y*y + 4

ni = 30; nj = 40; nk = 3
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
m = C.initVars(m,'F1',F,['x','y','z'])
m = C.initVars(m,'F2',0.); m = C.initVars(m,'F3',0.)

varname = ['F1','F2','F3']
p = P.computeCurl(m, varname) # defined on centers
p = C.center2Node(p) # back on init grid
p = C.addVars([m,p])
C.convertArrays2File([p], "out.plt")
