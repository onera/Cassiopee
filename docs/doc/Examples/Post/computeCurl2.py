# - computeCurl (array) -
import Converter as C
import Post as P

sol = C.convertFile2Arrays("naca.plt","bin_tp")
solp = []
varname = ['VelocityX','VelocityY','VelocityZ']
for i in sol :
    coord = C.extractVars(i,['x','y','z'])
    p = P.computeVariables(i,varname)
    p = C.addVars([coord,p])
    p = P.computeCurl(p,varname)
    coordc = P.node2Center(coord)
    p = C.addVars([coordc,p])
    solp.append(p)
C.convertArrays2File(solp, "out1.plt", "bin_tp")

# deuxieme test 2D :
import Transform as T

solp = []
for i in sol :
    i = T.subzone(i,(1,1,1),(i[2],i[3],1))
    coord = C.extractVars(i,['x','y','z'])
    p = P.computeVariables(i,varname)
    p0 = C.addVars([coord,p])
    p = P.computeCurl(p0,varname)
    coordc = P.node2Center(p0)
    p = C.addVars([coordc,p])
    solp.append(p)

C.convertArrays2File(solp, "out2.plt", "bin_tp")
#
# REFERENCE
ref = C.convertFile2Arrays("out1.plt", "bin_tp")
sol = []
for i in ref :
    i = T.subzone(i,(1,1,1),(i[2],i[3],1))
    sol.append(i)
C.convertArrays2File(sol, "ref.plt", "bin_tp")
