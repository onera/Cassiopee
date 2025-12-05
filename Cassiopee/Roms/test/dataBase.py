# dataBase
import Roms.DB.DataBase as DataBase

# create dataBase
db = DataBase.DataBase("BASEX", parameters=['Mach', 'AoA'])
q = db.query()
print(q)

# register without data
db.register("sample", {'Mach':0.5, 'AoA':0.0}, data=None)
db.register("sample", {'Mach':0.6, 'AoA':0.0}, data=None)
db.register("sample", {'Mach':0.7, 'AoA':1.0}, data=None)

# Create data (pytree)
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal
z = G.cart((0,0,0), (1,1,1), (10,10,10))
C._initVars(z, '{Pressure}={CoordinateX}')
t = C.newPyTree(['Base', z])
db.registerReference(t, "ref1")

# filter the relevant data before register
Internal._rmNodesFromName(t, 'GridCoordinates')

# register data (replace)
db.register("sample", {'Mach':0.5, 'AoA':0.0}, "ref1", data=t)
db.register("sample", {'Mach':0.6, 'AoA':0.0}, "ref1", data=t)
db.register("sample", {'Mach':0.7, 'AoA':1.0}, "ref1", data=t)

# query catalog
print(db.query())

# query com
q = db.query("id = 1")
print(q)
q = db.query("Mach = 0.7 AND AoA >= 0.0")
print(q)

# query point
q = db.query({'Mach':0.5, 'AoA':0.0})
print(q)

# exist point
ret = db.exist({'Mach':0.5, 'AoA':0.0})
print(ret)

# fetch tree
ts = db.fetchTree(q)[0]
#C.convertPyTree2File(ts, 'out.cgns')
Internal.printTree(ts)

# fetch param as vector
param = db.fetchParams(q)
print(param)

# fetch param as points
points = db.fetchPoints(q)
print(db.exist(points[0]))
print(points)

# fetch matrix
matrix = db.fetchMatrix(q, variables=['Pressure'])
print(matrix)

# delete query
q = db.query("Mach = 0.7 AND AoA = 1.0")
db.delete(q)
