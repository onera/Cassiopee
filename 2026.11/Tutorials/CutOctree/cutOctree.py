# cut an octree along X/Y/Z planes
#
import Converter.PyTree as C
import Transform.PyTree as T
import Generator.PyTree as G
import Geom.PyTree as D
import Post.PyTree as P
import Converter.Internal as Internal

# dictionary of the X/Y/Z planes
dico = {'CoordinateX':2.,'CoordinateZ':1.}

a = D.sphere((0,0,0), 1., 20); a = C.convertArray2Tetra(a); a = G.close(a)
surf = C.newPyTree(['SURF',2,a])

snear=[0.5]

# symetrize the bodies wrt. X/Y/Z planes
for coord in dico:
    for zone in Internal.getByType(surf,'Zone_t')[2]:
        if coord=='CoordinateX': surf[2][1][2].append(T.symetrize(zone,(dico[coord],0.,0.),(0,1,0),(0,0,1)))
        if coord=='CoordinateY': surf[2][1][2].append(T.symetrize(zone,(0.,dico[coord],0.),(1,0,0),(0,0,1)))
        if coord=='CoordinateZ': surf[2][1][2].append(T.symetrize(zone,(0.,0.,dico[coord]),(1,0,0),(0,1,0)))
    snear+=snear

octrbef = G.octree(surf, snear, dfar=5.)
octr = octrbef
# select cells below X/Y/Z planes (can be modified for cells beyond planes)
for coord in dico:
    octr = P.selectCells(octr, '{%s}<=%s'%(coord,dico[coord]), strict=1)
    
octr = G.octree2Struct(octr,vmin=11,ext=0,optimized=0,merged=1,sizeMax=250000)
carttree = C.newPyTree(['Cart',Internal.getZones(octr),
                        'CartBefCut',Internal.getZones(octrbef)])

C.convertPyTree2File(carttree,'cart.cgns')
