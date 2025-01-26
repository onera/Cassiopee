"""All modles of walls."""
import Geom as D
import Transform as T
import Generator as G
import Converter as C
from . import Boxes
import math

#==============================================================================
# IN: line: une courbe dans le plan (x,y)
# IN: Bx,By,Bz: taille des briques
# IN: nlayers: nbre de couches
# IN: chamfer: chanfrein
#==============================================================================
def wall(line,Bx,By,Bz,nlayers=1,chamfer=-1., shrink=1.):
    line = C.convertBAR2Struct(line)

    l = D.getLength(line)
    Nb = int(l/(Bx*shrink))+1
    distrib = G.cart((0,0,0), (1./(Nb-1),1,1), (Nb,1,1))
    line = G.map(line, distrib)

    posZ = 0.
    bx = Boxes.box((0,0,0), (Bx,By,Bz), chamfer)
    hbx = Boxes.box((0,0,0), (Bx*0.5,By,Bz), chamfer) # half box
    bricks = []
    for n in range(nlayers):
        for i in range(Nb-1):
            [x,y,z] = C.getValue(line, i)
            [xp,yp,zp] = C.getValue(line, i+1)

            if n%2 == 0:
                xi = 0.5*(x+xp); yi = 0.5*(y+yp)
                bx2 = T.translate(bx, (xi,yi,posZ))
                if xp-x>1.e-16: alpha = math.atan((yp-y)/(xp-x))
                else: alpha = math.pi*0.5
                bx2 = T.rotate(bx2, (xi,yi,posZ), (0,0,1), alpha*180/math.pi)
                bricks.append(bx2)
            else:
                if i > 0: [xm,ym,zm] = C.getValue(line, i-1)
                else: [xm,ym,zm] = C.getValue(line, 0)
                bx2 = T.translate(bx, (x,y,posZ))
                if xp-xm>1.e-12: alpha = math.atan((yp-ym)/(xp-xm))
                else: alpha = math.pi*0.5
                bx2 = T.rotate(bx2, (x,y,posZ), (0,0,1), alpha*180/math.pi)
                bricks.append(bx2)
        posZ += Bz
    o = T.join(bricks)
    o = G.close(o)
    return o
