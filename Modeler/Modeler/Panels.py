# - All panels -
import Geom as D
import Transform as T
import Generator as G
import Converter as C
import Boxes
import Cylinders

#==============================================================================
# Panneau avec un texte
# text: le texte a afficher (\n possible)
# scale: l'echelle
# h: hauteur du poteau
#==============================================================================
def panel(text="Hello", scale=0.1, h=1.):
    t = D.text2D(text, font='vera')
    t = T.rotate(t, (0,0,0), (1,0,0), 90.)
    t = T.scale(t, scale)
    bb = G.bbox(t)
    t = T.translate(t, (0,0,-bb[2]+h))
    bb = G.bbox(t)
    dt = 0.3
    box1 = Boxes.box((bb[0]-2*dt,bb[1]+0.001,bb[2]-2*dt), (bb[3]+2*dt,bb[4]+2*dt,bb[5]+2*dt), chamfer=0.25*dt)
    box1 = C.convertArray2Tetra(box1)
    bb = G.bbox(box1)

    dx = bb[3]-bb[0]; dy = bb[4]-bb[1]; dz = bb[5]-bb[2]
    poteau = Cylinders.cylinder(R1=0.75*dt,R2=0.75*dt,h=h)
    poteau = T.translate(poteau, (bb[0]+dx*0.5,bb[1]+dy*0.5,bb[2]-h))
    o = T.join([t,box1,poteau])
    o = G.close(o)
    return o
