# - all Hubs -
import Geom as D
import Transform as T
import Generator as G
import Post as P

#==============================================================================
# Hub sans mat correspondant a un petit cylindre
# IN: R: rayon du hub
# IN: zlow: z point bas
# IN: zup: z point haut
# IN: zcone: z du haut de la partie haute
#==============================================================================
def hub1(R, zlow, zup, zcone, NbPoints):
    """Cylindrical hub with a bump."""
    PtBasHub = (R, 0., zlow)
    PtHautHub = (R, 0., zup)
    PtHautCone = (0, 0, zcone)
    l1 = D.line((0.,0.,zlow), PtBasHub)
    l2 = D.line(PtBasHub, PtHautHub, N=10)
    PtC = (R, 0., zcone)
    c = D.polyline([PtHautHub,PtC,PtHautCone])
    l3 = D.spline(c, order=3, N=100)
    l4 = T.join([l1,l2,l3])
    l4 = D.uniformize(l4, NbPoints)
    Ntheta = int(NbPoints*0.7)
    if (Ntheta//2)*2 == Ntheta: Ntheta += 1
    l4 = D.axisym(l4, (0.,0.,0.), (0.,0.,1.), Ntheta=Ntheta)
    div = NbPoints//8
    l5 = T.subzone(l4,(1,1,1),(-div,1,-1))
    l5 = T.subzone(l5,(div,1,1),(-1,1,-1))
    l6 = P.exteriorFacesStructured(l5)
    z1 = l6[0]
    z2 = l6[1]
    z1 = G.TFIO(z1)
    z2 = G.TFIO(z2)
    proj1 = T.projectOrthoSmooth(z1, [l4], niter=3)
    proj2 = T.projectOrthoSmooth(z2, [l4], niter=3)
    return [l5]+proj1+proj2
