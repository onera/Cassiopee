# - fly mode autour d'un objet -
import Converter as C
import Generator as G
import CPlot

sigma = 10; beta = 8/3; ro = 28

# Vole jusqu'a la position suivante
def fly(a, bb, posCam, posEye, dirCam, step):
    CPlot.display(a, posCam=posCam, posEye=posEye, dirCam=dirCam)
    xc = 0.5*(bb[3]+bb[0])
    yc = 0.5*(bb[4]+bb[1])
    zc = 0.5*(bb[5]+bb[2])
    x = posCam[0]-xc; y = posCam[2]-yc; z = posCam[1]-zc
    x = 10*x /(bb[3]-bb[0]); y = 10*y /(bb[4]-bb[1]); z = 30*z /(bb[5]-bb[2]); 
    xp = x + step*sigma*(y - x)
    yp = y + step*(ro*x - y - x*z)
    zp = z + step*(x*y - beta*z)
    xp = xp * (bb[3]-bb[0])*0.1 + xc
    yp = yp * (bb[4]-bb[1])*0.1 + yc
    zp = zp * (bb[5]-bb[2])/30. + zc
    posCam[0] = xp
    posCam[1] = zp
    posCam[2] = yp
    return posCam

#a = C.convertFile2Arrays('fontaine.plt', 'bin_tp')
a = C.convertFile2Arrays('dauphin2.plt', 'bin_tp')
#a = C.convertFile2Arrays('interceptor.plt', 'bin_tp')
#a = C.convertFile2Arrays('tiep.plt', 'bin_tp')
#a = T.rotate(a, (0,0,0),(1,0,0),-90.)
CPlot.display(a, win=(700,500), displayBB=0, displayInfo=0,
              mode=1, solidStyle=0)

bb = G.bbox(a)
xc = 0.5*(bb[3]+bb[0])
yc = 0.5*(bb[4]+bb[1])
zc = 0.5*(bb[5]+bb[2])
posCam = [bb[0],bb[1],bb[2]]
posEye = [xc,yc,zc]
dirCam = [0,0,1]
i = 0
while (i != -1):
    posCam = fly(a, bb, posCam, posEye, dirCam, 0.003)
    i = i+1
