# - display (array) -
# display offscreen using OSMesa parallel
import CPlot
import Geom as D
import Converter.Mpi as Cmpi
import KCore.test as test

LOCAL = test.getLocal()

a = D.sphere((0,Cmpi.rank,0), 1, N=200)

# posCam, posEye, dirCam are compulsary for parallel export (otherwise autofit)
CPlot.display(a, mode=1, solidStyle=1, offscreen=7, export=LOCAL+'/out.png',
              posCam=(-3,0,0), posEye=(0,0,0), dirCam=(0,0,1))
CPlot.finalizeExport(7)
Cmpi.barrier()
if Cmpi.rank == 0: test.testF(LOCAL+'/out.png', 1)

# malheureusement, l'anaglyph ne fonctionne pas bien en parallele
#CPlot.display(a, mode=1, offscreen=7, export=LOCAL+'/out2.png',
#              posCam=(-3,0,0), posEye=(0,0,0), dirCam=(0,0,1),
#              stereo=1, stereoDist=0.55)
#CPlot.finalizeExport(7)
#Cmpi.barrier()
#if Cmpi.rank == 0: test.testF(LOCAL+'/out2.png', 2)
