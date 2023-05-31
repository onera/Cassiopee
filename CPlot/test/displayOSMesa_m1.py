# - display (array) -
# display offscreen using OSMesa parallel
import CPlot
import Transform as T
import Geom as D
import Converter.Mpi as Cmpi
import KCore.config
import KCore.test as test

if not KCore.config.CPlotOffScreen: 
    print("CPlot: no offscreen osmesa installed.")
    import sys; sys.exit()

a = D.sphere((0,Cmpi.rank,0), 1, N=200)

# posCam, posEye, dirCam are compulosary for parallel export (otherwise autofit)
CPlot.display(a, mode=1, offscreen=7, export='out.png', posCam=(-3,0,0), posEye=(0,0,0), dirCam=(0,0,1))
CPlot.finalizeExport(7)
Cmpi.barrier()
if Cmpi.rank == 0: test.testF('out.png', 1)

# idem with anaglyph
CPlot.display(a, mode=1, offscreen=7, export='out2.png',
              posCam=(-3,0,0), posEye=(0,0,0), dirCam=(0,0,1),
              stereo=1, stereoDist=0.55)
CPlot.finalizeExport(7)
Cmpi.barrier()
if Cmpi.rank == 0: test.testF('out2.png', 2)
