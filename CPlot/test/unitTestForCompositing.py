# - display (array) -
# using compositing
import Generator as G
import CPlot
import os

cube   = G.cart((-0.25,0.,-0.25),(1.,1.,1.),(2,2,2))
cube2  = G.cart((0.25,0.,0.25),(0.75,0.75,0.75),(2,2,2))

res = "1024x768"

posCamZ = (2.5,0.5,0.5)
posEyeZ = (0.5,0.5,0.5)
dirCamZ = (0.,0.,1.)

posCamY = (0.5,2.5,0.5)
posEyeY = (0.5,0.5,0.5)
dirCamY = (0.,0.,1.)

posCamX = (0.5,0.5,2.5)
posEyeX = (0.5,0.5,0.5)
dirCamX = (1.,0.,0.)

posCamDiag = (2.5,2.5,2.5)
posEyeDiag = (0.5,0.5,0.5)
dirCamDiag = (1.,0.,-1.)

print("Sauvegarde cube profondeur en z")
CPlot.display([cube,cube2], posCam=posCamZ, posEye=posEyeZ, dirCam=dirCamZ, exportResolution=res, export="ref_z.png", offscreen=2)
CPlot.finalizeExport()

print("Sauvegarde cube profondeur en y")
CPlot.display([cube,cube2], posCam=posCamY, posEye=posEyeY, dirCam=dirCamY, exportResolution=res, export="ref_y.png", offscreen=2)
CPlot.finalizeExport()

print("Sauvegarde cube profondeur en diag")
CPlot.display([cube,cube2], posCam=posCamDiag, posEye=posEyeDiag, dirCam=dirCamDiag, exportResolution=res, export="ref_diag.png", offscreen=2)
CPlot.finalizeExport()

# Accumulate in images
print("Composition en profondeur z")
CPlot.display([cube], posCam=posCamZ, posEye=posEyeZ, dirCam=dirCamZ,exportResolution=res,  export="composite_z.png", offscreen=3)
CPlot.finalizeExport(3)
CPlot.display([cube2], posCam=posCamZ, posEye=posEyeZ, dirCam=dirCamZ,exportResolution=res,  export="composite_z.png", offscreen=4)
CPlot.finalizeExport(4)

# Accumulation en y
print("Composition en profondeur y")
CPlot.display([cube], posCam=posCamY, posEye=posEyeY, dirCam=dirCamY, exportResolution=res, export="composite_y.png", offscreen=3)
CPlot.finalizeExport(3)
CPlot.display([cube2], posCam=posCamY, posEye=posEyeY, dirCam=dirCamY, exportResolution=res, export="composite_y.png", offscreen=4)
CPlot.finalizeExport(4)

print("Composition en profondeur diag")
CPlot.display([cube], posCam=posCamDiag, posEye=posEyeDiag, dirCam=dirCamDiag, exportResolution=res, export="composite_diag.png", offscreen=3)
CPlot.finalizeExport(3)
CPlot.display([cube2], posCam=posCamDiag, posEye=posEyeDiag, dirCam=dirCamDiag, exportResolution=res, export="composite_diag.png", offscreen=4)
CPlot.finalizeExport(4)

try:
	import matplotlib.pyplot as plt
	import matplotlib.image as mimage

	filenames = [ ('ref_y.png', 'composite_y.png'), ('ref_z.png','composite_z.png'), ('ref_diag.png', 'composite_diag.png') ]
	nrows, ncols = 3, 2
	fig = plt.figure()
	index=  1
	for files in filenames:
		a = fig.add_subplot(nrows, ncols, index)
		img1 = mimage.imread(files[0]);
		plt.imshow(img1)
		index += 1
		img2 = mimage.imread(files[1]);
		a = fig.add_subplot(nrows, ncols, index)
		plt.imshow(img2)
		index += 1
	plt.show()
except:
	print("matplotlib not installed")
os._exit(0)