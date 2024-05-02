# - nearestFaces (array) -
# doit etre exact
import Converter as C
import Generator as G
import numpy as np

eps = 1.e-30
N = 10
# structure
a = G.cart((1E8,0,0), (1E-6,1.e-6,1), (N,N,N))
b = G.cart((1E8,0,0), (1E-6,1.e-6,1), (N,N,N))
pert_x = 1.E-7*np.random.random((N*N*N))
pert_y = 1.E-7*np.random.random((N*N*N))
pert_z = 1.E-1*np.random.random((N*N*N))
a[1][0,:] += pert_x
a[1][1,:] += pert_y
a[1][2,:] += pert_z
b[1][:,:] = a[1][:,:]
hook = C.createHook(a, function='faceCenters')
faces = C.identifyFaces(hook,b,tol=eps)
ret = faces[faces<0].shape[0]
#if ret>0: print '1. identifyFaces (Struct) FAILED: FPU is not correct [Struct/FACES]. Check compilation options.'

# HEXA
a = G.cartHexa((1E8,0,0), (1E-6,1.e-6,1), (N,N,N))
b = G.cartNGon((1E8,0,0), (1E-6,1.e-6,1), (N,N,N))
pert_x = 1.E-7*np.random.random((N*N*N))
pert_y = 1.E-7*np.random.random((N*N*N))
pert_z = 1.E-1*np.random.random((N*N*N))
a[1][0,:] += pert_x
a[1][1,:] += pert_y
a[1][2,:] += pert_z
b[1][:,:] = a[1][:,:]
hook = C.createHook(a, function='faceCenters')
faces = C.identifyFaces(hook,b,tol=eps)
ret = faces[faces<0].shape[0]
#if ret>0: print 'identifyFaces (BE) FAILED: FPU is not correct [BE/FACES]. Check compilation options.'

# NGON
a = G.cartNGon((1E8,0,0), (1E-6,1.e-6,1), (N,N,N))
pert_x = 1.E-7*np.random.random((N*N*N))
pert_y = 1.E-7*np.random.random((N*N*N))
pert_z = 1.E-1*np.random.random((N*N*N))
a[1][0,:] += pert_x; a[1][1,:] += pert_y; a[1][2,:] += pert_z
b = a[:]
hook = C.createHook(a, function='faceCenters')
faces = C.identifyFaces(hook,b,tol=eps)
ret = faces[faces<0].shape[0]
#if ret >0: print 'identifyFaces (NGON) FAILED: FPU is not correct [NGON/FACES]. Check compilation options.'
print('done.')
