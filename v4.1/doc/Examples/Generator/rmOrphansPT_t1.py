# - rmOrphans (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import KCore.test as test

insertPos = [3, 10, 22, 57, 91]
norphans = len(insertPos)

def _addOrphans(t, api):
    import numpy as np
    np.random.seed(42)

    gc = Internal.getNodeFromType2(t, 'GridCoordinates_t')
    n_x = Internal.getNodeFromName1(gc, 'CoordinateX')
    n_y = Internal.getNodeFromName1(gc, 'CoordinateY')
    n_z = Internal.getNodeFromName1(gc, 'CoordinateZ')
    x = n_x[1]; y = n_y[1]; z = n_z[1]
    ngon = Internal.getNGonNode(t)
    ec = Internal.getNodeFromName1(ngon, 'ElementConnectivity')[1]

    npts = x.shape[0]
    npts2 = npts + norphans
    indir = np.full(npts2, -1, dtype=Internal.E_NpyInt)
    mask = np.ones(npts2, dtype=bool)
    mask[insertPos] = False
    indir[mask] = np.arange(npts)
    origPts = indir >= 0

    # New coords: copy existing coords and add new ones
    x2 = np.empty(npts2, dtype=x.dtype)
    y2 = np.empty(npts2, dtype=y.dtype)
    z2 = np.empty(npts2, dtype=z.dtype)

    x2[mask] = x[indir[mask]]
    y2[mask] = y[indir[mask]]
    z2[mask] = z[indir[mask]]

    x2[~mask] = np.random.uniform(10., 11., norphans)
    y2[~mask] = np.random.uniform(10., 11., norphans)
    z2[~mask] = np.random.uniform(10., 11., norphans)

    Internal.setValue(n_x, x2)
    Internal.setValue(n_y, y2)
    Internal.setValue(n_z, z2)

    if api != 3:
        mask = np.zeros(ec.size, dtype=bool)
        i = 0
        while i < ec.size:
            nvert = ec[i]
            mask[i+1:i+1+nvert] = True
            i += nvert + 1
    else:
        mask = np.ones(ec.size, dtype=bool)

    invIndir = np.empty(npts, dtype=Internal.E_NpyInt)
    invIndir[indir[origPts]] = np.flatnonzero(origPts)
    ec[mask] = invIndir[ec[mask] - 1] + 1

    t[1][0][0] += norphans

# add and delete orphans - api 1
api = 1
a = G.cartNGon((0,0,0), (1,1,1), (10,10,10), api=api)
test.testT(a, 1)

_addOrphans(a, api=api)
G._rmOrphans(a)
test.testT(a, 2)  # same as ref1 but non-compact

# add and delete orphans - api 3
api = 3
a = G.cartNGon((0,0,0), (1,1,1), (10,10,10), api=api)
test.testT(a, 3)  # create orphan-free ref

_addOrphans(a, api=api)
G._rmOrphans(a)
test.testT(a, 3)
