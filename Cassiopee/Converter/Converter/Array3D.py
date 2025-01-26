# Conversion from/to Array / 3D Array
# Written by P. Ginibre

import numpy as N

# Convert arrays in 3D arrays
def convertArrays2Arrays3D(CArrays, VERBOSE=None):
    """Convert a standard array to a 3D array."""
    Blocks = []
    if VERBOSE:
        print(len(CArrays))
    for b in range(len(CArrays)):
        dimI = CArrays[b][2]
        dimJ = CArrays[b][3]
        dimK = CArrays[b][4]
        if VERBOSE:
            print(dimI,dimJ,dimK)
        var = CArrays[b][0].split(',')
        if VERBOSE:
            print(CArrays[b][1].shape)
        Var = [N.swapaxes(CArrays[b][1][i,:].reshape((dimK,dimJ,dimI)),0,2) for i in range(CArrays[b][1].shape[0])]
        Blocks.append([var,Var])
    return Blocks

# Convert 3D arrays in arrays
def convertArrays3D2Arrays(CArrays):
    """Convert a 3D array to a standard array."""
    a = []
    for i in CArrays:
        a.append(convertArray3D2Array(i[0], i[1]))
    return a

# Convert a 3D array in an array
def convertArray3D2Array(var, Var):
    v = ','.join(var)
    if len(var) != len(Var):
        raise Exception("Les nombres de variables et de tableaux ne correspondent pas %d <=> %d"%(len(var),len(Var)))
    Shape0=Var[0].shape
    for n in range(1,len(Var)):
        if Var[n].shape != Shape0:
            raise Exception("Dimensions Differentes des variables non autorise !")

    V = N.zeros((len(var),Var[0].size))
    for n in range(len(var)):
        #print Var[n]
        if len(Shape0) == 2:
            VV = N.ravel(N.swapaxes(Var[n].reshape(Var[n].shape[0],Var[n].shape[1],1),2,0))
        elif len(Shape0) == 3:
            VV = N.ravel(N.swapaxes(Var[n],0,2))
        else:
            VV = Var[n]
        #print VV
        V[n,:] = VV[:]
    CArray = [v, V, Shape0[0], Shape0[1], Shape0[2]]
    return CArray
