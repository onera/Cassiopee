from XCore import xcore

from mpi4py import MPI
import numpy as np

comm = MPI.COMM_WORLD
nb_procs = comm.size
rank     = comm.rank

fich = open("Sortie{:02}.txt".format(rank), 'w')
splitted = None
if nb_procs == 1:
    elt2vert = np.array([1,2,3,4,3,2,4,5], np.int32)
    coord    = np.array([[0.,0.,1.,1./3.],
                         [0.,1.,0.,1./3.],
                         [1.,0.,0.,1./3.]
                         [1.,1.,1.,1./3.]], np.double)
    zone = ("TETRA", elt2vert, coord)
    zones = [([zone,], 5, 2),]
    #
    splitted = xcore.splitElements(zones)
else:
    elt2vert_glob = np.empty(8*nb_procs, np.int32)
    for i in range(0,2*nb_procs):
        elt2vert_glob[4*i:4*(i+1)] = [i+1,i+2,i+3,i+4]
        if (i%2 == 1) :
            elt2vert_glob[4*i+1],elt2vert_glob[4*i+2] = elt2vert_glob[4*i+2],elt2vert_glob[4*i+1]
    nb_vert_glob  = elt2vert_glob[-1]
    coords_glob   = np.empty((3,nb_vert_glob), np.double)
    for i in range(0,elt2vert_glob[-1]):
        coords_glob[0,i] = (i+1.)/2.
        coords_glob[1,i] = (i+2.)/2.
        coords_glob[2,i] = (i+3.)/3.
    elt2vert = np.empty(8, np.int32)
    nb_vert_loc = nb_vert_glob//nb_procs
    if rank < nb_vert_glob%nb_procs :
        nb_vert_loc += 1
    beg_vert = rank*nb_vert_loc
    if rank >= nb_vert_glob%nb_procs :
        beg_vert += nb_vert_glob%nb_procs
    coords = np.empty((3,nb_vert_loc), np.double)
    elt2vert[0:4] = elt2vert_glob[4*rank:4*(rank+1)]
    elt2vert[4:8] = elt2vert_glob[4*(rank+nb_procs):4*(rank+nb_procs+1)]
    coords[:,:]   = coords_glob[:,beg_vert:beg_vert+nb_vert_loc]
    fich.write("elt2vert : {}\n".format(elt2vert))
    fich.write("Coords ({}) =>\n{}\n".format(coords.shape, coords))
    zone = ("TETRA", elt2vert, coords)
    zones = [([zone,], elt2vert_glob[-1], 2*  nb_procs),]
    #
    splitted = xcore.splitElements(zones)
fich.write("Splitted data :\n{}".format(splitted))
fich.close()
