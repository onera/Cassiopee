# - loadSkeleton -
import Converter.Filter as Filter
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Mpi as CMpi
import sys

# Create test case
#fich = open(f'sortie{CMpi.rank}.txt', mode='w')
if (CMpi.rank == 0):
    a = G.cartHexa((0,0,0), (1,1,1), (8,8,8))
    #fich.write(f"a={a}\n")
    #print(f"a={a}", flush=True)
    C.convertPyTree2File(a, 'file.hdf')
CMpi.barrier()
# Create a handle on a CGNS file
h = Filter.Handle('file.hdf')

# Load skeleton
#print("splitted and distribute",flush=True)

a = h.distributedLoadAndSplitSkeleton()
#fich.write(f"a = {a}\n")
#fich.close()
#print(f"a = {a}", flush=True)
C.convertPyTree2File(a, f'file{CMpi.rank}.hdf')
