# - convertArrays2File (bin_df3) -
import Converter as C
import Generator as G

a = G.cart( (0,0,0), (1,1,1), (11,11,11) )
a = C.addVar(a, 'ro')

# Convert density field to povray density file
C.convertArrays2File([a], 'out.df3', 'bin_df3')
