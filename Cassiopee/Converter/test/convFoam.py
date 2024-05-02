# - convertArrays2File (Foam format) -
import Converter as C
import Generator as G

# Write NGon file
a = G.cartNGon((0,0,0), (1,1,1), (21,21,2))

# Write mesh file
C.convertArrays2File(a, "out.foam", "fmt_foam")

# Reread
b = C.convertFile2Arrays("out.foam", "fmt_foam")
print(b)

# Rewrite in plt format
C.convertArrays2File(b, 'foamNGon.plt')
