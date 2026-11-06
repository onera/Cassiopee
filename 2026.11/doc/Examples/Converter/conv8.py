# - convertFile2Arrays -
# Read structured blocks from a file
# Read unstructured blocks from a file
# Write a global file with all blocks
import Converter as C

arrays = C.convertFile2Arrays('dauphin.plt')
unsArrays = C.convertFile2Arrays('Data/unstr.plt')
arrays = arrays + unsArrays
C.convertArrays2File(arrays, 'out.plt')
