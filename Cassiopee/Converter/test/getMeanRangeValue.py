# - getMeanRangeValue (array) -
import Converter as C
import Generator as G

a = G.cart((0,0,0), (1.,1.,1.), (11,1,1))
# return the mean of the 30% smallest values
meanval = C.getMeanRangeValue(a, 'x', 0., 0.3); print(meanval)
# >> 0.75
