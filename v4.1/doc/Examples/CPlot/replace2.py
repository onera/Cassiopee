# - replace (array) -
import Generator as G
import CPlot
import time
dt = 0.1

a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
b = G.cartTetra( (-10,0,0), (1,1,1), (10,10,10) )
c = G.cart( (0,10,0), (1,1,1), (10,10,10) )
A = [a,b,c]
CPlot.display(A); time.sleep(dt)

# Replace first zone
for i in range(5):
    b = G.cart( (i,0,0), (1,1,1), (10,10,10) )
    CPlot.replace(A, 0, b); CPlot.render(); time.sleep(dt)

# Replace second zone
for i in range(5):
    b = G.cartTetra( (-10+i,0,0), (1,1,1), (10,10,10) )
    CPlot.replace(A, 1, b); CPlot.render(); time.sleep(dt)

# Replace third zone
for i in range(5):
    b = G.cart( (i,0,0), (1,1,1), (10,10,10) )
    CPlot.replace(A, 2, b); CPlot.render(); time.sleep(dt)

# Change nature
b = G.cart( (-20,0,0), (1,1,1), (10,10,10) )
CPlot.replace(A, 1, b); CPlot.render(); time.sleep(dt)

# Change nature
b = G.cartTetra( (-20,0,0), (1,1,1), (10,10,10) )
CPlot.replace(A, 0, b); CPlot.render(); time.sleep(dt)
