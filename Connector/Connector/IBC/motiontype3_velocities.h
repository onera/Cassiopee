uGrid_local = 0.;
vGrid_local = 0.;
wGrid_local = 0.;

//distance
cmx = tmp_x - axispntPtrX[noind+ideb];
cmy = tmp_y - axispntPtrY[noind+ideb];
cmz = tmp_z - axispntPtrZ[noind+ideb];

// k x cm
kvcmx =  axisvecPtrY[noind+ideb]*cmz-axisvecPtrZ[noind+ideb]*cmy;
kvcmy =  axisvecPtrZ[noind+ideb]*cmx-axisvecPtrX[noind+ideb]*cmz;
kvcmz =  axisvecPtrX[noind+ideb]*cmy-axisvecPtrY[noind+ideb]*cmx;

uGrid_local = transpeedPtrX[noind+ideb] + omgPtr[noind+ideb]*kvcmx;
vGrid_local = transpeedPtrY[noind+ideb] + omgPtr[noind+ideb]*kvcmy;
wGrid_local = transpeedPtrZ[noind+ideb] + omgPtr[noind+ideb]*kvcmz;


