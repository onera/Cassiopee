if (motionType==3){
  tmp_x=xPC[noind+ideb];
  tmp_y=yPC[noind+ideb];
  tmp_z=zPC[noind+ideb];
  
#include "IBC/motiontype3_velocities.h"
  
  uOut[indR]+=uGrid_local;
  vOut[indR]+=vGrid_local;
  wOut[indR]+=wGrid_local;
 }
