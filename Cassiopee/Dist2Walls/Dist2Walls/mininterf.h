      dx = xw2[indw2]-pt[0];
      dy = yw2[indw2]-pt[1];
      dz = zw2[indw2]-pt[2];
      dist = dx*dx + dy*dy + dz*dz;
      distancep[ind] = sqrt(dist);
      distmin = dist;
