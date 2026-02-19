# include <math.h>
# include <stdio.h>
# include "CompGeom/compGeom.h"

// ============================================================================
// Make a rotation of a mesh 
// IN: dim          ! mesh size
// IN: teta         ! angle
// IN: center       ! center of rotation
// IN: axis         ! rotation Vector
// OUT: xo, yo, zo  ! rotated mesh point
// ============================================================================
void K_COMPGEOM::rotateMesh(
  const E_Int dim, const E_Float teta,
  const E_Float* center, const E_Float* axis, 
  E_Float* xo, E_Float* yo, E_Float* zo
)
{
  E_Float unx, uny, unz;
  E_Float norm;
  E_Float e0, e1, e2, e3;
  E_Float a1, steta, stetas2;

  // nx,ny,nz must be unit vector
  norm = axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2];
  if (norm <= 1.e-12)
  {
    printf("Error: rotateMesh: nx,ny,nz has null norm.\n");
    return;
  }

  norm = 1.0/sqrt(norm);
  unx = axis[0]*norm;
  uny = axis[1]*norm;
  unz = axis[2]*norm;

  steta = sin(teta);
  stetas2 = sin(0.5*teta);

  // quaternion
  e0 = cos(0.5*teta);
  e1 = -unx*stetas2;
  e2 = -uny*stetas2;
  e3 = -unz*stetas2;
  a1 = e0*e0 - e1*e1 - e2*e2 - e3*e3;

  #pragma omp parallel
  {
    E_Float rx, ry, rz;
    E_Float a2;
    E_Float px, py, pz;

    #pragma omp for 
    for (E_Int ind = 0; ind < dim; ind++)
    {
      rx = xo[ind] - center[0];
      ry = yo[ind] - center[1];
      rz = zo[ind] - center[2];

      a2 = e1*rx + e2*ry + e3*rz;
      px = a1*rx + 2*e1*a2 - (ry*unz - rz*uny)*steta;
      py = a1*ry + 2*e2*a2 - (rz*unx - rx*unz)*steta;
      pz = a1*rz + 2*e3*a2 - (rx*uny - ry*unx)*steta;

      xo[ind] = center[0] + px;
      yo[ind] = center[1] + py;
      zo[ind] = center[2] + pz;
    }
  }
}

// ============================================================================
// Same function as previous, but different interface 
// IN: npts            ! mesh size
// IN: teta            ! angle
// IN: x, y, z         ! mesh coordinates
// IN: xc, yc, zc      ! center of rotation
// IN: nx, ny, nz      ! rotation vector  -- nx,ny,nz must be unit vector
// OUT: xo, yo, zo     ! rotated mesh
// ============================================================================
void K_COMPGEOM::rotateMesh2(
  const E_Int npts, const E_Float teta,
  const E_Float xc, const E_Float yc, const E_Float zc,
  const E_Float nx, const E_Float ny, const E_Float nz,
  const E_Float* x, const E_Float* y, const E_Float* z,
  E_Float* xo, E_Float* yo, E_Float* zo
)
{
  E_Float unx, uny, unz;
  E_Float norm;
  E_Float e0, e1, e2, e3;
  E_Float a1, sinteta, sinteta5;

  // nx,ny,nz must be unit vector
  norm = nx*nx + ny*ny + nz*nz;
  if (norm <= 1.e-12)
  {
    printf("Error: rotate: nx,ny,nz has null norm.\n");
    return;
  }

  norm = 1.0/sqrt(norm);
  unx = nx*norm;
  uny = ny*norm;
  unz = nz*norm;

  sinteta = sin(teta);
  sinteta5 = sin(0.5*teta);

  // quaternion
  e0 = cos(0.5*teta);
  e1 = -unx*sinteta5;
  e2 = -uny*sinteta5;
  e3 = -unz*sinteta5;
  a1 = e0*e0 - e1*e1 - e2*e2 - e3*e3;
  
  #pragma omp parallel
  {
    E_Float rx, ry, rz;
    E_Float a2;
    E_Float px, py, pz;
    #pragma omp for 
    for (E_Int ind = 0; ind < npts; ind++)
    {
      rx = x[ind] - xc;
      ry = y[ind] - yc;
      rz = z[ind] - zc;
        
      a2 = e1*rx + e2*ry + e3*rz;
      px = a1*rx + 2*e1*a2 - (ry*unz - rz*uny)*sinteta;
      py = a1*ry + 2*e2*a2 - (rz*unx - rx*unz)*sinteta;
      pz = a1*rz + 2*e3*a2 - (rx*uny - ry*unx)*sinteta;
        
      xo[ind] = xc + px;
      yo[ind] = yc + py;
      zo[ind] = zc + pz;
    }
  }
}
