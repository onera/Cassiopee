# include <math.h>
# include <stdio.h>
# include "CompGeom/compGeom.h"

// ============================================================================
// Make a rotation of a mesh 
// ============================================================================
//dim        ! [IN]  mesh size
//teta       ! [IN]  angle
//center     ! [IN]  center of rotation
//axis       ! [IN]  rotation Vector
//x0, y0, z0 ! [OUT] rotated mesh point
void K_COMPGEOM::rotateMesh(const E_Int dim, const E_Float teta,
			    const E_Float* center, const E_Float* axis, 
			    E_Float* x0, E_Float* y0, E_Float* z0)
{
  // LOCAL
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
#pragma omp parallel default(shared)
  {
    E_Float rx, ry, rz;
    E_Float a2;
    E_Float px, py, pz;
#pragma omp for 
    for (E_Int ind = 0; ind < dim; ind++)
      {
	rx = x0[ind] - center[0];
	ry = y0[ind] - center[1];
	rz = z0[ind] - center[2];

	a2 = e1*rx + e2*ry + e3*rz;
	px = a1*rx + 2*e1*a2 - (ry*unz - rz*uny)*steta;
	py = a1*ry + 2*e2*a2 - (rz*unx - rx*unz)*steta;
	pz = a1*rz + 2*e3*a2 - (rx*uny - ry*unx)*steta;

	x0[ind] = center[0] + px;
	y0[ind] = center[1] + py;
	z0[ind] = center[2] + pz;
      }
  }
}

// ============================================================================
// Same function as previous, but different interface 
// ============================================================================
//npts            ! [IN]  mesh size
//teta            ! [IN]  angle
//x, y, z         ! [IN]  mesh coordinates
//xc, yc, zc      ! [IN]  center of rotation
//nx, ny, nz      ! [IN]  rotation vector  -- nx,ny,nz must be unit vector
//x0, y0, z0      ! [OUT] rotated mesh
void K_COMPGEOM::rotateMesh2(const E_Int npts, const E_Float teta,
			     const E_Float xc, const E_Float yc, const E_Float zc,
			     const E_Float nx, const E_Float ny, const E_Float nz,
			     const E_Float* x, const E_Float* y, const E_Float* z,
			     E_Float* x0, E_Float* y0, E_Float* z0)
{
  // LOCAL
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
  
#pragma omp parallel default(shared)
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
	  
	x0[ind] = xc + px;
	y0[ind] = yc + py;
	z0[ind] = zc + pz;
      }
  }
}
