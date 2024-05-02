#include "proto.h"

void compute_face_centers_and_areas
(
  K_FLD::FldArrayI &cn, E_Float *x, E_Float *y, E_Float *z,
  E_Float *fcenters, E_Float *fareas
)
{
  E_Int nfaces = cn.getNFaces();
  E_Int *ngon = cn.getNGon();
  E_Int *indPG = cn.getIndPG();

  for (E_Int i = 0; i < nfaces; i++) {
    E_Int np = -1;
    E_Int *pn = cn.getFace(i, np, ngon, indPG);
    K_METRIC::compute_face_center_and_area(i, np, pn, x, y, z,
      &fcenters[3*i], &fareas[3*i]);
  }
}

void compute_cell_centers_and_vols
(
  K_FLD::FldArrayI &cn, E_Float *x, E_Float *y, E_Float *z,
  E_Int *owner, E_Int *neigh, E_Float *fcenters, E_Float *fareas,
  E_Float *centers, E_Float *volumes
)
{
  E_Int *nface = cn.getNFace();
  E_Int *indPH = cn.getIndPH();
  E_Int nfaces = cn.getNFaces();
  E_Int ncells = cn.getNElts();
  
  E_Float *vols;
  if (volumes) vols = volumes;
  else vols = (E_Float *)XMALLOC(ncells * sizeof(E_Float)); 

  // Estimate cell centers as average of face centers
  std::vector<E_Float> cEst(3*ncells, 0);
  for (E_Int i = 0; i < nfaces; i++) {
    E_Float *fc = &fcenters[3*i];
    E_Float *cE;
    
    E_Int own = owner[i];
    cE = &cEst[3*own];

    for (E_Int j = 0; j < 3; j++) cE[j] += fc[j];

    E_Int nei = neigh[i];
    if (nei == -1) continue;
    
    cE= &cEst[3*nei];

    for (E_Int j = 0; j < 3; j++) cE[j] += fc[j];
  }

  for (E_Int i = 0; i < ncells; i++) {
    E_Int stride = -1;
    cn.getElt(i, stride, nface, indPH);
    E_Float *cE = &cEst[3*i];
    for (E_Int j = 0; j < 3; j++)
      cE[j] /= stride;
  }

  memset(vols, 0, ncells*sizeof(E_Float));
  memset(centers, 0, 3*ncells*sizeof(E_Float));

  for (E_Int i = 0; i < nfaces; i++) {
    E_Float *fa = &fareas[3*i];
    E_Float *fc = &fcenters[3*i];
    E_Float pyr3vol, pc[3], *cE, d[3], *cc;

    E_Int own = owner[i];
    cE = &cEst[3*own];
    cc = &centers[3*own];
    for (E_Int j = 0; j < 3; j++) d[j] = fc[j]-cE[j];
    pyr3vol = dot(fa, d, 3);
    for (E_Int j = 0; j < 3; j++) pc[j]= 0.75*fc[j] + 0.25*cE[j];
    for (E_Int j = 0; j < 3; j++) cc[j] += pyr3vol*pc[j];
    vols[own] += pyr3vol;

    E_Int nei = neigh[i];
    if (nei == -1) continue;

    cE = &cEst[3*nei];
    cc = &centers[3*nei];
    for (E_Int j = 0; j < 3; j++) d[j] = cE[j]-fc[j];
    pyr3vol = dot(fa, d, 3);
    for (E_Int j = 0; j < 3; j++) pc[j]= 0.75*fc[j] + 0.25*cE[j];
    for (E_Int j = 0; j < 3; j++) cc[j] += pyr3vol*pc[j];
    vols[nei] += pyr3vol;
  }

  for (E_Int i = 0; i < ncells; i++) {
    E_Float *ptr = &centers[3*i];
    for (E_Int j = 0; j < 3; j++) ptr[j] /= vols[i];
    vols[i] /= 3.0;
    assert(vols[i] > 0.0);
  }

  if (!volumes) XFREE(vols);
}
