// 0: algo orig (ortho)
// 1: mixte
// 2: au centre
#define ALGO 0

/* algo de distance ortho */

    indicesBB.clear(); candidates.clear();
    distmin = distancep[ind];

    // Point P
    pt[0] = xt[ind]; pt[1] = yt[ind]; pt[2] = zt[ind];
    // recherche du sommet P' des parois le plus proche de P
    //indw2 = kdt.getClosest(pt);
    //rx = xw2[indw2]-pt[0]; ry = yw2[indw2]-pt[1]; rz = zw2[indw2]-pt[2];
    //dist = rx*rx + ry*ry + rz*rz;

    indw2 = kdt.getClosest(pt, dist);
    if (indw2 == IDX_NONE) indw2 = 0;
    rx = xw2[indw2]-pt[0]; ry = yw2[indw2]-pt[1]; rz = zw2[indw2]-pt[2];

    rad = sqrt(dist);
    rmax = lmaxp[indw2];


    E_Int isDoOrtho=1;
    if (isminortho==1)
    {
    if (dist < distmin) { distancep[ind] = dist; distmin = dist; }
    E_Float h_target;
	E_Float d_target=dTarget;
	E_Float dhx =abs(xt[ind+1]-xt[ind]);
	E_Float dhx2=abs(xt[ind-1]-xt[ind]);
	E_Float dhy =abs(yt[ind+1]-yt[ind]);
	E_Float dhy2=abs(yt[ind-1]-yt[ind]);
	E_Float dhz =abs(zt[ind+1]-zt[ind]);
	E_Float dhz2=abs(zt[ind-1]-zt[ind]);
	if (ind==npts-1) { h_target= max({dhx2,dhy2,dhz2}); }
	else if (ind==0) { h_target= max({dhx,dhy,dhz}); }
	else { h_target= max({min(dhx,dhx2),min(dhy,dhy2),min(dhz,dhz2)}); }
	if (isIBM_F1 == 1) { d_target = sqrt(3)*4*h_target; }
	else
    {
	  if (d_target>999) { d_target = sqrt(3)*10*h_target; }
	}
#include "mininterf_ortho.h"
        if (ret != -1){
          dx = xp_local-pt[0]; dy = yp_local-pt[1]; dz = zp_local-pt[2];
          dist = dx*dx + dy*dy + dz*dz;
          if (dist < distmin) { distancep[ind] = dist; distmin = dist; }
        }
	if (sqrt(distmin)> d_target) isDoOrtho=0;
    }


    // distance au centre le plus proche
#if ALGO == 2
    if (dist < distmin) { distancep[ind] = dist; distmin = dist; }
#endif

    // distance mixte
#if ALGO == 1
    if (rad > 2*rmax)
    {
        if (dist < distmin) { distancep[ind] = dist; distmin = dist; }
    }
    else if (rad > rmax)
    {
        // find the right wall
        ret = K_COMPGEOM::projectOrthoPrecond(pt[0], pt[1], pt[2], xw, yw, zw,
                                              candidates, cnloc, xp, yp, zp,
                                              p0, p1, p2, p);
        if (ret != -1)
        {
          dx = xp-pt[0]; dy = yp-pt[1]; dz = zp-pt[2];
          dist = dx*dx + dy*dy + dz*dz;
          if (dist < distmin) { distancep[ind] = dist; distmin = dist; }
        }
    }
    else
    {
        // calcul de la bounding box de la sphere de rayon PP'
        A = 1./(10.*rmax);
        rad2 = exp(-A*rad);
        alpha = 1.-rad2;
        R = rad*rad2;
        xQ = pt[0] + alpha*rx;
        yQ = pt[1] + alpha*ry;
        zQ = pt[2] + alpha*rz;
        minB[0] = xQ-R; minB[1] = yQ-R; minB[2] = zQ-R;
        maxB[0] = xQ+R; maxB[1] = yQ+R; maxB[2] = zQ+R;

        // calcul des cellules intersectantes
        for (E_Int now = 0; now < nwalls; now++)
        {
            indicesBB.clear(); candidates.clear();
            K_SEARCH::BbTree3D* bbtree = vectOfBBTrees[now];
            bbtree->getOverlappingBoxes(minB, maxB, indicesBB);
            FldArrayF* fieldv = fieldsw[now];
            E_Int posxw = posxv[now]; E_Int posyw = posyv[now]; E_Int poszw = poszv[now]; E_Int poscw = poscv[now];
            E_Float* xw = fieldv->begin(posxw);
            E_Float* yw = fieldv->begin(posyw);
            E_Float* zw = fieldv->begin(poszw);
            E_Float* cellnw = fieldv->begin(poscw);
            FldArrayI& cnloc = *cntw[now];
            E_Int nbb = indicesBB.size();
            E_Int nvert = cnloc.getNfld();
            E_Float prodCellN2 = pow(2.,nvert);
            for (E_Int i = 0; i < nbb; i++)
            {
                et = indicesBB[i];
                prod = 1.;
                for (E_Int novert = 1; novert <= nvert; novert++)
                {
                    ind10 = cnloc(et, novert)-1;
                    prod = prod*cellnw[ind10];
                }
                if (prod != 0. && prod != prodCellN2) candidates.push_back(et);
            }
            ret = K_COMPGEOM::projectOrthoPrecond(pt[0], pt[1], pt[2], xw, yw, zw,
                                                candidates, cnloc, xp, yp, zp,
                                                p0, p1, p2, p);
            if (ret != -1)
            {
                dx = xp-pt[0]; dy = yp-pt[1]; dz = zp-pt[2];
                dist = dx*dx + dy*dy + dz*dz;
                if (dist < distmin) { distancep[ind] = dist; distmin = dist; }
            }
        }
    }
#endif

    // distance ortho
#if ALGO == 0
    if (isDoOrtho==1){
        // calcul de la bounding box de la sphere de rayon PP'
        A = 1./(10.*rmax);
        rad2 = exp(-A*rad);
        alpha = 1.-rad2;
        R = rad*rad2;

        // calcul de la bounding box de la sphere de rayon PP'
        xQ = pt[0] + alpha*rx;
        yQ = pt[1] + alpha*ry;
        zQ = pt[2] + alpha*rz;
        minB[0] = xQ-R; minB[1] = yQ-R; minB[2] = zQ-R;
        maxB[0] = xQ+R; maxB[1] = yQ+R; maxB[2] = zQ+R;
        //if (fabs(pt[1])<1.e-10) { printf("%f %f R=%f, delta=%f\n",pt[0],pt[2],R,rad); }

        if (dist < distmin) { distancep[ind] = dist; distmin = dist; }

        // calcul des cellules intersectantes
        for (E_Int now = 0; now < nwalls; now++)
        {
            indicesBB.clear(); candidates.clear();
            K_SEARCH::BbTree3D* bbtree = vectOfBBTrees[now];
            bbtree->getOverlappingBoxes(minB, maxB, indicesBB);
            FldArrayF* fieldv = fieldsw[now];
            posxw = posxv[now]; posyw = posyv[now]; poszw = poszv[now]; poscw = poscv[now];
            xw = fieldv->begin(posxw);
            yw = fieldv->begin(posyw);
            zw = fieldv->begin(poszw);
            cellnw = fieldv->begin(poscw);
            FldArrayI& cnloc = *cntw[now];
            nbb = indicesBB.size();
            nvert = cnloc.getNfld();
            E_Float prodCellN2 = pow(2.,nvert);
            for (E_Int i = 0; i < nbb; i++)
            {
              et = indicesBB[i];
              prod = 1.;
              for (E_Int novert = 1; novert <= nvert; novert++)
              {
                ind10 = cnloc(et, novert)-1;
                prod = prod*cellnw[ind10];
              }
              if (prod != 0. && prod != prodCellN2) candidates.push_back(et);
            }
            ret = K_COMPGEOM::projectOrthoPrecond(pt[0], pt[1], pt[2], xw, yw, zw,
                                                  candidates, cnloc, xp, yp, zp,
                                                  p0, p1, p2, p);
            if (ret != -1)
            {
              dx = xp-pt[0]; dy = yp-pt[1]; dz = zp-pt[2];
              dist = dx*dx + dy*dy + dz*dz;
              if (dist < distmin) { distancep[ind] = dist; distmin = dist; }
            }
          }
    }
#endif
