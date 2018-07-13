// calcul de la bounding box de la sphere de rayon PP'
rx = pt[0]-xb2[indp]; ry = pt[1]-yb2[indp]; rz = pt[2]-zb2[indp];
rad = sqrt(rx*rx+ry*ry+rz*rz);
minB[0] = pt[0]-rad; minB[1] = pt[1]-rad; minB[2] = pt[2]-rad;
maxB[0] = pt[0]+rad; maxB[1] = pt[1]+rad; maxB[2] = pt[2]+rad;
dist2 = K_CONST::E_MAX_FLOAT;
for(E_Int nos = 0; nos < nbodies; nos++)
{
    K_SEARCH::BbTree3D* bbtree = vectOfBodyBBTrees[nos];
    bbtree->getOverlappingBoxes(minB, maxB, indicesBB);
    if ( indicesBB.size() != 0 ) 
    {
        E_Float* xs = unstrbF[nos]->begin(posxb[nos]);
        E_Float* ys = unstrbF[nos]->begin(posyb[nos]);
        E_Float* zs = unstrbF[nos]->begin(poszb[nos]);
        notri = K_COMPGEOM::projectOrthoPrecond(xc0, yc0, zc0, xs, ys, zs, 
                                                indicesBB, *cnb[nos], 
                                                xw0, yw0, zw0);
        indicesBB.clear();
        if ( notri > -1)
        {
            distl = (xc0-xw0)*(xc0-xw0)+(yc0-yw0)*(yc0-yw0)+(zc0-zw0)*(zc0-zw0);
            if (distl < dist2) 
            {
                dist2 = distl;
                xsav = xw0; ysav = yw0; zsav = zw0;            
                ok = 1;//one projection found
                noibctype = nos;
            }
        } 
    }
}
