// calcul de la bounding box de la sphere de rayon PP'
rx = pt[0]-xf2[indp]; ry = pt[1]-yf2[indp]; rz = pt[2]-zf2[indp];
rad = sqrt(rx*rx+ry*ry+rz*rz);
minB[0] = pt[0]-rad; minB[1] = pt[1]-rad; minB[2] = pt[2]-rad;
maxB[0] = pt[0]+rad; maxB[1] = pt[1]+rad; maxB[2] = pt[2]+rad;
dist2 = K_CONST::E_MAX_FLOAT;
for(E_Int nos = 0; nos < nfronts; nos++)
{
    K_SEARCH::BbTree3D* bbtree = vectOfFrontBBTrees[nos];
    bbtree->getOverlappingBoxes(minB, maxB, indicesBB);
    if ( indicesBB.size() != 0) 
    {
        E_Float* xs = unstrfF[nos]->begin(posxf[nos]);
        E_Float* ys = unstrfF[nos]->begin(posyf[nos]);
        E_Float* zs = unstrfF[nos]->begin(poszf[nos]);
        notri = K_COMPGEOM::projectOrthoPrecond(xc0, yc0, zc0, xs, ys, zs, 
                                                indicesBB, *cnf[nos], 
                                                xw0, yw0, zw0);
        indicesBB.clear();
        if ( notri > -1)
        {
            cnVert1 = cnf[nos]->begin(1);
            cnVert2 = cnf[nos]->begin(2);
            cnVert3 = cnf[nos]->begin(3);
            distl = (xc0-xw0)*(xc0-xw0)+(yc0-yw0)*(yc0-yw0)+(zc0-zw0)*(zc0-zw0);
            if (distl < dist2) 
            {
                dist2 = distl;
                xsav = xw0; ysav = yw0; zsav = zw0;            
                ok = 1;//one projection found

                indvert1 = cnVert1[notri]-1;
                indvert2 = cnVert2[notri]-1;
                indvert3 = cnVert3[notri]-1;
                edgeLen1 = (xs[indvert1]-xs[indvert2])*(xs[indvert1]-xs[indvert2])+
                (ys[indvert1]-ys[indvert2])*(ys[indvert1]-ys[indvert2])+
                (zs[indvert1]-zs[indvert2])*(zs[indvert1]-zs[indvert2]);

                edgeLen2 = (xs[indvert3]-xs[indvert2])*(xs[indvert3]-xs[indvert2])+
                (ys[indvert3]-ys[indvert2])*(ys[indvert3]-ys[indvert2])+
                (zs[indvert3]-zs[indvert2])*(zs[indvert3]-zs[indvert2]);
                snearloc = K_FUNC::E_max(edgeLen1, edgeLen2);
            }
        } 
    }
}
