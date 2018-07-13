pr1[0] = xc0; pr2[0]=xc0+dirx0; 
pr1[1] = yc0; pr2[1]=yc0+diry0; 
pr1[2] = zc0; pr2[2]=zc0+dirz0;
dist2 = K_CONST::E_MAX_FLOAT;
xsav = xc0; ysav = yc0; zsav = zc0; 
ok = -1;
for(E_Int nos = 0; nos < nfronts; nos++)
{
    K_SEARCH::BbTree3D* bbtree = vectOfFrontBBTrees[nos];
    bbtree->getIntersectingBoxes(pr1, pr2, indicesBB, tol);
    if ( indicesBB.size() != 0) 
    { 
        E_Float* xs = unstrfF[nos]->begin(posxf[nos]);
        E_Float* ys = unstrfF[nos]->begin(posyf[nos]);
        E_Float* zs = unstrfF[nos]->begin(poszf[nos]);
        notri = K_COMPGEOM::projectDir(xc0, yc0, zc0, dirx0, diry0, dirz0,
                                       xs, ys, zs, indicesBB, *cnf[nos], 
                                       xi0, yi0, zi0, oriented);                
        indicesBB.clear();

        if (notri > -1)
        {
            cnVert1 = cnf[nos]->begin(1);
            cnVert2 = cnf[nos]->begin(2);
            cnVert3 = cnf[nos]->begin(3);
            distl = (xc0-xi0)*(xc0-xi0)+(yc0-yi0)*(yc0-yi0)+(zc0-zi0)*(zc0-zi0);
            if (distl < dist2) 
            {
                dist2 = distl;
                xsav = xi0; ysav = yi0; zsav = zi0;            
                ok = 1;//one projection found

                // indvert1 = cnVert1[notri]-1;
                // indvert2 = cnVert2[notri]-1;
                // indvert3 = cnVert3[notri]-1;
                // edgeLen1 = (xs[indvert1]-xs[indvert2])*(xs[indvert1]-xs[indvert2])+
                // (ys[indvert1]-ys[indvert2])*(ys[indvert1]-ys[indvert2])+
                // (zs[indvert1]-zs[indvert2])*(zs[indvert1]-zs[indvert2]);

                // edgeLen2 = (xs[indvert3]-xs[indvert2])*(xs[indvert3]-xs[indvert2])+
                // (ys[indvert3]-ys[indvert2])*(ys[indvert3]-ys[indvert2])+
                // (zs[indvert3]-zs[indvert2])*(zs[indvert3]-zs[indvert2]);
                // snearloc0 = K_FUNC::E_max(edgeLen1, edgeLen2);
            }
        }         
    }
}
