pr1[0] = xc0; pr2[0]=xc0+dirx0; 
pr1[1] = yc0; pr2[1]=yc0+diry0; 
pr1[2] = zc0; pr2[2]=zc0+dirz0;
dist2 = K_CONST::E_MAX_FLOAT;
xsf2 = xc0; ysf2 = yc0; zsf2 = zc0; 
ok = -1;
for (E_Int nos = 0; nos < nfront2s; nos++)
{
    K_SEARCH::BbTree3D* bbtree = vectOfFront2BBTrees[nos];
    bbtree->getIntersectingBoxes(pr1, pr2, indicesBB, tol);
    if (indicesBB.size() != 0) 
    { 
        E_Float* xs = unstrf2F[nos]->begin(posxf2[nos]);
        E_Float* ys = unstrf2F[nos]->begin(posyf2[nos]);
        E_Float* zs = unstrf2F[nos]->begin(poszf2[nos]);
        notri = K_COMPGEOM::projectDir(xc0, yc0, zc0, dirx0, diry0, dirz0,
                                       xs, ys, zs, indicesBB, *cnf2[nos], 
                                       xi0, yi0, zi0, oriented);                
        indicesBB.clear();

        if (notri > -1)
        {
            cnVert1 = cnf2[nos]->begin(1);
            cnVert2 = cnf2[nos]->begin(2);
            cnVert3 = cnf2[nos]->begin(3);
            distl = (xc0-xi0)*(xc0-xi0)+(yc0-yi0)*(yc0-yi0)+(zc0-zi0)*(zc0-zi0);
            if (distl < dist2) 
            {
                dist2 = distl;
                xsf2 = xi0; ysf2 = yi0; zsf2 = zi0;            
                ok = 1;//one projection found
            }
        }         
    }
}
