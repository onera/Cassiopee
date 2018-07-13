pr1[0]=xc0; pr2[0]=xc0+dirx0; 
pr1[1]=yc0; pr2[1]=yc0+diry0; 
pr1[2]=zc0; pr2[2]=zc0+dirz0;
dist2 = K_CONST::E_MAX_FLOAT;
xsav = xc0; ysav = yc0; zsav = zc0; 
ok = -1;

for(E_Int nos = 0; nos < nbodies; nos++)
{
    K_SEARCH::BbTree3D* bbtree = vectOfBodyBBTrees[nos];
    bbtree->getIntersectingBoxes(pr1, pr2, indicesBB, tol);
    if ( indicesBB.size() != 0 ) 
    {
        E_Float* xs = unstrbF[nos]->begin(posxb[nos]);
        E_Float* ys = unstrbF[nos]->begin(posyb[nos]);
        E_Float* zs = unstrbF[nos]->begin(poszb[nos]);
        notri = K_COMPGEOM::projectDir(xc0, yc0, zc0, dirx0, diry0, dirz0,
                                       xs, ys, zs, indicesBB, *cnb[nos], 
                                       xi0, yi0, zi0, oriented);                
        indicesBB.clear();
        if ( notri > -1)
        {
            distl = (xc0-xi0)*(xc0-xi0)+(yc0-yi0)*(yc0-yi0)+(zc0-zi0)*(zc0-zi0);
            if (distl < dist2) 
            {
                dist2 = distl;
                xsav = xi0; ysav = yi0; zsav = zi0;            
                ok = 1;//one projection found
                noibctype = nos;
            }
        }         
    }
}
