 if (isminortho == 1)
 {
   npts_walls_limit.push_back(nptsmax-1);
   // calcul de la connectivite vertex->elements
   // give a vertex & get the elements
   npts_local = fieldsw[v]->getSize();
   vector< vector<E_Int>  > cVE_local(npts_local);
   K_CONNECT::connectEV2VE(*cntw[v], cVE_local);
   cVE_all.push_back(cVE_local);
 }
