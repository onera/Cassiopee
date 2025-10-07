      for (vw = 0; vw < nwalls; vw++)
      {
        if (indw2 <= npts_walls_limit[vw]) break;
      }
      if (vw > 0) indw2 = indw2-(npts_walls_limit[vw-1]+1);

      FldArrayF* fieldv_local = fieldsw[vw];
      E_Int posxw_local = posxv[vw];
      E_Int posyw_local = posyv[vw];
      E_Int poszw_local = poszv[vw];
      E_Float* xw_local = fieldv_local->begin(posxw_local);
      E_Float* yw_local = fieldv_local->begin(posyw_local);
      E_Float* zw_local = fieldv_local->begin(poszw_local);
      E_Float xp_local, yp_local, zp_local;
      FldArrayI& cnloc_local = *cntw[vw];
      vector< vector<E_Int> >& cVE_11_local = cVE_all[vw];
      vector<E_Int> & cVE_1_local = cVE_11_local[indw2];

      ret = K_COMPGEOM::projectOrthoPrecond(pt[0], pt[1], pt[2], xw_local, yw_local, zw_local,
                          cVE_1_local, cnloc_local, xp_local, yp_local, zp_local,
                          p0, p1, p2, p);
