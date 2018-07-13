        //
        // Calcul dimension zone receveuse
        //
        PyObject* solR;
        if  (loc==0) { solR = K_PYTREE::getNodeFromName1(zoneR , "FlowSolution"        ); }
        else  { solR = K_PYTREE::getNodeFromName1(zoneR , "FlowSolution#Centers"); }
        t = K_PYTREE::getNodeFromName1(solR, varname);
        iptroR = K_PYTREE::getValueAF(t, hook);

        // get type
        t =  K_PYTREE::getNodeFromName1(zoneR, "ZoneType");
        type =  K_PYTREE::getValueS(t, s, hook);
        // get dims zone receveuse
        d  =  K_PYTREE::getValueAI(zoneR, s0, s1, hook);

        if  (K_STRING::cmp(type, s, "Structured") == 0)
         {
          E_Int shift = 0; if(loc == 1) shift = 3;
          if (s0 == 1) { ndimdxR= d[0+shift]; }
          else if (s0 == 2) { ndimdxR= d[0+shift]*d[1+shift]; } 
          else if (s0 == 3) { ndimdxR= d[0+shift]*d[1+shift]*d[2+shift]; } 
         }
        else // non structure
         {
             ndimdxR= d[0]* d[1]; // npoint, nelements
         }
