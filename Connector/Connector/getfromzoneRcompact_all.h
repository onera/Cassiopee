       //
       // Calcul dimension zone receveuse
       //
       PyObject* sol;  PyObject* t;

       sol = K_PYTREE::getNodeFromName1(zoneR , "FlowSolution#Centers");
       if(sol != NULL)
       {  
         ipt_roR[nd]   = NULL;
         t                           = K_PYTREE::getNodeFromName1(sol, varname );
         if(t != NULL) ipt_roR[nd]   = K_PYTREE::getValueAF(t, hook);
       }

       sol = K_PYTREE::getNodeFromName1(zoneR , "FlowSolution");
       if(sol != NULL)
       {  
         ipt_roR_vert[nd] = NULL;
         t                              = K_PYTREE::getNodeFromName1(sol, varname );
         if(t != NULL) ipt_roR_vert[nd] = K_PYTREE::getValueAF(t, hook);
       }


        //if( loc==0) { solR = K_PYTREE::getNodeFromName1(zoneR , "FlowSolution"        ); }
        //else        { solR = K_PYTREE::getNodeFromName1(zoneR , "FlowSolution#Centers"); }
        //t             = K_PYTREE::getNodeFromName1(solR, varname );
        //iptroR        = K_PYTREE::getValueAF(t, hook);

        char* type; E_Int s, s0, s1;  E_Int* d;
        // get type
        t    =  K_PYTREE::getNodeFromName1(zoneR, "ZoneType");
        type =  K_PYTREE::getValueS(t, s, hook);
        // get dims zone receveuse
        d    =  K_PYTREE::getValueAI(zoneR, s0, s1, hook);

        if (K_STRING::cmp(type, s, "Structured") == 0)
         {

          if      (s0 == 1){ ipt_ndimdxR[nd + nidomR]= d[0]          ; ipt_ndimdxR[nd]= d[3];           }
          else if (s0 == 2){ ipt_ndimdxR[nd + nidomR]= d[0]*d[1]     ; ipt_ndimdxR[nd]= d[3]*d[4];      }
          else if (s0 == 3){ ipt_ndimdxR[nd + nidomR]= d[0]*d[1]*d[2]; ipt_ndimdxR[nd]= d[3]*d[4]*d[5]; }
          
          /*E_Int shift = 0; if(loc == 1) shift =1;
          if      (s0 == 1) { ndimdxR= d[0+shift];          }
          else if (s0 == 2) { ndimdxR= d[0+shift]*d[1+shift];     } 
          else if (s0 == 3) { ndimdxR= d[0+shift]*d[1+shift]*d[2+shift];} */
         }

        else // non structure
         {
              ipt_ndimdxR[nd        ]= d[0]* d[1]; // npoint, nelements
              ipt_ndimdxR[nd+ nidomR]= d[0]* d[1]; // pourquoi pas de dependance a la localisation en non structure??
         }
