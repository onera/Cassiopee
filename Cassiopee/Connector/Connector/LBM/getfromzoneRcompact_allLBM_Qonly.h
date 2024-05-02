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

	 ipt_distR[nd]  = NULL;
	 t              = K_PYTREE::getNodeFromName1(sol,"TurbulentDistance" );
	 if(t != NULL) ipt_distR[nd]= K_PYTREE::getValueAF(t, hook);
 
       }

       PyObject* own      = K_PYTREE::getNodeFromName1(zoneR , ".Solver#ownData");
               t          = K_PYTREE::getNodeFromName1(own, "Parameter_int");
       ipt_param_intR[nd] = K_PYTREE::getValueAI(t, hook);

               t          = K_PYTREE::getNodeFromName1(own, "Parameter_real");
       ipt_param_realR[nd]= K_PYTREE::getValueAF(t, hook);


