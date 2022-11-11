  double* x = zonep->x;
  double* y = zonep->y;
  double* z = zonep->z;
  E_Int* connect = zonep->connect;
  // Grid dimensions
  E_Int ne = zonep->ne;
  E_Int ne2 = 2*ne;
  E_Int ne3 = 3*ne;
  E_Int ne4 = 4*ne;
  E_Int ne5 = 5*ne;
  E_Int ne6 = 6*ne;
  E_Int ne7 = 7*ne;
  E_Int ne8 = 8*ne;
  E_Int ne9 = 9*ne;
  E_Int ne10= 10*ne;
  E_Int ne11= 11*ne;
  E_Int eltType = zonep->eltType;
  E_Int nbNodesPerPatch = zonep->eltSize;
  glPatchParameteri( GL_PATCH_VERTICES, 3 );
  glBegin(GL_PATCHES);
  //glBegin(GL_LINES);
  if (zonep->blank == 0)
  {
    switch (eltType)
    {
      case 2: // TRI
        if (nbNodesPerPatch==6)
        {
          for (i = 0; i < ne; i++)
          {
            n1 = connect[i]-1;
            n2 = connect[i+ne ]-1;
            n3 = connect[i+ne3]-1;
            PLOTHO;
            
            n1 = connect[i+ne]-1;
            n2 = connect[i+ne2]-1;
            n3 = connect[i+ne4]-1;
            PLOTHO;
          
            n1 = connect[i+ne2]-1;
            n2 = connect[i]-1;
            n3 = connect[i+ne5]-1;
            PLOTHO;
          }
        }
        else if ( (nbNodesPerPatch==9) or (nbNodesPerPatch==10) )
        {
          for (i = 0; i < ne; i++)
          {
            n1 = connect[i]-1;
            n2 = connect[i+ne3 ]-1;
            n3 = connect[i+ne4]-1;
            PLOTHO;
            
            n1 = connect[i+ne3]-1;
            n2 = connect[i+ne ]-1;
            n3 = connect[i+ne4]-1;
            PLOTHO;
 
            n1 = connect[i+ne ]-1;
            n2 = connect[i+ne5]-1;
            n3 = connect[i+ne6]-1;
            PLOTHO;
 
            n1 = connect[i+ne5]-1;
            n2 = connect[i+ne2]-1;
            n3 = connect[i+ne6]-1;
            PLOTHO;

            n1 = connect[i+ne2]-1;
            n2 = connect[i+ne7]-1;
            n3 = connect[i+ne8]-1;
            PLOTHO;

            n1 = connect[i+ne7]-1;
            n2 = connect[i]-1;
            n3 = connect[i+ne8]-1;
            PLOTHO;
          }
        }
        break;
        
      case 3: // QUAD
        if ( (nbNodesPerPatch == 8) or (nbNodesPerPatch == 9) )
        {
          for (i = 0; i < ne; i++)
          {
            n1 = connect[i]-1;
            n2 = connect[i+ne ]-1;
            n3 = connect[i+ne4]-1;
            PLOTHO;
          
            n1 = connect[i+ne ]-1;
            n2 = connect[i+ne2]-1;
            n3 = connect[i+ne5]-1;
            PLOTHO;
          
            n1 = connect[i+ne2]-1;
            n2 = connect[i+ne3]-1;
            n3 = connect[i+ne6]-1;
            PLOTHO;
          
            n1 = connect[i+ne3]-1;
            n2 = connect[i]-1;
            n3 = connect[i+ne7]-1;
            PLOTHO;
          }// for
        }// if ( nbNodesPerPatch)
        else if ( (nbNodesPerPatch == 12) or (nbNodesPerPatch == 16) )
        {
          for (i = 0; i < ne; i++)
          {
            n1 = connect[i]-1;
            n2 = connect[i+ne5]-1;
            n3 = connect[i+ne4]-1;
            PLOTHO;

            n1 = connect[i+ne4]-1;
            n2 = connect[i+ne ]-1;
            n3 = connect[i+ne5]-1;
            PLOTHO;

            n1 = connect[i+ne ]-1;
            n2 = connect[i+ne7]-1;
            n3 = connect[i+ne6]-1;
            PLOTHO;

            n1 = connect[i+ne6]-1;
            n2 = connect[i+ne2]-1;
            n3 = connect[i+ne7]-1;
            PLOTHO;

            n1 = connect[i+ne2]-1;
            n2 = connect[i+ne9]-1;
            n3 = connect[i+ne8]-1;
            PLOTHO;

            n1 = connect[i+ne8]-1;
            n2 = connect[i+ne3]-1;
            n3 = connect[i+ne9]-1;
            PLOTHO;

            n1 = connect[i+ne3 ]-1;
            n2 = connect[i+ne11]-1;
            n3 = connect[i+ne10]-1;
            PLOTHO;

            n1 = connect[i+ne10]-1;
            n2 = connect[i     ]-1;
            n3 = connect[i+ne11]-1;
            PLOTHO;
          }// for

        }// End if
        break;
    }// End switch
  }
  else // With blanking
  {
    switch (eltType)
    {
        
      case 2: // TRI
        assert(nbNodesPerPatch == 6);
        for (i = 0; i < ne; i++)
        {
          n1 = connect[i]-1;
          n2 = connect[i+ne ]-1;
          n3 = connect[i+ne3]-1;
          PLOTBHO;
          
          n1 = connect[i+ne]-1;
          n2 = connect[i+ne2]-1;
          n3 = connect[i+ne5]-1;
          PLOTBHO;
          
          n1 = connect[i+ne2]-1;
          n2 = connect[i]-1;
          n3 = connect[i+ne4]-1;
          PLOTBHO;
        }
        break;
        
      case 3: // QUAD
        if ( (nbNodesPerPatch == 8) or (nbNodesPerPatch == 9) )
        {
          for (i = 0; i < ne; i++)
          {
            n1 = connect[i]-1;
            n2 = connect[i+ne ]-1;
            n3 = connect[i+ne4]-1;
            PLOTBHO;
          
            n1 = connect[i+ne ]-1;
            n2 = connect[i+ne2]-1;
            n3 = connect[i+ne5]-1;
            PLOTBHO;
          
            n1 = connect[i+ne2]-1;
            n2 = connect[i+ne3]-1;
            n3 = connect[i+ne6]-1;
            PLOTBHO;
          
            n1 = connect[i+ne3]-1;
            n2 = connect[i]-1;
            n3 = connect[i+ne7]-1;
            PLOTBHO;
          }
        }
    }// switch eltType
  }// else blanking
  glEnd();
