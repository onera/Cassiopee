  double* x = zonep->x;
  double* y = zonep->y;
  double* z = zonep->z;
  int* connect = zonep->connect;
  // Grid dimensions
  int ne = zonep->ne;
  int ne2 = 2*ne;
  int ne3 = 3*ne;
  int ne4 = 4*ne;
  int ne5 = 5*ne;
  int ne6 = 6*ne;
  int ne7 = 7*ne;
  int eltType = zonep->eltType;
  int nd, l;
  int nbNodesPerPatch = zonep->eltSize;
  glPatchParameteri( GL_PATCH_VERTICES, 3 );
  glBegin(GL_PATCHES);
  //glBegin(GL_LINES);
  if (zonep->blank == 0)
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
        break;
        
      case 3: // QUAD
        if ( ( nbNodesPerPatch == 8 ) or ( nbNodesPerPatch == 9 ) )
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
      break;
    }
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
        if ( ( nbNodesPerPatch == 8 ) or ( nbNodesPerPatch == 9 ) )
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
