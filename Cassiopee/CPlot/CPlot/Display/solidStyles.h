  switch (ptrState->solidStyle)
  {
    case 0:
      // Couleur uniforme bleue
      color1[0] = 0.5; color1[1] = 0.5; color1[2] = 1.;
      color2[0] = 0.;  color2[1] = 0.;  color2[2] = 1;
      break;
      
    case 1:
      // Multicolor suivant le no de la zone
      getrgb(this, zone*nz, &r, &g, &b);
      color1[0] = r; color1[1] = g; color1[2] = b;
      if (b > 0.8 && r < 0.2 && g < 0.2) 
      { color1[0] = r; color1[1] = b; color1[2] = g; }
      color2[0] = 0.1; color2[1] = 0.1; color2[2] = 1;
      break;
      
    case 2:
      // Couleur uniforme blanche
      color1[0] = 1.; color1[1] = 1.; color1[2] = 1.;
      color2[0] = 1.; color2[1] = 1.; color2[2] = 1;
      break;

    case 3:
      // Couleur uniforme cyan
      color1[0] = 0.9; color1[1] = 1.; color1[2] = 1.;
      color2[0] = 0.1; color2[1] = 0.1; color2[2] = 1;
      break;

    case 4:
      // multicolor + sobel
      getrgb(this, zone*nz, &g, &r, &b);
      color1[0] = r; color1[1] = g; color1[2] = b;
      if (b > 0.8 && r < 0.2 && g < 0.2) 
      { color1[0] = r; color1[1] = b; color1[2] = g; }
      color2[0] = 0.1; color2[1] = 0.1; color2[2] = 1;
      break;

    default:
      // Couleur uniforme bleue
      color1[0] = 0.5; color1[1] = 0.5; color1[2] = 1.;
      color2[0] = 0.;  color2[1] = 0.;  color2[2] = 1;
      break;
  }
  // Ecrasement si renderTag
  if (zonep->colorR > -0.5)
  {color1[0] = zonep->colorR; color1[1] = zonep->colorG; color1[2] = zonep->colorB;}

  // specifique aux zones 1D
  if (is1D && ptrState->mode == RENDER) glLineWidth(1.+5*zonep->shaderParam1);
  else if (is1D) glLineWidth(3.);
  else glLineWidth(1.);
