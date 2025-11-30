  void complementColor(float r, float g, float b, float& ro, float& go, float& bo);
  double resf = std::max(_view.w / 1080., 1.);

  switch (ptrState->meshStyle)
  {
    case 0:
      // Wires rouge
      color1[0] = 0.95; color1[1] = 0.; color1[2] = 0.;
      color2[0] = 0.1; color2[1] = 0.1; color2[2] = 1.;
      // Ecrasement si render tag
      if (zonep->meshColorR != -1.)
      { color1[0] = zonep->meshColorR; color1[1] = zonep->meshColorG; color1[2] = zonep->meshColorB; }
      else if (zonep->colorR != -1.)
      { complementColor(zonep->colorR, zonep->colorG, zonep->colorB, color1[0],color1[1],color1[2]); }
      if (zonep->meshWidth != -1.)
      { glLineWidth(zonep->meshWidth*resf); }
      else { glLineWidth(1.*resf); }
      break;

    case 1:
      // Wires multicolor rgb
      getrgb(this, zone*nz, &r, &g, &b);
      color1[0] = r; color1[1] = g; color1[2] = b;
      if (r > 0.8 && g < 0.2 && b < 0.2)
      { color1[0] = r; color1[1] = b; color1[2] = g; }
      color2[0] = 0.1; color2[1] = 0.1; color2[2] = 1.;
      // Ecrasement si render tag
      if (zonep->meshColorR != -1.)
      { color1[0] = zonep->meshColorR; color1[1] = zonep->meshColorG; color1[2] = zonep->meshColorB; }
      else if (zonep->colorR != -1.)
      { color1[0] = zonep->colorR; color1[1] = zonep->colorG; color1[2] = zonep->colorB; }
      if (zonep->meshWidth != -1.)
      { glLineWidth(zonep->meshWidth*resf); }
      else { glLineWidth(1.*resf); }
      break;

    case 2:
      // Wires multicolor gbr
      getrgb(this, zone*1./_numberOfZones, &r, &g, &b);
      color1[0] = g; color1[1] = b;  color1[2] = r;
      if (r > 0.8 && g < 0.2 && b < 0.2)
      { color1[0] = g; color1[1] = r; color1[2] = b; }
      //color2[0] = 0.4;  color2[1] = 0.4;  color2[2] = 1.;
      color2[0] = 0.7;  color2[1] = 0.88;  color2[2] = 1.;
      // Ecrasement si render tag
      if (zonep->meshColorR != -1.)
      { color1[0] = zonep->meshColorR; color1[1] = zonep->meshColorG; color1[2] = zonep->meshColorB; }
      else if (zonep->colorR != -1.)
      { complementColor(zonep->colorR,zonep->colorG,zonep->colorB, color1[0],color1[1],color1[2]); }  
      if (zonep->meshWidth != -1.)
      { glLineWidth(zonep->meshWidth*resf); }
      else { glLineWidth(1.*resf); }
      break;

   case 3:
      // Wires noirs
      color1[0] = 0.; color1[1] = 0.; color1[2] = 0.;
      color2[0] = 0.7; color2[1] = 0.88; color2[2] = 1.;
      // Ecrasement si render tag
      if (zonep->meshColorR != -1.)
      { color1[0] = zonep->meshColorR; color1[1] = zonep->meshColorG; color1[2] = zonep->meshColorB; }
      else if (zonep->colorR != -1.)
      { complementColor(zonep->colorR,zonep->colorG,zonep->colorB, color1[0],color1[1],color1[2]); }  
      if (zonep->meshWidth != -1.)
      { glLineWidth(zonep->meshWidth*resf); }
      else { glLineWidth(2.*resf); }
      break;

    case 4:
      // Wires noirs fins      
      color1[0] = 0.; color1[1] = 0.; color1[2] = 0.;
      color2[0] = 0.7; color2[1] = 0.88; color2[2] = 1.;
      // Ecrasement si render tag
      if (zonep->meshColorR != -1.)
      { color1[0] = zonep->meshColorR; color1[1] = zonep->meshColorG; color1[2] = zonep->meshColorB; }
      else if (zonep->colorR != -1.)
      { complementColor(zonep->colorR,zonep->colorG,zonep->colorB, color1[0],color1[1],color1[2]); }  
      if (zonep->meshWidth != -1.)
      { glLineWidth(zonep->meshWidth*resf); }
      else { glLineWidth(0.4*resf); }
      break;

    default:
      color1[0] = 0.95; color1[1] = 0.95; color1[2] = 1.;
      color2[0] = 0.1; color2[1] = 0.1; color2[2] = 1.;
      glLineWidth(1.*resf);
      break;
  }

  // 1D specific
  if (is1D) 
  { 
    if (zonep->colorR != -1.)
    { color1[0] = zonep->colorR; color1[1] = zonep->colorG; color1[2] = zonep->colorB; }
    color2[0] = 0.1; color2[1] = 0.1; color2[2] = 1.; 
    glLineWidth(3.);
  }
  