// uscaln: u au point cible scalaire la normale au pt cible
uscaln = u*n0 + v*n1 + w*n2;
  
//composante normale de la vitesse
un = uscaln*n0;
vn = uscaln*n1;
wn = uscaln*n2;

//composante tangentielle de la vitesse au pt cible
ut = u-un;
vt = v-vn;
wt = w-wn;

// uext: norme de la composante tangentielle de la vitesse externe
uext = sqrt(ut*ut+vt*vt+wt*wt);
uext = std::max(uext, 1.e-12);
