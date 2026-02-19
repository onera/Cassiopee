/*    
    Copyright 2013-2026 ONERA.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/

// Formatted ensight routines - by R. Barrier (not usable)

# include <stdio.h>

//===========================================================================
// ecriture_fmt_ensight : ecrit un champ structure au format ensight5 
// entrees: 
//    - nom_fichier : le nom generique du fichier de sortie (il sera suivi de .geo pour le maillage, .res pour le champ a�ro,
//		    .scl pour le nom des variables aero)
//    - titre : un titre dans l'entete du fichier
//    - nom_variable : le tableau du nom des variables
//    - endi, endj, endk: dimensions du champ a�ro
//    - nombre_variable: nombre de variables dans le champ a�ro initial (sans les variables compos�es)
//    - var: champ aero range dans un tableau 1D par i puis j puis k puis numero_variable croissants
// ===========================================================================
void ecriture_fmt_ensight(string nom_fichier, string titre, 
                          string* nom_variable, long endi, 
                          long endj, long endk, long nombre_variable, 
                          double* var)
{
   // maillage
   nom_fichier=nom_fichier+".geo";
   string entete1, entete2;
   entete1="Fichier geometrie ensight5, export de Zeppelin v1.2";   
   entete2="TITRE : "+titre;
   
   long num_elem=1;
    
   cout << "   ecriture -- " << nom_fichier << endl;    
   FILE* pFile = fopen (nom_fichier.c_str(),"wt");   
   if (pFile!=NULL)
   {	  
      fprintf (pFile , "%s\n", entete1.c_str());      
      fprintf (pFile , "%s\n", entete2.c_str());
      fprintf (pFile , "node id given\n");
      fprintf (pFile , "element id given\n");
      fprintf (pFile , "coordinates\n");
      fprintf (pFile , "%8d\n", endi*endj*endk);
      
      if (nombre_variable>=3)
      {
         for (long i=0;i<endi*endj*endk;i=i+1)
         {
   	    fprintf (pFile , "%8d%12.5e%12.5e%12.5e\n", i+1, var[i], var[i+endi*endj*endk], var[i+2*endi*endj*endk]);   	 
         }  
      }
      if (nombre_variable==2)
      {
         for (long i=0;i<endi*endj*endk;i=i+1)
         {
   	    fprintf (pFile , "%8d%12.5e%12.5e%12.5e\n", i+1, (double)i, var[i], var[i+endi*endj*endk]);   	 
         }  
      }
      if (nombre_variable==1)
      {
         for (long i=0;i<endi*endj*endk;i=i+1)
         {
   	    fprintf (pFile , "%8d%12.5e%12.5e%12.5e\n", i+1, (double)i, var[i], 0.0);   	 
         }  
      }
      
      fprintf (pFile , "part 1\n");
      fprintf (pFile , "description line\n");
      // hexa8
      if (endi>1 && endj>1 && endk>1)
      { 
         fprintf (pFile , "hexa8\n");
         fprintf (pFile , "%8d\n", (endi-1)*(endj-1)*(endk-1));
	 for (long k=0;k<(endk-1);k=k+1)
         {
   	    for (long j=0;j<(endj-1);j=j+1)
            {
	       for (long i=0;i<(endi-1);i=i+1)
               {
	          fprintf (pFile , "%8d%8d%8d%8d%8d%8d%8d%8d%8d\n", num_elem,
		  						    i+1+(j*endi)+(k*endi*endj), i+2+(j*endi)+(k*endi*endj),
								    i+2+((j+1)*endi)+(k*endi*endj), i+1+((j+1)*endi)+(k*endi*endj),
								    i+1+(j*endi)+((k+1)*endi*endj), i+2+(j*endi)+((k+1)*endi*endj),
								    i+2+((j+1)*endi)+((k+1)*endi*endj), i+1+((j+1)*endi)+((k+1)*endi*endj));	
		  num_elem=num_elem+1;
	       }
	    }						     	 
         }         
         fclose (pFile);
      }
      // quad4      
      else if (endi==1 && endj>1 && endk>1)
      {
         fprintf (pFile , "quad4\n");
         fprintf (pFile , "%8d\n", (endj-1)*(endk-1));	 
	 for (long j=0;j<(endk-1);j=j+1)
         { 
   	    for (long i=0;i<(endj-1);i=i+1)
	    {
	       fprintf (pFile , "%8d%8d%8d%8d%8d\n", num_elem, i+1+(j*endj), i+2+(j*endj), i+2+endj+(j*endj), i+1+endj+(j*endj));	   
	       num_elem=num_elem+1;   
	    }
         }  
      }
      else if (endj==1 && endi>1 && endk>1)
      {   
	 fprintf (pFile , "quad4\n");
         fprintf (pFile , "%8d\n", (endi-1)*(endk-1));
	 for (long j=0;j<(endk-1);j=j+1)
         { 
   	    for (long i=0;i<(endi-1);i=i+1)
	    {
	       fprintf (pFile , "%8d%8d%8d%8d%8d\n", num_elem, i+1+(j*endi), i+2+(j*endi), i+2+endi+(j*endi), i+1+endi+(j*endi));
	       num_elem=num_elem+1; 	      
	    }
         } 
      }
      else if (endk==1 && endi>1 && endj>1)
      {      
	 fprintf (pFile , "quad4\n");
         fprintf (pFile , "%8d\n", (endi-1)*(endj-1));
	 for (long j=0;j<(endj-1);j=j+1)
         { 
   	    for (long i=0;i<(endi-1);i=i+1)
	    {
	       fprintf (pFile , "%8d%8d%8d%8d%8d\n", num_elem, i+1+(j*endi), i+2+(j*endi), i+2+endi+(j*endi), i+1+endi+(j*endi));	
	       num_elem=num_elem+1;       
	    }
         } 
      }
      // bar2
      else if (endi==1 && endj==1 && endk>1)
      { 
         fprintf (pFile , "bar2\n");
         fprintf (pFile , "%8d\n", (endk-1));
	 for (long i=0;i<(endk-1);i=i+1)
         { 
   	    fprintf (pFile , "%8d%8d%8d\n", num_elem, i+1, i+2);      
	    num_elem=num_elem+1;     
         } 
      }
      else if (endj==1 && endk==1 && endi>1)
      { 
         fprintf (pFile , "bar2\n");
         fprintf (pFile , "%8d\n", (endi-1));
	 for (long i=0;i<(endi-1);i=i+1)
         { 
   	    fprintf (pFile , "%8d%8d%8d\n", num_elem, i+1, i+2);		   
	    num_elem=num_elem+1;
         } 
      } 
      else if (endk==1 && endi==1 && endj>1)
      { 
         fprintf (pFile , "bar2\n");
         fprintf (pFile , "%8d\n", (endj-1));
	 for (long i=0;i<(endj-1);i=i+1)
         { 
   	    fprintf (pFile , "%8d%8d%8d\n", num_elem, i+1, i+2);	
	    num_elem=num_elem+1;     
         } 
      }   
      //point
      else if (endi==1 && endj==1 && endk==1)
      { 
         fprintf (pFile , "point\n");
         fprintf (pFile , "%8d\n", 1);
	 fprintf (pFile , "%8d\n", 1);
      }   
      
   }
   else
   {
      cout << "ecriture_fmt_ensight : echec lors de l'ouverture de " << nom_fichier << endl; 
   }
   
   
   if (nombre_variable>3)
   {
      // resultat
      string nom_resultat;
      nom_fichier=nom_fichier.substr(0,nom_fichier.length()-3);
      nom_fichier=nom_fichier+"res";  
      cout << "   ecriture -- " << nom_fichier << endl; 
      FILE* pFile = fopen (nom_fichier.c_str(),"wt");   
      if (pFile!=NULL)
      {	  
         fprintf (pFile , SF_3D_ "\n", nombre_variable-3, 0, 0);
         fprintf (pFile , "1\n");
         fprintf (pFile , "0.0\n");
         for (long nvar=0;nvar<nombre_variable-3;nvar=nvar+1)
         {
   	    nom_resultat=nom_fichier.substr(0,nom_fichier.length()-3);
            nom_resultat=nom_resultat+"scl"+entierAchaine(nvar);  
	    fprintf (pFile , "%s %s\n", nom_resultat.c_str(), nom_variable[3+nvar].c_str());   	 
         } 
         fclose (pFile);
      }
      else
      {
         cout << "ecriture_fmt_ensight : echec lors de l'ouverture de " << nom_fichier << endl; 
      }
      
      
      // scalaires
      for (long nvar=0;nvar<nombre_variable-3;nvar=nvar+1)
      {
         nom_resultat=nom_fichier.substr(0,nom_fichier.length()-3);
         nom_resultat=nom_resultat+"scl"+entierAchaine(nvar);  
	 cout << "   ecriture -- " << nom_resultat << endl; 
         FILE* pFile = fopen (nom_resultat.c_str(),"wt");   
         if (pFile!=NULL)
         {	  
            fprintf (pFile , "Fichier scalaire ensight5, variable : %s \n", nom_variable[3+nvar].c_str());
	    
	    for (long i=0;i<(endi*endj*endk);i=i+1)
            {
   	       fprintf (pFile , "%12.5e", var[i+(endi*endj*endk)*(nvar+3)]);
	       if ((i+1)%6==0)
   	       {
   	          fprintf (pFile , "\n");
   	       }	 
            } 
            fclose (pFile);
         }
         else
         {
            cout << "ecriture_fmt_ensight : echec lors de l'ouverture de " << nom_resultat << endl; 
         }
      }
   }
   
}
