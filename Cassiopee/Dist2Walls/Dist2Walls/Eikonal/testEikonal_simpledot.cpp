# include <cstdlib>
# include <string>
# include <cstdio>
# include <cmath>
# include <ctime>
# include <chrono>
# include <iostream>
# include "Eikonal.hpp"

bool isSource( double x, double y, double z )
{
  double x1 = x+0.5, y1 = y, z1 = z;
  double x2 = x-0.5, y2 = y, z2 = z;
  if ( ( x1*x1+y1*y1+z1*z1 < 1.E-2) || ( std::fabs(x2*x2+y2+z2) < 1.E-2 ) ) return true;
  return false;
}

double speedFunction( double x, double y, double z )
{
  if ( x+y < 0.5 ) return 0.75;
  return 1.;
}

void saveVTK( const std::string& fileName, unsigned n, double lbx, double lby, double lbz,
	      double h, const double* sol )
{
  FILE* fich = fopen(fileName.c_str(), "w");
  fprintf(fich, "# vtk DataFile Version 2.0\n");
  fprintf(fich, "Distance des noeuds a l'origine\n");
  fprintf(fich, "ASCII\n");
  fprintf(fich, "DATASET STRUCTURED_GRID\n");
  fprintf(fich, "DIMENSIONS %d %d %d\n",n,n,n);
  fprintf(fich, "POINTS %d double\n", n*n*n);
  double cz = 0.;
  for ( unsigned k = 0; k < n; ++k ) {
    double cy = 0.;
    for ( unsigned j = 0; j < n; ++j ) {
      double cx = 0.;
      for ( unsigned i = 0; i < n; ++i ) {
	fprintf(fich, "%lg %lg %lg\n", cx, cy, cz );
	cx += h;
      }
      cy += h;
    }
    cz += h;
  }
  fprintf(fich, "POINT_DATA %d\n",n*n*n);
  fprintf(fich, "SCALARS Distance double\n");
  fprintf(fich, "LOOKUP_TABLE default\n");

  for ( unsigned i = 0; i < n*n*n; ++i )
    fprintf(fich, "%lg\t", h*sol[i]);
  fprintf(fich, "\n");
  fclose(fich);
}

int main()
{
  unsigned n = 512;// Nombre de points par direction
  double lbx = -1, lby = -1, lbz = -1;
  double h = 2./n;
  double* sol = new double[n*n*n];
  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::cout << "Commence algo par point sequentiel\n";
  start = std::chrono::system_clock::now();
  solveEikonalOnIsotropGrid( n, n, n, lbx, lby, lbz, h, isSource, speedFunction, sol );
  end = std::chrono::system_clock::now();
  std::cout << "Fini algo par point sequentiel\n";

  std::chrono::duration<double> elapsed_seconds = end-start;
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);

  std::cout << "finished sequential Eikonal in " << std::ctime(&end_time)
	    << "elapsed time: " << elapsed_seconds.count() << "s\n";
  // On sauvegarde sol dans un gnuplot :
  //saveVTK( "Distance_simpledot.vtk", n, lbx, lby, lbz, h, sol );
  std::cout << "Commence algo par bloc\n";
  start = std::chrono::system_clock::now();
  blockFIM( n, n, n, lbx, lby, lbz, h, 8, 8, 8, 5, isSource,
	    speedFunction, sol );
  end = std::chrono::system_clock::now();
  std::cout << "Fini algo par bloc\n";
  elapsed_seconds = end-start;
  end_time = std::chrono::system_clock::to_time_t(end);

  std::cout << "finished Block Parallel Eikonal in " << std::ctime(&end_time)
	    << "elapsed time: " << elapsed_seconds.count() << "s\n";
  //saveVTK( "Distance_blockdot.vtk", n, lbx, lby, lbz, h, sol );
  delete [] sol;
  return EXIT_SUCCESS;
}
