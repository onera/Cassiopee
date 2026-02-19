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
//Authors : Sam Landier (sam.landier@onera.fr)

#include "Nuga/include/Plaster.h"
#include "Nuga/include/Triangle.h"
#include "Nuga/include/EltAlgo.h"
#include "Nuga/include/MeshTool.h"
#include "Nuga/include/FittingBox.h"
#include "Nuga/include/T3Mesher.h"
#include "Nuga/include/KdTree.h"
#include "Nuga/include/BbTree.h"
#ifdef WIN32
#ifdef E_DEBUG
#include "meshIO/meshIO.h"
#endif
#endif

using namespace NUGA;

#ifdef WIN32
#ifdef E_DEBUG
int Plaster::_count = 0;
#endif
#endif

using namespace NUGA;


Plaster::Plaster(void)
{
}

Plaster::~Plaster(void)
{
}

///
E_Int 
Plaster::make
(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect,
 K_FLD::FloatArray& plaster, E_Int & ni, E_Float bump_factor)
{
  K_FLD::FloatArray       P(3,3), iP(3,3);// Transform matrix from fitting space to real space and inverse.
  K_FLD::FloatArray       posE2(pos), pos2D, plaster2D;
  K_FLD::IntArray         connectE2(connect);
  int_vector_type         nodesE2, /*hN(hard_nodes),*/ new_IDs;
  //size_t                  nb_nodes;
  E_Int                   nj, err, NIJMAX(1000);
  E_Float                 minB[3], maxB[3], z0;

  // hack
  bool refine = true;
  if (ni ==-1) refine = false;

  ni = 0;
  plaster.clear();

  if (pos.cols() == 0 || connect.cols() < 3 )
  {
    ni = 1; // to prevent floating point exception (dividing by it) in fittinPlaster.cpp)
    return 0;
  }

  bump_factor = K_FUNC::E_max(bump_factor, -1.); // factor must be in [-1., 1.]
  bump_factor = K_FUNC::E_min(bump_factor, 1.);

  // Work only on connect points.
  NUGA::MeshTool::compact_to_mesh(posE2, connectE2, new_IDs);
  //update hN
//  for (size_t i = 0; i < hN.size(); ++i)
//    hN[i] = new_IDs[hN[i]];

  connectE2.uniqueVals(nodesE2);
  //nodesE2.insert(nodesE2.end(), hN.begin(), hN.end());
  //nb_nodes = nodesE2.size();

  /* Computes the fitting box coordinate system optimizing the view over the contour*/
  // cela fait appel au mailleur ??
  err = FittingBox::computeOptimalViewFrame(posE2, connectE2, P);
  if (err) return err;
  
  iP = P;
  K_FLD::FloatArray::inverse3(iP);
  NUGA::transform(posE2, iP);// Now we are in the fitting coordinate system.

  NUGA::MeshTool::boundingBox(posE2, connectE2, minB, maxB);

#ifdef WIN32
#ifdef E_DEBUG
  drawBox(posE2, connectE2, minB, maxB);
  meshIO::write("transformed.mesh", posE2, connectE2);
#endif
#endif

  // Project along Z(0,0,1) the contour nodes on the box' top side (z = maxB[2])
  // i.e where the point of view of the contour is maximum. 
  z0 = maxB[2];
  std::vector<E_Float> zE2(*std::max_element(nodesE2.begin(), nodesE2.end())+1, 0.);
  pos2D = posE2;
  pos2D.resize(2, pos2D.cols());
  E_Float devmin(NUGA::FLOAT_MAX);
  E_Float devmax(0.);        // pos2D is the projection on any plane normal to z.
  for (E_Int i = 0; i < posE2.cols(); ++i){ // Set the z for the contour nodes.
    zE2[i] = z0 - posE2(2,i);
    devmax = ( devmax < zE2[i] )? zE2[i] : devmax;
    devmin = (zE2[i] < devmin) ? zE2[i] : devmin;
  }

  //std::cout << "devmin : " << devmin << std::endl;
  //std::cout << "devmax : " << devmax << std::endl;
  bool is_planar = (std::max(fabs(devmin), fabs(devmax)) < EPSILON);
  //std::cout << "is planar ? " << is_planar << std::endl;
  
  // if the contour is planar no need for a fine patch
  E_Float dx = 0.2 *K_FUNC::E_min((maxB[0] - minB[0]), (maxB[1] - minB[1]));
  //std::cout << "refine : " << refine << std::endl;
  //std::cout << "bump_factor : " << bump_factor << std::endl;
  //std::cout << "is_planar : " << is_planar << std::endl;
  if (refine || (bump_factor != 0.) || !is_planar)
  {
    dx = __computeCharacteristicLength(pos2D, connectE2);
  }

  /*std::cout << "dx : " << dx << std::endl;
  std::cout << "dx0 : " << maxB[0] - minB[0] << std::endl;
  std::cout << "dy0 : " << maxB[1] - minB[1] << std::endl;*/
  
  // Compute ni and nj
  maxB[0] += 2. * dx; // Enlarge the plaster to ensure to have to rank of nodes outside the domain.
  maxB[1] += 2. * dx;
  minB[0] -= 2. * dx;
  minB[1] -= 2. * dx;  
  
  E_Float nif = 1. + (maxB[0] - minB[0]) / dx;
  E_Float njf = 1. + (maxB[1] - minB[1]) / dx;
  nif *= fabs(bump_factor) + 1.; // 2 times more if factor is 1 or -1.
  njf *= fabs(bump_factor) + 1.;
#ifdef E_ADOLC
  ni = E_Int(nif.value());
  nj = E_Int(njf.value());
#else
  ni = E_Int(nif);
  nj = E_Int(njf);
#endif
  ni = std::min(ni, NIJMAX);
  nj = std::min(nj, NIJMAX);

  // Generate the plaster (a cartesian mesh) on the top side 
  minB[2] = z0;
  __cartesian(minB, maxB, ni, nj, plaster2D); // Create the plaster.
  plaster2D.resize(2, plaster2D.cols());

  //std::cout << "plaster 9" << std::endl;

  // Triangulate the projected contour.
  DELAUNAY::MesherMode mode;
  mode.mesh_mode = mode.TRIANGULATION_MODE;
  DELAUNAY::T3Mesher<E_Float> mesher(mode);
  //K_FLD::IntArray connectT3;
  DELAUNAY::MeshData dataTmp(pos2D, connectE2);
  //data.hardNodes = hN;
  mesher.seed_random(1);
  err = mesher.run(dataTmp);
  if (err) return err;

  // Initialize the plaster field.
  std::vector<E_Float> z(ni*nj, -NUGA::FLOAT_MAX);
  __initializePlaster(plaster2D, ni, pos2D, connectE2, /*hN,*/ zE2, dataTmp.connectM, z, bump_factor);

  //std::cout << "plaster 10 : z field " << z.size() << std::endl;

  // Smooth the z-field.
  __smooth(z, ni, bump_factor, 1.e-6);

  //std::cout << "plaster 11" << std::endl;

  // Set the correct z to each plaster node.
  plaster = plaster2D;
  plaster.resize(3, plaster.cols(), &z0);

#ifdef WIN32
#ifdef E_DEBUG2
  K_FLD::FloatArray pp = posE2;
  pp.pushBack(plaster);
  meshIO::write("init_plaster.mesh", pp, connectE2);
#endif
#endif

  for (E_Int i = 0; i < plaster.cols(); ++i)
    plaster(2, i) -= z[i];

  // Transform back to the real space.
  NUGA::transform(plaster, P);

  //std::cout << "plaster 12" << std::endl;

#ifdef WIN32
#ifdef E_DEBUG
  {
    NUGA::transform(posE2, P);
    K_FLD::FloatArray pp = posE2;
    NUGA::int_vector_type new_IDs;
    NUGA::MeshTool::compact_to_mesh(pp, connectE2, new_IDs);
    pp.pushBack(plaster);
    meshIO::write("plaster.mesh", pp, connectE2);
  }
#endif
#endif

  //std::cout << "plaster 13" << std::endl;

  return 0;
}

///
void
Plaster::__cartesian
(const E_Float* minB, const E_Float* maxB, E_Int ni, E_Int nj, K_FLD::FloatArray& cart)
{
  E_Float I[] = {1.,0.,0.} , J[] = {0.,1.,0.}, K[] = {0., 0., 1.};
  E_Float* Xi, *Xj, dx, dy, Pi[3];

  cart.clear();

  if (minB[0] == maxB[0])
  {
    //std::cout << "I" << std::endl;
    Xi = J;
    dx = maxB[1] - minB[1];
    Xj = K;
    dy = maxB[2] - minB[2];
  }
  else if (minB[1] == maxB[1])
  {
    //std::cout << "J" << std::endl;
    Xi = I;
    dx = maxB[0] - minB[0];
    Xj = K;
    dy = maxB[2] - minB[2];
  }
  else if (minB[2] == maxB[2])
  {
    //std::cout << "K" << std::endl;
    Xi = I;
    dx = maxB[0] - minB[0];
    Xj = J;
    dy = maxB[1] - minB[1];
  }
  else
    return; // error because not handled yet...

  dx /= ni;
  dy /= nj;

  for (E_Int k = 0; k < 3; ++k)
  {
    Xi[k] *= dx;
    Xj[k] *= dy;
  }

  cart.reserve(3, ni*nj);

  for (E_Int j = 0; j < nj; ++j)
  {
    for (E_Int i = 0; i < ni; ++i)
    {
      for (E_Int k = 0; k < 3; ++k)
        Pi[k] = minB[k] + i * Xi[k] + j * Xj[k];

      cart.pushBack(Pi, Pi+3);
    }
  }
}

///
void
Plaster::__smooth
(std::vector<E_Float>& z, E_Int ni, E_Float bump_factor, E_Float tol)
{
  if (K_FUNC::fEqualZero(bump_factor))
    __smooth_1(z, ni, tol);
  else
    __smooth_2(z, ni, tol);
}

///
void
Plaster::__smooth_1
(std::vector<E_Float>& z, E_Int ni, E_Float tol)
{
  E_Int iter(0), iterMax(5000);
  E_Int indH, indB, indD, indG;
  E_Int NBPOINTS(E_Int(z.size())), ind, J;
  E_Float threshold(tol), dMax, d, q;
  bool carry_on(true);
  bool_vector_type processed(NBPOINTS, false);
  
  // Reset nodes to be computed to 0.
  for (E_Int i = 0; i < NBPOINTS; ++i)
  {
    processed[i] = (z[i] != -NUGA::FLOAT_MAX);
    if (!processed[i]) z[i] = 0.;
  }

  while (carry_on)
  {
    dMax = -NUGA::FLOAT_MAX;

    for (ind = 0; ind < NBPOINTS; ++ind)
    {
      if (processed[ind]) continue;

      J = ind % (ni);

      indH = (ind + ni) < NBPOINTS ? ind + ni : ind;
      indB = (ind - ni) >= 0 ? ind - ni : ind;
      indG = (J > 0) ? ind - 1 : ind ;
      indD = (J < ni-1) ? ind + 1 : ind;

      q = 0.25 * ( z[indH] + z[indB] + z[indG] + z[indD] );
      
      d = fabs(z[ind] - q);
      dMax = (dMax < d) ? d : dMax;
      z[ind] = q;
    }

    carry_on = (++iter < iterMax) && (dMax > threshold);
  }
}


///
void
Plaster::__smooth_2
(std::vector<E_Float>& z, E_Int ni, E_Float tol)
{
  E_Int iter(0), iterMax(5000);
  E_Int indH, indB, indD, indG, indH2, indB2, indD2, indG2;
  E_Int NBPOINTS(E_Int(z.size())), ind, J;
  E_Float threshold(tol), dMax, d, q;
  bool carry_on(true);
  bool_vector_type processed(NBPOINTS, false);

   E_Float zh, zb, zg, zd, k(0.333333), k1(1.+k);
  
  // Reset nodes to be computed to 0.
  for (E_Int i = 0; i < NBPOINTS; ++i)
  {
    processed[i] = (z[i] != -NUGA::FLOAT_MAX);
    if (!processed[i]) z[i] = 0.;
  }

  while (carry_on)
  {
    dMax = -NUGA::FLOAT_MAX;

    for (ind = 0; ind < NBPOINTS; ++ind)
    {
      if (processed[ind]) continue;

      J = ind % (ni);

      indH = (ind + ni) < NBPOINTS ? ind + ni : ind;
      indB = (ind - ni) >= 0 ? ind - ni : ind;
      indG = (J > 0) ? ind - 1 : ind ;
      indD = (J < ni-1) ? ind + 1 : ind;

      indH2 = (indH + ni) < NBPOINTS ? indH + ni : indH;
      indB2 = (indB - ni) >= 0 ? indB - ni : indB;
      indG2 = ((indG % (ni)) > 0) ? indG - 1 : indG ;
      indD2 = ((indD % (ni)) < ni-1) ? indD + 1 : indD;

      zh = (k1 * z[indH]) - (k * z[indH2]);  
      zb = (k1 * z[indB]) - (k * z[indB2]);
      zg = (k1 * z[indG]) - (k * z[indG2]);
      zd = (k1 * z[indD]) - (k * z[indD2]);

      q = 0.25 * (zh +zb +zg +zd);

      d = fabs(z[ind] - q);
      dMax = (dMax < d) ? d : dMax;
      z[ind] = q;
    }

    carry_on = (++iter < iterMax) && (dMax > threshold);
  }

  //std::cout << "iter : " << iter << std::endl;
}

///
void
Plaster::__initializePlaster
(const K_FLD::FloatArray& plaster2D, E_Int ni, const K_FLD::FloatArray& pos2D,
 const K_FLD::IntArray& connectE2/*, const std::vector<E_Int>& hard_nodes*/, const std::vector<E_Float>& zE2,
 const K_FLD::IntArray& connectT3, std::vector<E_Float>& z, E_Float bump_factor)
{
  //Fast return
  if (plaster2D.cols() == 0) return;
  if (ni == 0)               return;
  if (pos2D.cols() == 0)     return;
  if (connectE2.cols() == 0) return;
  if (connectT3.cols() == 0) return;
  if (zE2.empty())           return;

  // Mask the plaster nodes that are inside (strictly) the domain.
  // O = inside / 1 = outside.
  bool_vector_type mask;
  __mask(pos2D, plaster2D, connectT3, mask);
  
  // Block the plaster nodes surrounding the contour.
  NUGA::int_set_type onodes;
  __blockSurroundingNodes(plaster2D, ni, mask, pos2D, connectE2, zE2, z, onodes);

  // Block the inner nodes close to har nodes
  //__blockInsideNodes(plaster2D, ni, pos2D, hard_nodes, zE2, z);
  
  // Bump the plaster.
  __bumpPlaster(plaster2D, ni, mask, bump_factor, onodes, z);
}

///
void
Plaster::__blockSurroundingNodes
(const K_FLD::FloatArray& plaster2D, E_Int ni, 
 const NUGA::bool_vector_type& mask,
 const K_FLD::FloatArray& pos2D,
 const K_FLD::IntArray& connectE2, const std::vector<E_Float>& zE2,
 std::vector<E_Float>& z, NUGA::int_set_type& onodes)
{
  // Get the plaster nodes that are immediately surrounding the domain.
  __getPlasterBoundary(mask, ni, /*outside*/true, onodes);

#ifdef WIN32
#ifdef E_DEBUG
  K_FLD::FloatArray blocked = pos2D;
  for (std::set<E_Int>::const_iterator i = onodes.begin(); i != onodes.end(); ++i)
  {blocked.pushBack(plaster2D.col(*i), plaster2D.col(*i)+2);}
  meshIO::write("blocked.mesh", blocked, connectE2);
#endif
#endif

  // Set the z value for those nodes based on an interpolation of the closest contour's edge.
  __blockNodes(pos2D, plaster2D, connectE2, zE2, onodes, z);
}
/*
///
void
Plaster::__blockInsideNodes
(const K_FLD::FloatArray& plaster2D, E_Int ni, 
 const K_FLD::FloatArray& pos2D, const std::vector<E_Int>& hard_nodes,
 const std::vector<E_Float>& zE2, std::vector<E_Float>& z)
{
  if (hard_nodes.empty())
    return;

  K_FLD::ArrayAccessor<K_FLD::FloatArray> pA(plaster2D);
  K_SEARCH::KdTree<> ptree(pA);

  for (size_t i = 0; i < hard_nodes.size(); ++i)
  {
    E_Int N = ptree.getClosest(pos2D.col(hard_nodes[i]));
    z[N] = zE2[hard_nodes[i]];
  }
}
*/
///
void
Plaster::__bumpPlaster
(const K_FLD::FloatArray& plaster2D, E_Int ni,
 const NUGA::bool_vector_type& mask,
 E_Float bump_factor, const NUGA::int_set_type& onodes,
 std::vector<E_Float>& z)
{
  bump_factor = K_FUNC::E_max(bump_factor, -1.); // factor must be in [-1., 1.]
  bump_factor = K_FUNC::E_min(bump_factor, 1.);

  if (K_FUNC::fEqualZero(bump_factor)) return;

  const E_Float BUMP_ANGLE_MAX = 1.5 * NUGA::PI_4; //3PI/8
  E_Float ta = tan(bump_factor * BUMP_ANGLE_MAX);

  std::vector<E_Int> oonodes;
  oonodes.insert(oonodes.end(), onodes.begin(), onodes.end());

  K_FLD::ArrayAccessor<K_FLD::FloatArray> pA(plaster2D);
  K_SEARCH::KdTree<> tree(pA, oonodes);

  int_set_type inodes;
  __getPlasterBoundary(mask, ni, false/*inside*/, inodes);

  for (int_set_type::iterator it = inodes.begin(); it != inodes.end(); ++it)
  {
    E_Int N = tree.getClosest(plaster2D.col(*it));
    E_Float dx = sqrt(NUGA::sqrDistance(plaster2D.col(N), plaster2D.col(*it), 2));
    z[*it] = z[N] + (ta * dx);
  }
}


///
bool Plaster::__IsStrictlyInT3
(const E_Float* P, const E_Float* P0, const E_Float* P1, const E_Float* P2)
{
  E_Float s;
  
  s = K_MESH::Triangle::surface(P0, P1, P, 2);
  if (s < -EPSILON) return false;

  s = K_MESH::Triangle::surface(P1, P2, P, 2);
  if (s < -EPSILON) return false;
  
  s = K_MESH::Triangle::surface(P2, P0, P, 2);
  if (s < -EPSILON) return false;

  return true;
}

///
void
Plaster::__mask
(const K_FLD::FloatArray& pos2D, const K_FLD::FloatArray& plaster2D,
 const K_FLD::IntArray& connectT3, NUGA::bool_vector_type& mask)
{
  const E_Int NB_POINTS(plaster2D.cols());
  E_Int S;
  std::vector<E_Int> T3s;
  const E_Float *P, *P0, *P1, *P2;
  K_FLD::IntArray::const_iterator pS;

  typedef K_SEARCH::BoundingBox<2> BBox2DType;

  mask.clear();
  mask.resize(NB_POINTS, true);

  if (pos2D.cols() == 0) return;
  if (connectT3.cols() == 0) return;
  if (plaster2D.cols() == 0) return;

  // Build the boxes;
  std::vector<BBox2DType*> boxes(connectT3.cols());
  for (size_t i = 0; i < boxes.size(); ++i)
    boxes[i] = new BBox2DType(pos2D, connectT3.col(i), 3);
  // Build the box tree.
  K_SEARCH::BbTree2D tree(boxes);

  for (E_Int i = 0; i < NB_POINTS; ++i)
  {
    T3s.clear();
    P = plaster2D.col(i);
    tree.getOverlappingBoxes(P, P, T3s);

    for (size_t k = 0; (k < T3s.size()) && mask[i]; ++k)
    {
      S = T3s[k];
      pS = connectT3.col(S);
      P0 = pos2D.col(*pS);
      P1 = pos2D.col(*(pS+1));
      P2 = pos2D.col(*(pS+2));

      mask[i] = !__IsStrictlyInT3(P, P0, P1, P2);
    }
  }

  // Final cleaning.
  for (size_t i = 0; i < boxes.size(); ++i) delete boxes[i];

#ifdef WIN32
#ifdef E_DEBUG
  K_FLD::FloatArray pp;
  for (size_t i = 0; i < mask.size(); ++i)
    if (mask[i])
      pp.pushBack(plaster2D.col(i), plaster2D.col(i)+2);

  meshIO::write("mask.mesh", pp);
  meshIO::write("tria.mesh", pos2D, connectT3);
#endif
#endif
}

///
void
Plaster::__getPlasterBoundary
(const NUGA::bool_vector_type& mask, E_Int ni, bool outside,
 NUGA::int_set_type& nodes)
{
  E_Int ind, indH, indB, indG, indD, J, NB_POINTS(mask.size());
  for (ind = 0; ind < NB_POINTS; ++ind)
  {
    if (mask[ind] == outside) continue;

    J = ind % (ni);

    indH = (ind + ni) < NB_POINTS ? ind + ni : ind;
    indB = (ind - ni) >= 0 ? ind - ni : ind;
    indG = (J > 0) ? ind - 1 : ind ;
    indD = (J < ni-1) ? ind + 1 : ind;

    if (mask[indH] == outside)
      nodes.insert(indH);
    if (mask[indB] == outside)
      nodes.insert(indB);
    if (mask[indG] == outside)
      nodes.insert(indG);
    if (mask[indD] == outside)
      nodes.insert(indD);
  }
}

void
Plaster::__blockNodes
(const K_FLD::FloatArray& pos2D, const K_FLD::FloatArray& plaster2D,
 const K_FLD::IntArray& connectE2, const std::vector<E_Float>& zE2,
 const NUGA::int_set_type& nodes, std::vector<E_Float>& z)
{
  typedef K_SEARCH::BoundingBox<2>  BBox2DType;
  
  // Build the boxes;
  std::vector<BBox2DType*> E2boxes(connectE2.cols());
  for (size_t i = 0; i < E2boxes.size(); ++i)
    E2boxes[i] = new BBox2DType(pos2D, connectE2.col(i), 2);
  // Build the box tree.
  K_SEARCH::BbTree2D E2tree(E2boxes);

  E_Float dx = sqrt(2. * NUGA::sqrDistance(plaster2D.col(0), plaster2D.col(1), 2));
  std::vector<E_Int> E2s;
  E_Float mB[2], MB[2], min_d, de, lambda;
  E_Int N;
  const E_Float* P;
  for (std::set<E_Int>::const_iterator i = nodes.begin(); i != nodes.end(); ++i)
  {
    E2s.clear();
    N = *i;
    P = plaster2D.col(N);
    mB[0] = plaster2D(0, N) - dx;
    mB[1] = plaster2D(1, N) - dx;
    MB[0] = plaster2D(0, N) + dx;
    MB[1] = plaster2D(1, N) + dx;

    E2tree.getOverlappingBoxes(mB, MB, E2s);

    assert(!E2s.empty());

    min_d = NUGA::FLOAT_MAX;
    for (size_t e = 0; e < E2s.size(); ++e)
    {
      E_Int Ni = connectE2(0, E2s[e]);
      E_Int Nj = connectE2(1, E2s[e]);
      de = K_MESH::Edge::edgePointMinDistance<2>(pos2D.col(Ni), pos2D.col(Nj), P, lambda);
      if (de < min_d)
      {
        min_d = de;
        z[N] = (1. - lambda) * zE2[Ni] + lambda * zE2[Nj];
      }
    }
  }

  // Final cleaning.
  for (size_t i = 0; i < E2boxes.size(); ++i)
    delete E2boxes[i];
}

E_Float
Plaster::__computeCharacteristicLength
(const K_FLD::FloatArray& pos2D, const K_FLD::IntArray& connectE2)
{
  E_Float min_d = NUGA::FLOAT_MAX, max_d = -NUGA::FLOAT_MAX, L,perimeter(0.);
  K_FLD::FloatArray Lengths;
  NUGA::MeshTool::computeEdgesSqrLengths<2>(pos2D, connectE2, Lengths);
  for (E_Int l = 0; l < Lengths.cols(); ++l)
  {
    L = sqrt(Lengths(0,l));
    min_d = K_FUNC::E_min(min_d, L);
    max_d = K_FUNC::E_max(max_d, L);
    perimeter += L;
  }
  return perimeter/connectE2.cols()/*0.5*(min_d+max_d)*/;
}

#ifdef WIN32
#ifdef E_DEBUG

void Plaster::make_box(const E_Float* minB, const E_Float* maxB, K_FLD::FloatArray& boxPs, K_FLD::IntArray& boxC)
{
  boxPs.resize(3, 8);
  boxC.resize(2, 12);

  boxPs(0,0) = minB[0];boxPs(1,0) = minB[1];boxPs(2,0) = minB[2];
  boxPs(0,1) = maxB[0];boxPs(1,1) = minB[1];boxPs(2,1) = minB[2];
  boxPs(0,2) = maxB[0];boxPs(1,2) = maxB[1];boxPs(2,2) = minB[2];
  boxPs(0,3) = minB[0];boxPs(1,3) = maxB[1];boxPs(2,3) = minB[2];
  
  boxPs(0,4) = minB[0];boxPs(1,4) = minB[1];boxPs(2,4) = maxB[2];
  boxPs(0,5) = maxB[0];boxPs(1,5) = minB[1];boxPs(2,5) = maxB[2];
  boxPs(0,6) = maxB[0];boxPs(1,6) = maxB[1];boxPs(2,6) = maxB[2];
  boxPs(0,7) = minB[0];boxPs(1,7) = maxB[1];boxPs(2,7) = maxB[2];

  boxC(0,0) = 0;boxC(1,0) = 1;
  boxC(0,1) = 1;boxC(1,1) = 2;
  boxC(0,2) = 2;boxC(1,2) = 3;
  boxC(0,3) = 3;boxC(1,3) = 0;
  
  boxC(0,4) = 4;boxC(1,4) = 5;
  boxC(0,5) = 5;boxC(1,5) = 6;
  boxC(0,6) = 6;boxC(1,6) = 7;
  boxC(0,7) = 7;boxC(1,7) = 4;

  boxC(0,8) = 0;boxC(1,8) = 4;
  boxC(0,9) = 1;boxC(1,9) = 5;
  boxC(0,10) = 2;boxC(1,10) = 6;
  boxC(0,11) = 3;boxC(1,11) = 7;
}

void Plaster::drawBox(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, const E_Float* mB, const E_Float* MB)
{
  K_FLD::FloatArray pp1(pos), pp;
  K_FLD::IntArray cc1(connect), cc;
  make_box(mB, MB, pp, cc);
  cc.shift(pos.cols());
  cc1.pushBack(cc);
  pp1.pushBack(pp);

  std::ostringstream o;
  o << "box_" << _count++ << ".mesh";

  meshIO::write(o.str().c_str(), pp1, cc1);
}
#endif
#endif

