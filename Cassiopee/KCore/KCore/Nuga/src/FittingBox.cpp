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

#include "Nuga/include/FittingBox.h"
#include "Nuga/include/defs.h"
#include "Nuga/include/MeshTool.h"
#include "Nuga/include/T3Mesher.h"
#ifdef DEBUG_FITTINGBOX
#include "IO/DynArrayIO.h"
#endif

#define ROUND(x) (((x<EPSILON) && (x>-EPSILON)) ? 0.: x)

using namespace NUGA;

///
void
FittingBox::computeNormalToContour
(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, E_Float* W)
{
  NUGA::int_vector_type nodes;
  E_Int nb_nodes, nb_elts(connect.cols()), Ni, Nj;

  K_FLD::FloatArray wPos(pos);

  connect.uniqueVals(nodes);
  nb_nodes = nodes.size();

  // Compute the Barycenter
  E_Float G[] = {0., 0., 0.};
  for (E_Int i = 0; i < nb_nodes; ++i)
    NUGA::sum<3>(pos.col(nodes[i]), G, G);

  for (E_Int i = 0; i < 3; ++i)
    G[i] /= nb_nodes;

  // Compute an approximate normal W to the contour's surface (oriented toward outside).
  E_Float V1[3], V2[3], w[3];
  for (E_Int i = 0; i < 3; ++i) W[i] = 0.;

  for (E_Int c = 0; c < nb_elts; ++c)
  {
    Ni = connect(0, c);
    Nj = connect(1, c);

    NUGA::diff<3>(pos.col(Ni), G, V1);
    NUGA::diff<3>(pos.col(Nj), pos.col(Ni), V2);
    
    // prevent numerical error when computing cross product. fixme : should be done evrywhere a cross product or determinant is done ?
    for (size_t i=0; i < 3; ++i)
    {
      V1[i]=ROUND(V1[i]);
      V2[i]=ROUND(V2[i]);
    }

    NUGA::crossProduct<3>(V1, V2, w);
    NUGA::sum<3>(w, W, W);
  }
  NUGA::normalize<3>(W);
}

///
void
FittingBox::computeFittingFrame
(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect,
 const E_Float* W, K_FLD::FloatArray& P)
{
  P.clear();

  // Compute a coord system with W as Z-axis.
  computeAFrame(W, P);

  // Transform the working coordinates.
  K_FLD::FloatArray iP(P);
  K_FLD::FloatArray::inverse3(iP);

  K_FLD::FloatArray wPos(pos);
  transform (wPos, iP);
  //E_Float minB[3], maxB[3];
  //NUGA::MeshTool::boundingBox(wPos, connect, minB, maxB);
  //drawBox(wPos, connect, minB, maxB);

  // Now rotate the contour around W to reduce the box.
  __fitByRotating(wPos, connect, P);
}

///
E_Int
FittingBox::computeOptimalViewFrame
(const K_FLD::FloatArray& posE2, const K_FLD::IntArray& connectE2, K_FLD::FloatArray& P)
{
  K_FLD::FloatArray posFrame, iP(3,3), R(3, 3);
  E_Float W[3];
  E_Int err;
  E_Float alphamax = NUGA::PI_2;
  E_Float alpha(2.* alphamax/E_Float(_maxStep));

  FittingBox::computeNormalToContour(posE2, connectE2, W);

  // Rotation around x-axis
  R(0,0) = R(1,0) = R(2,1) = R(2,2) = 0.;
  R(2,0) = 1.;
  R(0, 1) = R(1, 2) = cos(alpha);
  R(1, 1) = sin(alpha);
  R(0, 2) = - R(1, 1);
  err = __computeOptimalViewFrame(posE2, connectE2, W, R, P); // Try rotating aound X-axis
  if (err)
  {
    // Rotation around y-axis
    R(0,1) = R(1,1) = R(2,0) = R(2,2) = 0.;
    R(2,1) = 1.;
    R(0, 0) = R(1, 2) = cos(alpha);
    R(1, 0) = sin(alpha);
    R(0, 2) = - R(1, 0);
    err = __computeOptimalViewFrame(posE2, connectE2, W, R, P); // Try rotating around y-axis
  }
  if (err) return err;

  for (E_Int k = 0; k < 3; ++k) // Best direction found.
    W[k] = P(k,2);

  FittingBox::computeFittingFrame(posE2, connectE2, W, P); // Reduce the bounding box.

  return 0;
}

///
E_Int
FittingBox::__computeOptimalViewFrame
(const K_FLD::FloatArray& posE2, const K_FLD::IntArray& connectE2, const E_Float* W0, const K_FLD::FloatArray& R, K_FLD::FloatArray& P)
{
  E_Float W[3];
  bool carry_on = true;
  E_Int error(0), err, iter(0);
  K_FLD::FloatArray posFrame, iP(3,3);

  DELAUNAY::MesherMode mode;
  mode.mesh_mode = mode.TRIANGULATION_MODE;
  DELAUNAY::T3Mesher<E_Float> mesher(mode);

  for (E_Int k = 0; k < 3; ++k) W[k] = W0[k];

  while (carry_on)
  {
    error = 0;
    posFrame = posE2;
    NUGA::computeAFrame(W, P);
    iP = P;
    K_FLD::FloatArray::inverse3(iP);
    NUGA::transform(posFrame, iP); // Transform to computed frame.

#ifdef DEBUG_FITTINGBOX
    K_CONVERTER::DynArrayIO::write("transformed.mesh", posFrame, connectE2, "BAR");
#endif

    posFrame.resize(2, posFrame.cols()); // project.

#ifdef DEBUG_FITTINGBOX
    K_CONVERTER::DynArrayIO::write("Proj_Contour.mesh", posFrame, connectE2, "BAR");
#endif

    // Triangulate the projection.
    DELAUNAY::MeshData dataTmp(posFrame, connectE2);
    //DELAUNAY::iodata::write("out.dat", dataTmp);
    //data.hardNodes = hN;
    mesher.seed_random(1);
    err = mesher.run(dataTmp);

    if (err || (dataTmp.connectM.cols() == 0))
    {
      error = 1;
      // Rotate
      NUGA::transform(P, R);
      W[0] = P(0, 2);
      W[1] = P(1, 2);
      W[2] = P(2, 2);
    }

    carry_on = (error && iter++ < _maxStep);
  }

  return error;
}


///
void
FittingBox::__fitByRotating
(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, K_FLD::FloatArray& P)
{
  E_Float             minB[3], maxB[3], s0,s1, min_s;
  E_Float             alpha(NUGA::PI_4/_maxStep);
  E_Float             G[] = {0., 0., 0.};
  K_FLD::FloatArray   wpos0(pos), wpos;
  E_Float             sign = +1., a(0);
  K_FLD::FloatArray   R(3,3, 0.), best_R(3,3,0.);

  /* Init */
  // Initial bounding box size.
  NUGA::MeshTool::boundingBox(wpos0, connect, minB, maxB);
  min_s = s0 = (maxB[1] - minB[1])*(maxB[0] - minB[0]);
  // Compute the Barycenter.
  for (E_Int i = 0; i < wpos0.cols(); ++i)
    NUGA::sum<3>(wpos0.col(i), G, G);
  // Translate.
  for (E_Int i = 0; i < wpos0.cols(); ++i)
    for (E_Int k = 0; k < 3; ++k)
      wpos0(k, i) -= G[k];

  best_R(0,0) = best_R(1,1) = best_R(2,2) = 1.0; // Identity matrix.
  R(2,2) = 1.;                                   // alpha-rotation matrix.
  R(0,0) = R(1,1) = cos(alpha);
  R(1,0) = sin(alpha);
  R(0,1) = - R(1,0);

  /* Find out the reducing direction. */
  wpos = wpos0;
  transform(wpos, R);
  NUGA::MeshTool::boundingBox(wpos, connect, minB, maxB);
  s1 = (maxB[1] - minB[1])*(maxB[0] - minB[0]);
  if (s1 > s0)
    sign = -1;

  /* get the rotation which minimize the bounding box. */
  for (E_Int i = 1; i < _maxStep; ++i)
  {
    wpos = wpos0;
    a = sign * alpha * i;

    // Rotation matrix.
    R(0,0) = R(1,1) = cos(a);
    R(1,0) = sin(a);
    R(0,1) = - R(1,0);

    // Rotate
    transform(wpos, R);

    // New bounding box.
    NUGA::MeshTool::boundingBox(wpos, connect, minB, maxB);
    s1 = (maxB[1] - minB[1])*(maxB[0] - minB[0]);

    if (s1 > min_s)
      break;

    min_s = s1;
    best_R = R;

    //NUGA::MeshTool::boundingBox(wpos, connect, minB, maxB);
    //drawBox(wpos, connect, minB, maxB);
  }

  best_R.transpose();
  P = P*best_R;
}

#ifdef E_DEBUG

void FittingBox::drawBox(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, const E_Float* mB, const E_Float* MB)
{
  K_FLD::FloatArray pp1(pos), pp;
  K_FLD::IntArray cc1(connect), cc;
  make_box(mB, MB, pp, cc);
  cc.shift(pos.cols());
  cc1.pushBack(cc);
  pp1.pushBack(pp);

  meshIO::write("box.mesh", pp1, cc1);
}

void FittingBox::make_box(const E_Float* minB, const E_Float* maxB, K_FLD::FloatArray& boxPs, K_FLD::IntArray& boxC)
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

#endif

