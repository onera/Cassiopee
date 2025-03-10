/*    
    Copyright 2013-2025 Onera.

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

#ifndef __MERGE_H__
#define __MERGE_H__

#include "Nuga/include/defs.h"
#include "Nuga/include/ArrayAccessor.h"
#include "Nuga/include/KdTree.h"
#include "Nuga/include/Triangle.h"
#include "Nuga/include/openMP.h"

//namespace K_CONNECT
//{

// remove duplicates without changing order, give the set of unique values upon exit
template <typename T>
inline void stable_unique_values
(const std::vector<T>& source, std::vector<T>& result)
{
  std::set<T> unics;
  result.clear();
  result.reserve(source.size());
  size_t sz = source.size();
  for (size_t i = 0; i < sz; ++i)
  {
    if (unics.insert(source[i]).second)
      result.push_back(source[i]);
  }
}


struct triplet_t
{
  triplet_t()=default;
  triplet_t(double id, E_Int iNi, E_Int iNj):d(id), Ni(iNi), Nj(iNj){}

  double d; E_Int Ni, Nj;
};

inline bool lower_than(const triplet_t& p1, const triplet_t& p2)
{return p1.d < p2.d;}

//==============================================================================
template <typename ArrayType>
void 
merge_no_order_omp
(const K_FLD::ArrayAccessor<ArrayType>& coordAcc, E_Float tol,
 NUGA::int_vector_type& new_IDs)
{
  K_SEARCH::KdTree<ArrayType> tree(coordAcc, EPSILON, true/*do_omp*/);
  E_Int npts = coordAcc.size();
  new_IDs.resize(npts);
  
#pragma omp parallel for default (shared)
  for (E_Int i = 0; i < npts; i++)
  {
    new_IDs[i] = i;
  }

  E_Float tol2 = tol*tol;

  // Premiere version - elle ne thread pas
  /*

//#pragma omp parallel default(shared)
    {
      std::vector<E_Int> onodes;
      E_Int m;
      E_Float Xn[3]; E_Float Xm[3];
      size_t size; E_Float d2;
      E_Int mID, iID;

//#pragma omp for schedule(dynamic)
      for (E_Int i = 0; i < npts; i++)
      {
        coordAcc.getEntry(i, Xn);
        //printf("pt %f %f %f\n", Xn[0], Xn[1], Xn[2]);
        onodes.clear();
        tree.getInSphere(Xn, tol, onodes);
        size = onodes.size();
        //printf("cloud size %d\n", size);
        if (size > 0)
        {
          //m = *std::min_element(onodes.begin(), onodes.end());
          //std::sort(onodes.begin(), onodes.end());
          iID = new_IDs[i];
         
          for (E_Int j = 0; j < size; j++)
          {
            m = onodes[j];
            mID = new_IDs[m]; // firewall needed
            if (i > m)
            {
              if (mID == m) { new_IDs[i] = m; break; }
              else if (mID != i)
              {
                coordAcc.getEntry(mID, Xm);
                d2 = (Xm[0]-Xn[0])*(Xm[0]-Xn[0])+
                  (Xm[1]-Xn[1])*(Xm[1]-Xn[1])+
                  (Xm[2]-Xn[2])*(Xm[2]-Xn[2]);
                if (d2 < tol2) { new_IDs[i] = mID; break; }
              }
            }
          }
          
        }
       

      }
    }
    printf("ok\n");
  */
    
    // nouvelle version --
    std::vector< std::vector<E_Int> > allnodes(npts);
#pragma omp parallel default(shared)
    {
      E_Float Xn[3]; 
     
#pragma omp for schedule(dynamic)
      for (E_Int i = 0; i < npts; i++)
      {
        coordAcc.getEntry(i, Xn);
        //printf("pt %f %f %f\n", Xn[0], Xn[1], Xn[2]);
        tree.getInSphere(Xn, tol, allnodes[i]);
      }
    }

    E_Int m; E_Float d2;
    E_Float Xn[3]; E_Float Xm[3];
    E_Int mID;  size_t size; 
    
    for (E_Int i = 0; i < npts; i++)
    {
      coordAcc.getEntry(i, Xn);
      std::vector<E_Int>& onodes = allnodes[i]; 
      size = onodes.size();
      //printf("cloud size %d\n", size);
      if (size > 0)
      {
        for (size_t j = 0; j < size; j++)
        {
          m = onodes[j];
          mID = new_IDs[m]; // firewall needed
          if (i > m)
          {
            if (mID == m) { new_IDs[i] = m; break; }
            else if (mID != i)
            {
              coordAcc.getEntry(mID, Xm);
              d2 = (Xm[0]-Xn[0])*(Xm[0]-Xn[0])+
                (Xm[1]-Xn[1])*(Xm[1]-Xn[1])+
                (Xm[2]-Xn[2])*(Xm[2]-Xn[2]);
              if (d2 < tol2) { new_IDs[i] = mID; break; }
            }
          }
        }
        
      }
    }
    
    // Fin deuxieme version
}

//==============================================================================
template <typename ArrayType>
E_Int
__merge
(const K_FLD::ArrayAccessor<ArrayType>& coordAcc, E_Float tol,
 NUGA::int_vector_type& umoving,
 NUGA::int_vector_type& utarget,
 NUGA::int_vector_type& new_IDs)
{
  size_t ssz = umoving.size();
  size_t tsz = utarget.size();
  size_t nsz = coordAcc.size();
  
  new_IDs.resize(nsz);
  for (size_t i = 0; i < nsz; ++i) new_IDs[i]=i;

  // Fast return
  if (ssz*tsz*nsz == 0) return 0;

  // build the KdTree on moving nodes.
  K_SEARCH::KdTree<ArrayType> moving_tree(coordAcc, umoving);
  
  //E_Float d2, tol2(tol*tol);
  typedef std::vector<triplet_t > palmares_t;
  palmares_t palma;

  std::vector<E_Int> onodes;
  std::vector<E_Float> dist2;

  // loop on target nodes and get all moving nodes in sphere of radius tol.
  size_t osz;
  for (size_t i = 0; i < tsz; ++i)
  {
    onodes.clear(); dist2.clear();

    E_Int& Fi = utarget[i];
    moving_tree.getInSphere(Fi, tol, onodes, dist2);
    osz = onodes.size();
    for (size_t j = 0; j < osz; ++j)
      palma.push_back(triplet_t(dist2[j], Fi, onodes[j]));
  }

  if (palma.empty()) return 0;

  // sort them by increasing distance (keep intial order in case of equality to priorize it).
  std::stable_sort(palma.begin(), palma.end(), lower_than);

  // move the node
  E_Int Fi, Mi, nb_merges(0);
  size_t psz = palma.size();
  for (size_t i = 0; i < psz; ++i)
  {
    Fi=palma[i].Ni;
    Mi=palma[i].Nj;

    // rule of merging : if target and moving has not been moved already
    if ((new_IDs[Fi] == Fi) && (new_IDs[Mi] == Mi))
    {
      new_IDs[Mi]=Fi;
      ++nb_merges;
    }
  }

  // update the pointers to point to the leaves
  for (size_t i =0; i < nsz; ++i)
  {
    Fi=new_IDs[i];
    while (Fi != new_IDs[Fi]) Fi=new_IDs[Fi];
    new_IDs[i]=Fi;
  }

  return nb_merges;
}

//==============================================================================
template <typename ArrayType>
E_Int
__merge_omp
(const K_FLD::ArrayAccessor<ArrayType>& coordAcc, E_Float tol,
  NUGA::int_vector_type& umoving,
  NUGA::int_vector_type& utarget,
  NUGA::int_vector_type& new_IDs)
{
  size_t ssz = umoving.size();
  size_t tsz = utarget.size();
  size_t nsz = coordAcc.size();

  new_IDs.resize(nsz);
  
  for (size_t i = 0; i < nsz; ++i) new_IDs[i] = i;

  // Fast return
  if (ssz*tsz*nsz == 0) return 0;

  // build the KdTree on moving nodes.
  K_SEARCH::KdTree<ArrayType> moving_tree(coordAcc, umoving, EPSILON, true /*do omp*/);

  using palmares_t = std::vector<triplet_t>;

  std::vector<palmares_t> palma_thrd(__NUMTHREADS__);
  palmares_t palma;

  std::vector<std::vector<E_Int>> onodes_thrd(__NUMTHREADS__);
  std::vector<std::vector<E_Float>> dist2_thrd(__NUMTHREADS__);

  E_Int id, Fi, Fj;
  size_t i, j, osz;
  double d2;

  // loop on target nodes and get all moving nodes in sphere of radius tol.
#pragma omp parallel shared(tsz, moving_tree, palma, utarget, tol, palma_thrd, dist2_thrd, onodes_thrd) private (i, j, id, osz, Fi, Fj, d2) default(none)
  {
    id = __CURRENT_THREAD__;

#pragma omp for schedule(static)
    for (i = 0; i < tsz; ++i)
    {
      onodes_thrd[id].clear(); dist2_thrd[id].clear();

      Fi = utarget[i];
      moving_tree.getInSphere(Fi, tol, onodes_thrd[id], dist2_thrd[id]);

      osz = onodes_thrd[id].size();

      for (j = 0; j < osz; ++j)
      {
        d2 = dist2_thrd[id][j];
        Fj = onodes_thrd[id][j];

        palma_thrd[id].push_back(triplet_t(d2, Fi, Fj));
      }
    }
  }

  size_t sz = 0;
  for (E_Int t = 0; t < __NUMTHREADS__; ++t) sz += palma_thrd[t].size();

  palma.reserve(sz);

  for (E_Int t = 0; t < __NUMTHREADS__; ++t)
  {
    //std::cout << "palma_thrd[t] sz : " << palma_thrd[t].size() << std::endl;
    if (!palma_thrd[t].empty())
      palma.insert(palma.end(), ALL(palma_thrd[t]));
  }

  if (palma.empty()) return 0;

  // sort them by increasing distance (keep intial order in case of equality to priorize it).
  std::stable_sort(palma.begin(), palma.end(), lower_than);

  // move the node
  E_Int Mi, nb_merges(0);
  size_t psz = palma.size();
  for (size_t i = 0; i < psz; ++i)
  {
    Fi = palma[i].Ni;
    Mi = palma[i].Nj;

    // rule of merging : if target and moving has not been moved already
    if ((new_IDs[Fi] == Fi) && (new_IDs[Mi] == Mi))
    {
      new_IDs[Mi] = Fi;
      ++nb_merges;
    }
  }
  
  // update the pointers to point to the leaves
  for (size_t i = 0; i < nsz; ++i)
  {
    Fi = new_IDs[i];
    while (Fi != new_IDs[Fi]) Fi = new_IDs[Fi];
    new_IDs[i] = Fi;
  }

  return nb_merges;
}

// // chunk version
// template <typename ArrayType>
// E_Int
// __merge_omp
// (const K_FLD::ArrayAccessor<ArrayType>& coordAcc, E_Float tol,
//   NUGA::int_vector_type& umoving,
//   NUGA::int_vector_type& utarget,
//   NUGA::int_vector_type& new_IDs)
// {
//   size_t ssz = umoving.size();
//   size_t tsz = utarget.size();
//   size_t nsz = coordAcc.size();

//   new_IDs.resize(nsz);
  
//   for (size_t i = 0; i < nsz; ++i) new_IDs[i] = i;

//   // Fast return
//   if (ssz*tsz*nsz == 0) return 0;

//   // build the KdTree on moving nodes.
//   K_SEARCH::KdTree<ArrayType> moving_tree(coordAcc, umoving, EPSILON, true /*do omp*/);

//   const int CHUNK_SZ = 50;
//   E_Int NB_CHUNKS = tsz / CHUNK_SZ;

//   using palmares_t = std::vector<triplet_t>;

//   std::vector<palmares_t> palma_thrd(NB_CHUNKS);
//   palmares_t palma;

//   std::vector<std::vector<E_Int>> onodes_thrd(NB_CHUNKS);
//   std::vector<std::vector<E_Float>> dist2_thrd(NB_CHUNKS);

//   std::vector<E_Int> ibeg(NB_CHUNKS), iend(NB_CHUNKS);
//   for (size_t k=0; k<NB_CHUNKS; ++k)
//   {
//     ibeg[k]=k*CHUNK_SZ;
//     iend[k]=std::min((k+1)*CHUNK_SZ, tsz);
//   }

//   E_Int id, Fi, Fj;
//   size_t i, j, osz, c;
//   double d2;

//   // loop on target nodes and get all moving nodes in sphere of radius tol.
// #pragma omp parallel shared(tsz, moving_tree, palma, utarget, tol, palma_thrd, dist2_thrd, onodes_thrd, NB_CHUNKS, ibeg, iend) private (c, i, j, id, osz, Fi, Fj, d2) default(none)
//   {
//     //id = __CURRENT_THREAD__;

// #pragma omp for schedule(dynamic)
//     for (c = 0; c < NB_CHUNKS; ++c)
//     {
//       id = c;

//       onodes_thrd[id].clear(); dist2_thrd[id].clear();

//       for (size_t i = ibeg[c]; i < iend[c]; ++i)
//       {

//         Fi = utarget[i];
//         moving_tree.getInSphere(Fi, tol, onodes_thrd[id], dist2_thrd[id]);

//         osz = onodes_thrd[id].size();

//         for (j = 0; j < osz; ++j)
//         {
//           d2 = dist2_thrd[id][j];
//           Fj = onodes_thrd[id][j];

//           palma_thrd[id].push_back(triplet_t(d2, Fi, Fj));
//         }
//       }
//     }
//   }

//   size_t sz = 0;
//   for (size_t t = 0; t < palma_thrd.size(); ++t) sz += palma_thrd[t].size();

//   palma.reserve(sz);

//   for (size_t t = 0; t < NB_CHUNKS; ++t)
//   {
//     //std::cout << "palma_thrd[t] sz : " << palma_thrd[t].size() << std::endl;
//     if (!palma_thrd[t].empty())
//       palma.insert(palma.end(), ALL(palma_thrd[t]));
//   }

//   if (palma.empty()) return 0;

//   // sort them by increasing distance (keep intial order in case of equality to priorize it).
//   std::stable_sort(palma.begin(), palma.end(), lower_than);

//   // move the node
//   E_Int Mi, nb_merges(0);
//   size_t psz = palma.size();
//   for (size_t i = 0; i < psz; ++i)
//   {
//     Fi = palma[i].Ni;
//     Mi = palma[i].Nj;

//     // rule of merging : if target and moving has not been moved already
//     if ((new_IDs[Fi] == Fi) && (new_IDs[Mi] == Mi))
//     {
//       new_IDs[Mi] = Fi;
//       ++nb_merges;
//     }
//   }
  
//   // update the pointers to point to the leaves
//   for (size_t i = 0; i < nsz; ++i)
//   {
//     Fi = new_IDs[i];
//     while (Fi != new_IDs[Fi]) Fi = new_IDs[Fi];
//     new_IDs[i] = Fi;
//   }

//   return nb_merges;
// }

template <typename ArrayType>
E_Int
__merge
(const K_FLD::ArrayAccessor<ArrayType>& coordAcc,
 const std::vector<E_Float>& nodal_metric2, 
 E_Float RTOL,
 NUGA::int_vector_type& umoving,
 NUGA::int_vector_type& utarget,
 NUGA::int_vector_type& new_IDs)
{
  size_t ssz = umoving.size();
  size_t tsz = utarget.size();
  size_t nsz = coordAcc.size();

  new_IDs.resize(nsz);
  for (size_t i = 0; i < nsz; ++i) new_IDs[i] = i;

  // Fast return
  if (ssz*tsz*nsz == 0) return 0;

  // build the KdTree on moving nodes.
  K_SEARCH::KdTree<ArrayType> moving_tree(coordAcc, umoving);

  //E_Float d2, tol2(tol*tol);
  typedef std::vector<triplet_t> palmares_t;
  palmares_t palma;

  std::vector<E_Int> onodes;
  std::vector<E_Float> dist2;

  // loop on target nodes and get all moving nodes in sphere of radius tol.
  size_t osz;
  //long int count{0};
  for (size_t i = 0; i < tsz; ++i)
  {
    onodes.clear(); dist2.clear();

    E_Int& Fi = utarget[i];
    assert (Fi >= 0 && (size_t)Fi < nodal_metric2.size());
    if (nodal_metric2[Fi] == NUGA::FLOAT_MAX)
    {
      //std::cout << "wrong TOLi" << std::endl;
      continue;
    }
    
    double TOLi = ::sqrt(nodal_metric2[Fi]) * RTOL;

    moving_tree.getInSphere(Fi, TOLi, onodes, dist2);
    osz = onodes.size();
    for (size_t j = 0; j < osz; ++j)
    {
      E_Int& Fj = onodes[j];

      if (nodal_metric2[Fj] == NUGA::FLOAT_MAX)
      {
        //std::cout << "wrong TOLj" << std::endl;
        continue;
      }

      double TOLj2 = nodal_metric2[Fj] * RTOL * RTOL;

      if (TOLj2 < dist2[j])
      {
        //++count;
        continue;
      }

      palma.push_back(triplet_t(dist2[j], Fi, Fj));
    }
  }

  //std::cout << "NB OF REJECTED : " << count << std::endl;

  if (palma.empty()) return 0;

  // sort them by increasing distance (keep intial order in case of equality to priorize it).
  std::stable_sort(palma.begin(), palma.end(), lower_than);

  // move the node
  E_Int Fi, Fj, nb_merges(0);
  size_t psz = palma.size();
  for (size_t i = 0; i < psz; ++i)
  {
    Fi = palma[i].Ni;
    Fj = palma[i].Nj;

    // rule of merging : if target and moving has not been moved already
    if ((new_IDs[Fi] == Fi) && (new_IDs[Fj] == Fj))
    {
      new_IDs[Fj] = Fi;
      ++nb_merges;
    }
  }

  // update the pointers to point to the leaves
  for (size_t i = 0; i < nsz; ++i)
  {
    Fi = new_IDs[i];
    while (Fi != new_IDs[Fi]) Fi = new_IDs[Fi];
    new_IDs[i] = Fi;
  }

  return nb_merges;
}

template <typename ArrayType>
E_Int
__merge_omp
(const K_FLD::ArrayAccessor<ArrayType>& coordAcc,
 const std::vector<E_Float>& nodal_metric2, 
 E_Float RTOL,
 NUGA::int_vector_type& umoving,
 NUGA::int_vector_type& utarget,
 NUGA::int_vector_type& new_IDs)
{
  size_t ssz = umoving.size();
  size_t tsz = utarget.size();
  size_t nsz = coordAcc.size();

  new_IDs.resize(nsz);
  for (size_t i = 0; i < nsz; ++i) new_IDs[i] = i;

  // Fast return
  if (ssz*tsz*nsz == 0) return 0;

  // build the KdTree on moving nodes.
  K_SEARCH::KdTree<ArrayType> moving_tree(coordAcc, umoving, EPSILON, true /*do omp*/);

  //const int NB_CANDS=50;

  using palmares_t = std::vector<triplet_t>;

  std::vector<palmares_t> palma_thrd(__NUMTHREADS__);
  palmares_t palma;

  std::vector<std::vector<E_Int>> onodes_thrd(__NUMTHREADS__);
  std::vector<std::vector<E_Float>> dist2_thrd(__NUMTHREADS__);

  // loop on target nodes and get all moving nodes in sphere of radius tol.
//#pragma omp parallel shared(tsz, moving_tree, palma, utarget, palma_thrd, dist2_thrd, onodes_thrd, RTOL, nodal_metric2) private (i, j, id, osz, Fi, Fj, d2, TOLi, TOLj2) default(none)
  #pragma omp parallel
  {
    E_Int id, Fi, Fj;
    size_t i, j, osz;
    double d2, TOLi, TOLj2;

    id = __CURRENT_THREAD__;

#pragma omp for schedule(static)
    for (i = 0; i < tsz; ++i)
    {
      onodes_thrd[id].clear(); dist2_thrd[id].clear();

      Fi = utarget[i];
      assert (Fi >= 0 && (size_t)Fi < nodal_metric2.size());
      if (nodal_metric2[Fi] == NUGA::FLOAT_MAX)
      {
        //std::cout << "wrong TOLi" << std::endl;
        continue;
      }
    
      TOLi = ::sqrt(nodal_metric2[Fi]) * RTOL;

      moving_tree.getInSphere(Fi, TOLi, onodes_thrd[id], dist2_thrd[id]);
    
      osz = onodes_thrd[id].size();

      for (j = 0; j < osz; ++j)
      {
        d2 = dist2_thrd[id][j];
        Fj = onodes_thrd[id][j];

        if (nodal_metric2[Fj] == NUGA::FLOAT_MAX)
        {
          //std::cout << "wrong metric" << std::endl;
          continue;
        }

        TOLj2 = nodal_metric2[Fj] * RTOL * RTOL;

        if (TOLj2 < dist2_thrd[id][j])
        {
          //++count;
          continue;
        }

        palma_thrd[id].push_back(triplet_t(d2, Fi, Fj));
      }
    }
  }

  size_t sz = 0;
  for (E_Int t = 0; t < __NUMTHREADS__; ++t) sz += palma_thrd[t].size();

  palma.reserve(sz);

  for (E_Int t = 0; t < __NUMTHREADS__; ++t)
  {
    //std::cout << "palma_thrd[t] sz : " << palma_thrd[t].size() << std::endl;
    if (!palma_thrd[t].empty())
      palma.insert(palma.end(), ALL(palma_thrd[t]));
  }

  if (palma.empty()) return 0;

  // sort them by increasing distance (keep intial order in case of equality to priorize it).
  std::stable_sort(palma.begin(), palma.end(), lower_than);

  // move the node
  E_Int nb_merges(0);
  E_Int Fi, Fj;
  size_t psz = palma.size();
  for (size_t i = 0; i < psz; ++i)
  {
    Fi = palma[i].Ni;
    Fj = palma[i].Nj;

    // rule of merging : if target and moving has not been moved already
    if ((new_IDs[Fi] == Fi) && (new_IDs[Fj] == Fj))
    {
      new_IDs[Fj] = Fi;
      ++nb_merges;
    }
  }

  // update the pointers to point to the leaves
  for (size_t i = 0; i < nsz; ++i)
  {
    Fi = new_IDs[i];
    while (Fi != new_IDs[Fi]) Fi = new_IDs[Fi];
    new_IDs[i] = Fi;
  }

  return nb_merges;
}

// Merge moving (mergeable) cloud into target cloud.
// targets might be moving too.
// Ni  Nj    Nk         Nl
// *   *     *          *
//    >-------< tol
// Ni -> Nj, Nj -> Nk
//
// so need to post process moving nodes to have all the node
// pointing to the end of the chain
template <typename ArrayType>
E_Int
merge
(const K_FLD::ArrayAccessor<ArrayType>& coordAcc, E_Float tol,
 const NUGA::int_vector_type& moving,
 const NUGA::int_vector_type& target,
 NUGA::int_vector_type& new_IDs, bool do_omp=false)
{
  // Fast return
  if (moving.size()*target.size()*coordAcc.size() == 0) return 0;
  
  // get unique lists.
  NUGA::int_vector_type umoving, utarget;
  stable_unique_values(moving, umoving);
  stable_unique_values(target, utarget);

  if (!do_omp)
    return ::__merge(coordAcc, tol, umoving, utarget, new_IDs);
  else
    return ::__merge_omp(coordAcc, tol, umoving, utarget, new_IDs);
}

// Merge points in the input coordinate matrix between them.
// Rather than doing the merge by modifying the coordinates, the new ids are given upon exit.
//==============================================================================
template <typename ArrayType>
E_Int
merge
(const K_FLD::ArrayAccessor<ArrayType>& coordAcc, E_Float tol,
 NUGA::int_vector_type& new_IDs, bool do_omp=false)
{
  NUGA::int_vector_type nodes(coordAcc.size());
  for (size_t i = 0; i < nodes.size(); ++i) nodes[i]=i;

  if (!do_omp)
    return ::__merge(coordAcc, tol, nodes, nodes, new_IDs);
  else
    return ::__merge_omp(coordAcc, tol, nodes, nodes, new_IDs);
}

template <typename ArrayType>
E_Int
merge
(const K_FLD::ArrayAccessor<ArrayType>& coordAcc, const std::vector<E_Float>&nodal_metric2,  E_Float RTOL,
  NUGA::int_vector_type& new_IDs, bool do_omp=false)
{
  NUGA::int_vector_type nodes(coordAcc.size());
  for (size_t i = 0; i < nodes.size(); ++i)nodes[i] = i;

  if (!do_omp)
    return ::__merge(coordAcc, nodal_metric2, RTOL, nodes, nodes, new_IDs);
  else
    return ::__merge_omp(coordAcc, nodal_metric2, RTOL, nodes, nodes, new_IDs);
}

/////////////////////////// temporary here to avoid cycling dep with EltAlgo
/// for BAR, TRI, QUAD, NGON-PG
template <typename Connectivity_t>
void getNodeToNodes
(const K_FLD::ArrayAccessor<Connectivity_t>& ELTContainer, std::map<E_Int, NUGA::int_vector_type > & bound_to_bounds)
{
  E_Int* pN;
  E_Int  ROWS(ELTContainer.stride()), COLS(ELTContainer.size()), s;
  
  pN = new E_Int[ROWS];

  for (E_Int Si = 0; Si < COLS; ++Si)
  {
    ELTContainer.getEntry(Si, pN);
    s=ELTContainer.stride(Si);
    
    for (E_Int n = 0; n < s; ++n)
    {
      bound_to_bounds[*(pN+n)].push_back(*(pN+(n+1)%s));
      bound_to_bounds[*(pN+(n+1)%s)].push_back(*(pN+n));
    }
  }
  
  delete [] pN;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Merge points in a mesh connectAcc (BAR, TRI, QUAD, NGON-PG) when the edge length is lower than tol
// Merge occurs only between connected nodes to preserve the topology
// Merge moving (mergeable) nodes into target nodes.
// targets might be moving too.
// Ni  Nj    Nk         Nl
// *   *     *          *
//    >-------< tol
// Ni -> Nj, Nj -> Nk
//
// so need to post process moving nodes to have all the node
// pointing to the end of the chain
template <typename CoordType, typename ConnectType>
E_Int
merge
(const K_FLD::ArrayAccessor<CoordType>& coordAcc, const K_FLD::ArrayAccessor<ConnectType>& connectAcc, E_Float tol,
 const NUGA::int_vector_type& moving, const NUGA::int_vector_type& target,
 NUGA::int_vector_type& new_IDs)
{  
  // get unique lists.
  NUGA::int_vector_type umoving, utarget;
  stable_unique_values(moving, umoving);
  stable_unique_values(target, utarget);

  size_t ssz = umoving.size();
  size_t tsz = utarget.size();
  size_t nsz = coordAcc.size();

  // Fast return
  if (ssz*tsz*nsz == 0) return 0;

  new_IDs.resize(nsz);
  for (size_t i = 0; i < nsz; ++i) new_IDs[i]=i;
  
  // get the connected nodes graph
  //typedef K_MESH::Triangle dummy;
  //typedef typename NUGA::EltAlgo<dummy> algo_t;
  //algo_t::NodeToNodesType n_to_ns;
  //algo_t::getNodeToNodes(connectAcc, n_to_ns);
  std::map<E_Int, NUGA::int_vector_type > n_to_ns;
  getNodeToNodes(connectAcc, n_to_ns);
  
  typedef std::vector<triplet_t > palmares_t;
  palmares_t palma;
  
  // loop on target nodes and get all connected nodes closer than tol.
  size_t osz;
  E_Float d2, tol2(tol*tol);
  for (size_t i = 0; i < tsz; ++i)
  {
    const E_Int& Fi = utarget[i];
    osz = n_to_ns[Fi].size();
    for (size_t j = 0; j < osz; ++j)
    {
      const E_Int& Mj = n_to_ns[Fi][j];
      d2 = coordAcc.dist2(Fi, Mj);
      if (d2 < tol2)
        palma.push_back(triplet_t(d2, Fi, Mj));
    }
  }

  // sort them by increasing distance (keep intial order in case of equality to priorize it).
  std::stable_sort(palma.begin(), palma.end(), lower_than);

  // move the node
  E_Int Fi, Mi, nb_merges(0);
  size_t psz = palma.size();
  for (size_t i =0; i < psz; ++i)
  {
    Fi=palma[i].Ni;
    Mi=palma[i].Nj;

    // rule of merging : if target and moving has not been moved already
    if ((new_IDs[Fi] == Fi) && (new_IDs[Mi] == Mi))
    {
      new_IDs[Mi]=Fi;
      ++nb_merges;
    }
  }

  // update the pointers to point to the leaves
  for (size_t i =0; i < nsz; ++i)
  {
    Fi=new_IDs[i];
    while (Fi != new_IDs[Fi])Fi=new_IDs[Fi];
    new_IDs[i]=Fi;
  }

  return nb_merges;
}

//}

#endif
