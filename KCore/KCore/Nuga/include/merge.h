/*



--------- NUGA v1.0



*/
//Authors : S�m Landier (sam.landier@onera.fr)

#ifndef __MERGE_H__
#define __MERGE_H__

#include "Nuga/include/defs.h"
#include "Nuga/include/ArrayAccessor.h"
#include "Nuga/include/KdTree.h"
#include "Nuga/include/Triangle.h"

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

inline bool lower_than(const std::pair<E_Float, NUGA::int_pair_type >& p1,
                      const std::pair<E_Float, NUGA::int_pair_type >& p2)
{return p1.first < p2.first;}

//==============================================================================
template <typename ArrayType>
void 
merge_no_order_omp
(const K_FLD::ArrayAccessor<ArrayType>& coordAcc, E_Float tol,
 NUGA::int_vector_type& new_IDs)
{
  K_SEARCH::KdTree<ArrayType> tree(coordAcc);
  size_t npts = coordAcc.size();
  new_IDs.resize(npts);
  printf("merge no order\n");

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
    std::vector< std::vector<E_Int>* > allnodes(npts);
#pragma omp parallel default(shared)
    {
      E_Float Xn[3]; 
     
#pragma omp for schedule(dynamic)
      for (E_Int i = 0; i < npts; i++)
      {
        std::vector<E_Int>* onodes = new std::vector<E_Int>;
        coordAcc.getEntry(i, Xn);
        //printf("pt %f %f %f\n", Xn[0], Xn[1], Xn[2]);
        tree.getInSphere(Xn, tol, *onodes);
        allnodes[i] = onodes;
      }
    }

    E_Int m; E_Float d2;
    E_Float Xn[3]; E_Float Xm[3];
    E_Int mID, iID;  size_t size; 
    
    for (E_Int i = 0; i < npts; i++)
    {
      std::vector<E_Int>& onodes = *allnodes[i]; 
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
    
    // free
    for (E_Int i = 0; i < npts; i++) delete allnodes[i];
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
  typedef std::vector<std::pair <E_Float, std::pair <E_Int, E_Int> > > palmares_t;
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
      palma.push_back(std::make_pair(dist2[j], std::make_pair(Fi, onodes[j]))); // (d2, Fi, Mi)
  }

  if (palma.empty()) return 0;

  // sort them by increasing distance (keep intial order in case of equality to priorize it).
  std::stable_sort(palma.begin(), palma.end(), lower_than);

  // move the node
  E_Int Fi, Mi, nb_merges(0);
  size_t psz = palma.size();
  for (size_t i = 0; i < psz; ++i)
  {
    Fi=palma[i].second.first;
    Mi=palma[i].second.second;

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
  E_Int ssz = (E_Int)umoving.size();
  E_Int tsz = (E_Int)utarget.size();
  E_Int nsz = (E_Int)coordAcc.size();

  new_IDs.resize(nsz);
  for (E_Int i = 0; i < nsz; ++i) new_IDs[i] = i;

  // Fast return
  if (ssz*tsz*nsz == 0) return 0;

  // build the KdTree on moving nodes.
  K_SEARCH::KdTree<ArrayType> moving_tree(coordAcc, umoving);

  //E_Float d2, tol2(tol*tol);
  typedef std::vector<std::pair <E_Float, std::pair <E_Int, E_Int> > > palmares_t;
  palmares_t palma;

  std::vector<E_Int> onodes;
  std::vector<E_Float> dist2;

  // loop on target nodes and get all moving nodes in sphere of radius tol.
  size_t osz;
  //long int count{0};
  for (E_Int i = 0; i < tsz; ++i)
  {
    onodes.clear(); dist2.clear();

    E_Int& Fi = utarget[i];
    assert (Fi >= 0 && Fi < nodal_metric2.size());
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

      palma.push_back(std::make_pair(dist2[j], std::make_pair(Fi, Fj))); // (d2, Fi, Fj)
    }
  }

  //std::cout << "NB OF REJECTED : " << count << std::endl;

  if (palma.empty()) return 0;

  // sort them by increasing distance (keep intial order in case of equality to priorize it).
  std::stable_sort(palma.begin(), palma.end(), lower_than);

  // move the node
  E_Int Fi, Fj, nb_merges(0);
  E_Int psz = palma.size();
  for (E_Int i = 0; i < psz; ++i)
  {
    Fi = palma[i].second.first;
    Fj = palma[i].second.second;

    // rule of merging : if target and moving has not been moved already
    if ((new_IDs[Fi] == Fi) && (new_IDs[Fj] == Fj))
    {
      new_IDs[Fj] = Fi;
      ++nb_merges;
    }
  }

  // update the pointers to point to the leaves
  for (E_Int i = 0; i < nsz; ++i)
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
 NUGA::int_vector_type& new_IDs)
{
  // Fast return
  if (moving.size()*target.size()*coordAcc.size() == 0) return 0;
  
  // get unique lists.
  NUGA::int_vector_type umoving, utarget;
  stable_unique_values(moving, umoving);
  stable_unique_values(target, utarget);

  return ::__merge(coordAcc, tol, umoving, utarget, new_IDs);
}

// Merge points in the input coordinate matrix between them.
// Rather than doing the merge by modifying the coordinates, the new ids are given upon exit.
//==============================================================================
template <typename ArrayType>
E_Int
merge
(const K_FLD::ArrayAccessor<ArrayType>& coordAcc, E_Float tol,
 NUGA::int_vector_type& new_IDs)
{
  NUGA::int_vector_type nodes(coordAcc.size());
  for (size_t i = 0; i < nodes.size(); ++i)nodes[i]=i;

  return ::__merge(coordAcc, tol, nodes, nodes, new_IDs);
}

template <typename ArrayType>
E_Int
merge
(const K_FLD::ArrayAccessor<ArrayType>& coordAcc, const std::vector<E_Float>&nodal_metric2,  E_Float RTOL,
  NUGA::int_vector_type& new_IDs)
{
  NUGA::int_vector_type nodes(coordAcc.size());
  for (size_t i = 0; i < nodes.size(); ++i)nodes[i] = i;

  return ::__merge(coordAcc, nodal_metric2, RTOL, nodes, nodes, new_IDs);
}

/////////////////////////// temporary here to avoid cycling dep with EltAlgo
/// for BAR, TRI, QUAD, NGON-PG
template <typename Connectivity_t>
void getNodeToNodes
(const K_FLD::ArrayAccessor<Connectivity_t>& ELTContainer, std::map<E_Int, NUGA::int_vector_type > & bound_to_bounds)
{
  E_Int*                          pN;
  E_Int                           ROWS(ELTContainer.stride()), COLS(ELTContainer.size()), s;
  
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
  
  typedef std::vector<std::pair <E_Float, std::pair <E_Int, E_Int> > > palmares_t;
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
        palma.push_back(std::make_pair(d2, std::make_pair(Fi, Mj))); // (d2, Fi, Mj)
    }
  }

  // sort them by increasing distance (keep intial order in case of equality to priorize it).
  std::stable_sort(palma.begin(), palma.end(), lower_than);

  // move the node
  E_Int Fi, Mi, nb_merges(0);
  size_t psz = palma.size();
  for (size_t i =0; i < psz; ++i)
  {
    Fi=palma[i].second.first;
    Mi=palma[i].second.second;

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
